"""veggenes data core — cBioPortal + Ensembl access shared by the web app.

Consolidates the verified logic from the two notebooks:
  * length-normalized mutation frequency  (mutated-genes/fetch + Ensembl canonical CDS)
  * co-occurrence / mutual exclusivity     (mutation/CNA/SV alteration sets + Fisher's exact)
"""
import time
import threading
import requests
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

CBIO = "https://www.cbioportal.org/api"
ENSEMBL = "https://rest.ensembl.org"
JSON_HDR = {"Accept": "application/json", "Content-Type": "application/json"}

_session = requests.Session()
_lock = threading.Lock()
_cache = {}  # tiny in-process memo cache


def _get(url, **kw):
    r = _session.get(url, headers={"Accept": "application/json"}, timeout=180, **kw)
    r.raise_for_status()
    return r.json()


def _post(url, payload, params=None):
    r = _session.post(url, headers=JSON_HDR, params=params or {}, json=payload, timeout=300)
    r.raise_for_status()
    return r.json()


# ---------------------------------------------------------------- studies / cohorts
def find_studies(keyword="", limit=40):
    studies = _cache.get("studies")
    if studies is None:
        studies = _get(f"{CBIO}/studies", params={"projection": "SUMMARY"})
        _cache["studies"] = studies
    k = (keyword or "").lower()
    rows = [s for s in studies
            if not k or k in s.get("name", "").lower() or k in s.get("studyId", "").lower()]
    rows.sort(key=lambda s: -(s.get("allSampleCount") or 0))
    return [{"studyId": s["studyId"], "name": s.get("name", ""),
             "allSampleCount": s.get("allSampleCount")} for s in rows[:limit]]


def _clinical_records(study_id, attribute_id="CANCER_TYPE_DETAILED"):
    key = ("clin", study_id, attribute_id)
    with _lock:
        if key in _cache:
            return _cache[key]
    out, page, size = [], 0, 50_000
    while True:
        cd = _get(f"{CBIO}/studies/{study_id}/clinical-data",
                  params={"attributeId": attribute_id, "clinicalDataType": "SAMPLE",
                          "projection": "SUMMARY", "pageSize": size, "pageNumber": page})
        if isinstance(cd, dict):
            raise RuntimeError(f"clinical-data error: {cd}")
        out += [(c["sampleId"], c.get("patientId"), c.get("value")) for c in cd]
        if len(cd) < size:
            break
        page += 1
    df = pd.DataFrame(out, columns=["sampleId", "patientId", "value"])
    with _lock:
        _cache[key] = df
    return df


def list_cancer_types(study_id, top=60):
    df = _clinical_records(study_id)
    s = df["value"].dropna().value_counts().head(top)
    return [{"value": v, "count": int(n)} for v, n in s.items()]


def get_cohort(study_id, cancer_type=None):
    df = _clinical_records(study_id)
    if cancer_type:
        df = df[df["value"] == cancer_type]
        if df.empty:
            raise ValueError(f"No samples with CANCER_TYPE_DETAILED == {cancer_type!r}")
    return df[["sampleId", "patientId"]].drop_duplicates().reset_index(drop=True)


def _profiles(study_id):
    key = ("prof", study_id)
    if key in _cache:
        return _cache[key]
    out = {}
    for p in _get(f"{CBIO}/studies/{study_id}/molecular-profiles"):
        t = p["molecularAlterationType"]
        if t == "MUTATION_EXTENDED":
            out["mutation"] = p["molecularProfileId"]
        elif t == "COPY_NUMBER_ALTERATION" and p.get("datatype") == "DISCRETE":
            out["cna"] = p["molecularProfileId"]
        elif t == "STRUCTURAL_VARIANT":
            out["sv"] = p["molecularProfileId"]
    _cache[key] = out
    return out


# ---------------------------------------------------------------- frequency tab
def _mutated_genes(study_id, cancer_type=None):
    body = {"studyIds": [study_id]}
    if cancer_type:
        body["clinicalDataFilters"] = [{"attributeId": "CANCER_TYPE_DETAILED",
                                        "values": [{"value": cancer_type}]}]
    data = _post(f"{CBIO}/mutated-genes/fetch", body)
    return data


def _ensembl_cds_lengths(genes, seq_type="cds"):
    """Hugo symbol -> canonical-transcript CDS length (bp), batched, cached."""
    genes = list(genes)
    out, todo = {}, []
    for g in genes:
        ck = ("cds", g)
        if ck in _cache:
            out[g] = _cache[ck]
        else:
            todo.append(g)

    def chunks(seq, n):
        for i in range(0, len(seq), n):
            yield seq[i:i + n]

    canonical = {}
    for batch in chunks(todo, 50):
        try:
            res = _post(f"{ENSEMBL}/lookup/symbol/homo_sapiens", {"symbols": batch})
        except Exception:
            continue
        for sym, info in res.items():
            ct = info.get("canonical_transcript")
            if ct:
                canonical[sym] = ct.split(".")[0]
        time.sleep(0.15)

    tx_to_genes = {}
    for sym, tx in canonical.items():
        tx_to_genes.setdefault(tx, []).append(sym)
    for batch in chunks(list(tx_to_genes), 50):
        try:
            res = _post(f"{ENSEMBL}/sequence/id", {"ids": batch}, params={"type": seq_type})
        except Exception:
            continue
        for rec in res:
            for sym in tx_to_genes[rec["query"]]:
                out[sym] = len(rec["seq"])
                _cache[("cds", sym)] = len(rec["seq"])
        time.sleep(0.15)
    return out


def frequencies(study_id, cancer_type=None, top_n=25):
    data = _mutated_genes(study_id, cancer_type)
    if not data:
        return {"n": 0, "genes": []}
    data.sort(key=lambda g: -g.get("numberOfAlteredCases", 0))
    top = data[:top_n]
    n = top[0]["numberOfProfiledCases"]
    genes = [g["hugoGeneSymbol"] for g in top]
    cds = _ensembl_cds_lengths(genes)

    out = []
    for g in top:
        sym = g["hugoGeneSymbol"]
        prof = g["numberOfProfiledCases"] or n
        freq = 100 * g["numberOfAlteredCases"] / prof
        bp = cds.get(sym)
        out.append({
            "gene": sym,
            "altered": g["numberOfAlteredCases"],
            "total_mut": g.get("totalCount"),
            "profiled": prof,
            "freq_pct": round(freq, 2),
            "cds_bp": bp,
            "freq_per_kb": round(freq / (bp / 1000.0), 3) if bp else None,
        })
    return {"n": int(n), "genes": out}


# ---------------------------------------------------------------- co-occurrence tab
def resolve_genes(symbols):
    res = _post(f"{CBIO}/genes/fetch", [s.upper() for s in symbols],
                params={"geneIdType": "HUGO_GENE_SYMBOL"})
    found = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in res}
    return found


def _altered_sample_sets(study_id, entrez_map, sample_ids, alterations):
    e2g = {v: k for k, v in entrez_map.items()}
    sset = set(sample_ids)
    sets = {g: set() for g in entrez_map}
    prof = _profiles(study_id)
    ids = list(entrez_map.values())

    muts = _post(f"{CBIO}/molecular-profiles/{prof['mutation']}/mutations/fetch",
                 {"entrezGeneIds": ids, "sampleIds": sample_ids}, params={"projection": "SUMMARY"})
    for m in muts:
        if m["entrezGeneId"] in e2g and m["sampleId"] in sset:
            sets[e2g[m["entrezGeneId"]]].add(m["sampleId"])

    if alterations == "altered":
        if "cna" in prof:
            cna = _post(f"{CBIO}/molecular-profiles/{prof['cna']}/molecular-data/fetch",
                        {"entrezGeneIds": ids, "sampleIds": sample_ids})
            for d in cna:
                if d.get("value") in (-2, 2) and d["entrezGeneId"] in e2g and d["sampleId"] in sset:
                    sets[e2g[d["entrezGeneId"]]].add(d["sampleId"])
        if "sv" in prof:
            sv = _post(f"{CBIO}/structural-variant/fetch",
                       {"molecularProfileIds": [prof["sv"]], "entrezGeneIds": ids})
            for x in sv:
                sid = x.get("sampleId")
                if sid in sset:
                    for sym in (x.get("site1HugoSymbol"), x.get("site2HugoSymbol")):
                        if sym in sets:
                            sets[sym].add(sid)
    return sets


def _binary_matrix(study_id, entrez_map, cohort, level, alterations):
    sample_ids = cohort["sampleId"].tolist()
    sets = _altered_sample_sets(study_id, entrez_map, sample_ids, alterations)
    s2p = dict(zip(cohort["sampleId"], cohort["patientId"]))
    genes = list(entrez_map)
    if level == "patient":
        index_ids = cohort["patientId"].drop_duplicates().tolist()
        to_unit = lambda s: s2p.get(s)
    else:
        index_ids = sample_ids
        to_unit = lambda s: s
    bm = pd.DataFrame(0, index=index_ids, columns=genes)
    for g in genes:
        units = {to_unit(s) for s in sets[g]} - {None}
        if units:
            bm.loc[list(units), g] = 1
    return bm


def cooccurrence(study_id, genes, cancer_type=None, level="sample", alterations="altered"):
    genes = [g.strip().upper() for g in genes if g.strip()]
    entrez_map = resolve_genes(genes)
    order = [g for g in genes if g in entrez_map]
    dropped = [g for g in genes if g not in entrez_map]
    if len(order) < 2:
        return {"error": "Need at least 2 valid genes.", "dropped": dropped}

    cohort = get_cohort(study_id, cancer_type)
    entrez_map = {g: entrez_map[g] for g in order}
    bm = _binary_matrix(study_id, entrez_map, cohort, level, alterations)
    n = bm.shape[0]

    cells = {}  # (a,b) -> [both, a_not_b, b_not_a, neither]
    count = [[0] * len(order) for _ in order]
    log2or = [[None] * len(order) for _ in order]
    pairs = []
    for i, gx in enumerate(order):
        for j, gy in enumerate(order):
            x, y = bm[gx], bm[gy]
            a = int(((x == 1) & (y == 1)).sum())
            count[j][i] = a
            if i == j:
                continue
            b = int(((x == 1) & (y == 0)).sum())
            c = int(((x == 0) & (y == 1)).sum())
            d = int(((x == 0) & (y == 0)).sum())
            _, p = fisher_exact([[a, b], [c, d]])
            if min(a, b, c, d) == 0:
                lor = np.log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
            else:
                lor = np.log2((a * d) / (b * c))
            log2or[j][i] = round(float(lor), 3)
            pairs.append((frozenset((gx, gy)), p, (gx, gy, a, b, c, d, float(lor))))

    seen, uniq = {}, []
    for key, p, info in pairs:
        if key not in seen:
            seen[key] = p
            uniq.append((key, p, info))
    qmap = {}
    if uniq:
        qs = multipletests([p for _, p, _ in uniq], method="fdr_bh")[1]
        qmap = {k: q for (k, _, _), q in zip(uniq, qs)}

    pair_rows = []
    for key, p, (gx, gy, a, b, c, d, lor) in uniq:
        q = float(qmap.get(key, np.nan))
        pair_rows.append({"A": gx, "B": gy, "both": a, "a_not_b": b, "b_not_a": c, "neither": d,
                          "log2_or": round(lor, 3), "p": p, "q": q,
                          "tendency": "co-occurrence" if lor > 0 else "mutual exclusivity",
                          "sig": bool(q < 0.05)})
    pair_rows.sort(key=lambda r: r["q"])

    # significance matrix (q<0.05) aligned with log2or
    qmat = [[None] * len(order) for _ in order]
    for i, gx in enumerate(order):
        for j, gy in enumerate(order):
            if i != j:
                qmat[j][i] = float(qmap.get(frozenset((gx, gy)), np.nan))

    return {"genes": order, "dropped": dropped, "n": int(n), "level": level,
            "alterations": alterations, "count": count, "log2or": log2or, "qmat": qmat,
            "pairs": pair_rows}
