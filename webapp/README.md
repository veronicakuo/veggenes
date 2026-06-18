# veggenes — cancer-genomics explorer

A simple two-tab web app over the cBioPortal + Ensembl REST APIs.

* **Tab 1 — Mutation frequency:** length-normalized (canonical CDS) mutation frequency for the
  top mutated genes in the chosen cohort, as a horizontal bar plot. Click a bar to add that gene
  to the co-occurrence basket.
* **Tab 2 — Gene co-occurrence:** co-occurrence / mutual-exclusivity matrix for a chosen gene set.
  Toggle between **log₂ odds ratio** (Fisher's exact, ★ = FDR<0.05) and **co-occurrence count**.
  Pick sample/patient level and mutation-only / altered (mutation+CNA+fusion).

Defaults reproduce cBioPortal's Mutual Exclusivity tab (sample-level, altered, plain odds ratio).

## Run

```bash
conda activate veggenes        # has flask, pandas, scipy, statsmodels, requests
cd veggenes/webapp
python app.py --port 8050      # then open http://127.0.0.1:8050
```

Or use the helper: `./run.sh`

## Typical workflow

1. Search and select a **study** (e.g. *MSK-IMPACT 50K*), then pick a **cancer type** (e.g. *Lung Adenocarcinoma*).
2. The **Mutation frequency** tab loads the length-normalized top genes.
3. Click genes of interest (or type them in the co-occurrence search) to build a gene set.
4. The **Gene co-occurrence** tab shows the matrix; toggle metric / level / alteration type and recompute.

## Files
- `app.py` — Flask server + JSON API (`/api/studies`, `/api/cancer_types`, `/api/frequencies`, `/api/cooccurrence`)
- `veggenes_core.py` — data layer (consolidated from the two notebooks; verified against cBioPortal)
- `static/` — `index.html`, `style.css`, `app.js` (Plotly.js via CDN)
