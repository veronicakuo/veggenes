"use strict";

const state = {
  study: null,
  studyName: null,
  cancerType: null,
  freq: null,          // last frequencies payload
  basket: [],          // genes for co-occurrence
  metric: "log2or",    // or "count"
  cooc: null,          // last cooccurrence payload
};

const $ = (id) => document.getElementById(id);
const spinner = $("spinner");
let busy = 0;
function setBusy(on) { busy += on ? 1 : -1; spinner.hidden = busy <= 0; }

async function api(url, opts) {
  setBusy(true);
  try {
    const r = await fetch(url, opts);
    const j = await r.json();
    if (j && j.error) throw new Error(j.error);
    return j;
  } finally { setBusy(false); }
}
function debounce(fn, ms) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }

/* ---------------- study search ---------------- */
const studyInput = $("studySearch"), studyResults = $("studyResults");

const searchStudies = debounce(async () => {
  const q = studyInput.value.trim();
  const rows = await api(`/api/studies?q=${encodeURIComponent(q)}&limit=25`);
  studyResults.innerHTML = "";
  rows.forEach((s) => {
    const li = document.createElement("li");
    li.innerHTML = `${s.name} <span class="sid">${s.studyId} · n=${s.allSampleCount}</span>`;
    li.onclick = () => selectStudy(s);
    studyResults.appendChild(li);
  });
  studyResults.hidden = rows.length === 0;
}, 220);

studyInput.addEventListener("input", searchStudies);
studyInput.addEventListener("focus", searchStudies);
document.addEventListener("click", (e) => {
  if (!studyInput.contains(e.target) && !studyResults.contains(e.target)) studyResults.hidden = true;
});

async function selectStudy(s) {
  state.study = s.studyId; state.studyName = s.name; state.cancerType = null;
  studyInput.value = s.name;
  studyResults.hidden = true;
  const sel = $("cancerType");
  sel.innerHTML = `<option value="">Whole study</option>`;
  sel.disabled = true;
  $("cohortInfo").innerHTML = `Loading cancer types for <b>${s.studyId}</b>…`;
  const cts = await api(`/api/cancer_types?study=${encodeURIComponent(s.studyId)}`);
  cts.forEach((c) => {
    const o = document.createElement("option");
    o.value = c.value; o.textContent = `${c.value} (n=${c.count})`;
    sel.appendChild(o);
  });
  sel.disabled = false;
  $("cohortInfo").innerHTML = `<b>${s.name}</b> — pick a cancer type or analyze the whole study.`;
  loadFrequencies();
}

$("cancerType").addEventListener("change", (e) => {
  state.cancerType = e.target.value || null;
  loadFrequencies();
});
$("topN").addEventListener("change", loadFrequencies);
$("normToggle").addEventListener("change", () => state.freq && renderFreq());

/* ---------------- frequency tab ---------------- */
async function loadFrequencies() {
  if (!state.study) return;
  $("freqEmpty").style.display = "none";
  const topN = $("topN").value;
  const url = `/api/frequencies?study=${encodeURIComponent(state.study)}` +
    (state.cancerType ? `&cancer_type=${encodeURIComponent(state.cancerType)}` : "") + `&top_n=${topN}`;
  state.freq = await api(url);
  const label = state.cancerType || "whole study";
  $("cohortInfo").innerHTML =
    `<b>${state.studyName}</b> · <b>${label}</b> — ${state.freq.n.toLocaleString()} profiled samples.`;
  renderFreq();
}

function renderFreq() {
  const norm = $("normToggle").checked;
  const key = norm ? "freq_per_kb" : "freq_pct";
  let genes = state.freq.genes.filter((g) => g[key] != null);
  genes = genes.slice().sort((a, b) => a[key] - b[key]);   // ascending -> largest on top in plotly
  const y = genes.map((g) => g.gene);
  const x = genes.map((g) => g[key]);
  const inBasket = genes.map((g) => state.basket.includes(g.gene));
  const colors = inBasket.map((b) => (b ? "#2f6f4f" : "#4c72b0"));
  const custom = genes.map((g) => [g.freq_pct, g.freq_per_kb, g.cds_bp, g.altered]);

  const trace = {
    type: "bar", orientation: "h", x, y,
    marker: { color: colors },
    customdata: custom,
    hovertemplate:
      "<b>%{y}</b><br>frequency: %{customdata[0]}%<br>" +
      "per kb CDS: %{customdata[1]}<br>CDS: %{customdata[2]} bp<br>" +
      "altered samples: %{customdata[3]}<extra></extra>",
  };
  const layout = {
    margin: { l: 70, r: 20, t: 16, b: 44 },
    height: Math.max(260, genes.length * 26 + 60),
    xaxis: { title: norm ? "% samples mutated per kb CDS" : "% samples mutated", zeroline: false, gridcolor: "#eee" },
    yaxis: { automargin: true },
    plot_bgcolor: "#fff", paper_bgcolor: "#fff",
    font: { family: "inherit", color: "#1a1a1a" },
  };
  Plotly.react("freqPlot", [trace], layout, { displayModeBar: false, responsive: true });
  $("freqPlot").on("plotly_click", (ev) => {
    const gene = ev.points[0].y;
    addGene(gene);
    flashTab("cooc");
  });
}

/* ---------------- gene basket ---------------- */
function addGene(g) {
  g = g.trim().toUpperCase();
  if (g && !state.basket.includes(g)) { state.basket.push(g); renderBasket(); if (state.freq) renderFreq(); }
}
function removeGene(g) {
  state.basket = state.basket.filter((x) => x !== g); renderBasket(); if (state.freq) renderFreq();
}
function renderBasket() {
  const box = $("geneBasket");
  box.innerHTML = "";
  state.basket.forEach((g) => {
    const chip = document.createElement("span");
    chip.className = "chip";
    chip.innerHTML = `${g} <button title="remove">×</button>`;
    chip.querySelector("button").onclick = () => removeGene(g);
    box.appendChild(chip);
  });
}

const geneInput = $("geneSearch"), geneResults = $("geneResults");
geneInput.addEventListener("input", () => {
  const q = geneInput.value.trim().toUpperCase();
  geneResults.innerHTML = "";
  if (!q || !state.freq) { geneResults.hidden = true; return; }
  const sugg = state.freq.genes.map((g) => g.gene).filter((g) => g.startsWith(q) && !state.basket.includes(g)).slice(0, 8);
  sugg.forEach((g) => {
    const li = document.createElement("li");
    li.textContent = g;
    li.onclick = () => { addGene(g); geneInput.value = ""; geneResults.hidden = true; };
    geneResults.appendChild(li);
  });
  geneResults.hidden = sugg.length === 0;
});
geneInput.addEventListener("keydown", (e) => {
  if (e.key === "Enter") { e.preventDefault(); addGene(geneInput.value); geneInput.value = ""; geneResults.hidden = true; }
});
document.addEventListener("click", (e) => {
  if (!geneInput.contains(e.target) && !geneResults.contains(e.target)) geneResults.hidden = true;
});

/* ---------------- tabs ---------------- */
document.querySelectorAll(".tab").forEach((t) => {
  t.onclick = () => flashTab(t.dataset.tab);
});
function flashTab(name) {
  document.querySelectorAll(".tab").forEach((t) => t.classList.toggle("active", t.dataset.tab === name));
  document.querySelectorAll(".tabpane").forEach((p) => p.classList.toggle("active", p.id === `tab-${name}`));
  if (name === "cooc" && state.cooc) Plotly.Plots.resize("coocPlot");
  if (name === "freq" && state.freq) Plotly.Plots.resize("freqPlot");
}

/* ---------------- metric toggle ---------------- */
document.querySelectorAll("#metricToggle .seg").forEach((b) => {
  b.onclick = () => {
    document.querySelectorAll("#metricToggle .seg").forEach((x) => x.classList.toggle("active", x === b));
    state.metric = b.dataset.metric;
    if (state.cooc) renderCooc();
  };
});

/* ---------------- co-occurrence tab ---------------- */
$("runCooc").addEventListener("click", runCooccurrence);

async function runCooccurrence() {
  if (!state.study) { alert("Pick a study first."); return; }
  if (state.basket.length < 2) { alert("Add at least two genes."); return; }
  $("coocEmpty").style.display = "none";
  const payload = {
    study: state.study, cancer_type: state.cancerType, genes: state.basket,
    level: $("levelSel").value, alterations: $("altSel").value,
  };
  const data = await api("/api/cooccurrence", {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload),
  });
  if (data.error) { alert(data.error); return; }
  state.cooc = data;
  if (data.dropped && data.dropped.length) {
    state.basket = state.basket.filter((g) => !data.dropped.includes(g)); renderBasket();
  }
  renderCooc();
}

function renderCooc() {
  const d = state.cooc, genes = d.genes;
  let trace, layout;
  const label = state.cancerType || "whole study";
  const title = `${state.studyName} · ${label} — ${d.alterations}, ${d.level}-level (n=${d.n.toLocaleString()})`;

  if (state.metric === "count") {
    trace = {
      type: "heatmap", x: genes, y: genes, z: d.count,
      colorscale: "Blues", colorbar: { title: `# ${d.level}s` },
      hovertemplate: "%{y} & %{x}<br>both = %{z}<extra></extra>",
    };
    layout = baseHeatLayout(title);
    layout.annotations = textAnnotations(genes, d.count, (v) => v);
  } else {
    const z = d.log2or;
    trace = {
      type: "heatmap", x: genes, y: genes, z,
      zmid: 0, colorscale: "RdBu", reversescale: true,
      colorbar: { title: "log₂ OR" },
      customdata: d.qmat, hoverongaps: false,
      hovertemplate: "%{y} & %{x}<br>log₂ OR = %{z:.2f}<br>q = %{customdata:.3g}<extra></extra>",
    };
    layout = baseHeatLayout(title);
    // ★ on significant cells
    const ann = [];
    for (let j = 0; j < genes.length; j++)
      for (let i = 0; i < genes.length; i++) {
        const q = d.qmat[j][i];
        if (q != null && q < 0.05) ann.push(cellAnn(genes[i], genes[j], "★", "#111", 16));
      }
    layout.annotations = ann;
  }
  Plotly.react("coocPlot", [trace], layout, { displayModeBar: false, responsive: true });
}

function baseHeatLayout(title) {
  return {
    title: { text: title, font: { size: 13 }, x: 0, xanchor: "left" },
    margin: { l: 70, r: 20, t: 40, b: 60 },
    height: 460,
    xaxis: { side: "bottom", automargin: true, constrain: "domain" },
    yaxis: { autorange: "reversed", automargin: true, scaleanchor: "x", constrain: "domain" },
    plot_bgcolor: "#fff", paper_bgcolor: "#fff", font: { family: "inherit", color: "#1a1a1a" },
  };
}
function cellAnn(x, y, text, color, size) {
  return { x, y, text: String(text), showarrow: false, font: { color, size } };
}
function textAnnotations(genes, mat, fmt) {
  const ann = [];
  for (let j = 0; j < genes.length; j++)
    for (let i = 0; i < genes.length; i++)
      ann.push(cellAnn(genes[i], genes[j], fmt(mat[j][i]), "#333", 11));
  return ann;
}
