#!/usr/bin/env python3
"""Build a self-contained, interactive cohort comparison report (docs/index.html).

Summarises a mtDNA-SV caller run into one offline HTML file: a caller x scenario
detection matrix, per-caller runtime, sensitivity, and call counts, with
interactive SVG figures (no external/CDN dependencies). Driven by the cohort
tables from postprocess.py plus the MitoHPC truth.

Usage:
  make_report.py --calls cohort_sv_calls.tsv --runtime cohort_runtime.tsv \\
                 --truth truth.tsv --samples "s1 s2 ..." --out docs/index.html \\
                 [--scope full] [--image ghcr.io/...] [--generated 2026-06-24]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from statistics import median

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import parsers  # noqa: E402
import sv_eval  # noqa: E402

GEN_TOL = 250
# Real-data expectations (not in the mock truth.tsv); committed real BAMs.
# 4-tuples (kind, bp5, bp3, expect) so they flow through sv_eval like mock truth.
REAL_TRUTH = {
    "spike_del4977_h20": [("del", 8469, 13447, "pass")],   # del4977 spiked @~20%
    "NA12718": [("none", None, None, "no_pass")],            # healthy 1000G
    "NA12748": [("none", None, None, "no_pass")],
    "NA12775": [("none", None, None, "no_pass")],
}


def _num(x):
    try:
        return int(round(float(x)))
    except (ValueError, TypeError):
        return None


def _f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def load_truth(path):
    """sample -> [(kind, bp5, bp3, expect)] from the MitoHPC truth.tsv
    (#sample kind bp5 bp3 svlen het depth expect)."""
    truth = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                expect = f[7] if len(f) > 7 else ""
                truth.setdefault(f[0], []).append((f[1], _num(f[2]), _num(f[3]), expect))
    truth.update(REAL_TRUTH)
    return truth


def load_scenarios(path):
    """scenarios.json -> (categories[list], descriptions{sample:popup})."""
    cats, desc = [], {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            d = json.load(fh)
        cats = d.get("categories", [])
        desc = d.get("descriptions", {})
    return cats, desc


def load_calls(path):
    calls = []
    if path and os.path.isfile(path):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["bp5"] = _num(row.get("bp5"))
                row["bp3"] = _num(row.get("bp3"))
                row["svlen"] = _num(row.get("svlen"))
                row["het"] = _f(row.get("het"))
                row["support"] = _f(row.get("support"))
                calls.append(row)
    return calls


def load_runtime(path):
    rt = []
    if path and os.path.isfile(path):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["seconds"] = _f(row.get("seconds"))
                rt.append(row)
    return rt


def _pctile(xs, q):
    if not xs:
        return None
    if len(xs) == 1:
        return xs[0]
    pos = (q / 100.0) * (len(xs) - 1)
    lo = int(pos); frac = pos - lo
    return xs[lo] * (1 - frac) + xs[min(lo + 1, len(xs) - 1)] * frac


def boxstats(vals, nd=1):
    """Five-number summary + Tukey whiskers/outliers + mean for a boxplot."""
    xs = sorted(v for v in vals if v is not None)
    n = len(xs)
    if n == 0:
        return None
    q1, med, q3 = _pctile(xs, 25), _pctile(xs, 50), _pctile(xs, 75)
    iqr = q3 - q1
    lof, hif = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inside = [x for x in xs if lof <= x <= hif]
    rnd = lambda v: round(v, nd)
    return {
        "n": n, "min": rnd(xs[0]), "q1": rnd(q1), "med": rnd(med),
        "q3": rnd(q3), "max": rnd(xs[-1]), "mean": rnd(sum(xs) / n),
        "wlo": rnd(min(inside) if inside else xs[0]),
        "whi": rnd(max(inside) if inside else xs[-1]),
        "outliers": [rnd(x) for x in xs if x < lof or x > hif][:40],
    }


def base_sample(s):
    return s[:-5] if s.endswith("_cram") else s


def _match(c5, c3, e5, e3, tol=GEN_TOL):
    if None in (c5, c3, e5, e3):
        return False
    return abs(c5 - e5) <= tol and abs(c3 - e3) <= tol


def build(calls, runtime, truth, samples, categories=None, descriptions=None):
    callers = list(parsers.CALLERS)
    calls_by_sample = {}
    for c in calls:
        calls_by_sample.setdefault(c["sample"], []).append(c)

    # ---- categorized scenario x caller trials + accuracy metrics (sv_eval) ----
    trials = sv_eval.build_trials(truth, samples, calls_by_sample, REAL_TRUTH)
    metrics = sv_eval.metrics_by_category(trials, callers)
    # matrix rows for the report (one per trial), with category + popup
    descriptions = descriptions or {}
    matrix = []
    for t in trials:
        matrix.append({
            "sample": t["sample"], "category": t["category"], "kind": t["kind"],
            "klass": t["klass"], "reason": t.get("reason", ""),
            "label": t["label"], "common": t["is_common"], "eval_only": t["eval_only"],
            "detected": t["detected"],
            "fp": {c: v for c, v in t["fp"].items()},
            "popup": descriptions.get(t["sample"], ""),
        })

    # ---- runtime per caller (ok runs only) ----
    rt_by_caller = {c: [] for c in callers}
    rt_per_sample = {c: {} for c in callers}
    status_by = {c: {} for c in callers}
    for r in runtime:
        c = r["caller"]
        if c not in rt_by_caller:
            rt_by_caller[c] = []; rt_per_sample[c] = {}; status_by[c] = {}
        status_by[c][r["sample"]] = r.get("status", "")
        if r["seconds"] is not None and r.get("status") == "ok":
            rt_by_caller[c].append(r["seconds"])
            rt_per_sample[c][r["sample"]] = r["seconds"]
    runtime_summary = {}
    for c in callers:
        secs = rt_by_caller.get(c, [])
        runtime_summary[c] = {
            "n": len(secs),
            "total": round(sum(secs), 1) if secs else 0.0,
            "mean": round(sum(secs) / len(secs), 1) if secs else None,
            "median": round(median(secs), 1) if secs else None,
            "box": boxstats(secs),
            "per_sample": rt_per_sample.get(c, {}),
        }

    # ---- per-caller "ran" flag (kept for the calls/scatter sections) ----
    det_summary = {}
    for c in callers:
        n_calls = sum(1 for x in calls if x["caller"] == c)
        det_summary[c] = {
            "n_calls": n_calls,
            "ran": runtime_summary[c]["n"] > 0 or n_calls > 0,
        }

    # categories present, ordered (Deletions first); report appends an "All" tab.
    cats_present = [cat for cat in sv_eval.CATEGORY_ORDER
                    if any(t["category"] == cat for t in trials)]
    cat_meta = {c["key"]: c for c in (categories or [])}
    cat_list = [{"key": k, "label": sv_eval.CATEGORY_LABEL.get(k, k),
                 "blurb": cat_meta.get(k, {}).get("blurb", ""),
                 "eval_only": k in sv_eval.EVAL_ONLY_CATS} for k in cats_present]

    return {
        "callers": callers,
        "samples": samples,
        "matrix": matrix,
        "categories": cat_list,
        "metrics": metrics,
        "runtime": runtime_summary,
        "detection": det_summary,
        "calls": [{k: c.get(k) for k in
                   ("sample", "caller", "sv_type", "bp5", "bp3", "svlen",
                    "support", "het", "common_deletion", "extra")} for c in calls],
    }


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
def render_html(data, meta):
    payload = json.dumps({"data": data, "meta": meta}).replace("</", "<\\/")
    return _HTML_TEMPLATE.replace("/*__PAYLOAD__*/", payload)


_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>mtDNA SV caller comparison</title>
<style>
:root{--bg:#0f1117;--panel:#171a23;--panel2:#1f2430;--ink:#e6e9ef;--mut:#9aa4b2;
--line:#2a2f3a;--accent:#6ea8fe;--good:#2fbf71;--bad:#e5534b;--warn:#e3a008;
--ref:#b692f6;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
a{color:var(--accent)}
.wrap{max-width:1100px;margin:0 auto;padding:28px 20px 80px}
h1{font-size:24px;margin:0 0 4px} h2{font-size:18px;margin:34px 0 12px;border-bottom:1px solid var(--line);padding-bottom:6px}
.sub{color:var(--mut);margin:0 0 18px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:18px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card .k{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em}
.card .v{font-size:22px;font-weight:600;margin-top:4px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:12px 0}
.legend{color:var(--mut);font-size:12px;margin:6px 0 2px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 8px;text-align:left;border-bottom:1px solid var(--line);white-space:nowrap}
th{color:var(--mut);font-weight:600;position:sticky;top:0;background:var(--panel)}
.matrix td.c{text-align:center;font-weight:700}
.cell-yes{background:rgba(47,191,113,.18);color:var(--good)}
.cell-no{color:#5b6472}
.cell-fp{background:rgba(229,83,75,.20);color:var(--bad)}
.tag{display:inline-block;padding:1px 7px;border-radius:999px;font-size:11px;border:1px solid var(--line)}
.tag.ref{color:var(--ref);border-color:var(--ref)}
.controls{margin:8px 0;display:flex;gap:8px;flex-wrap:wrap}
button.f{background:var(--panel2);color:var(--ink);border:1px solid var(--line);
border-radius:8px;padding:5px 10px;cursor:pointer;font-size:12px}
button.f.on{border-color:var(--accent);color:var(--accent)}
input[type=text]{background:var(--panel2);border:1px solid var(--line);color:var(--ink);
border-radius:8px;padding:5px 10px;font-size:12px;min-width:200px}
svg text{fill:var(--ink);font-size:11px}
svg .ax{stroke:var(--line)} svg .gl{stroke:var(--line);stroke-dasharray:2 3;opacity:.5}
.bar{cursor:pointer} .bar:hover{opacity:.85}
#tip{position:fixed;pointer-events:none;background:#000d;border:1px solid var(--line);
border-radius:8px;padding:6px 9px;font-size:12px;color:#fff;opacity:0;transition:opacity .1s;z-index:9}
.foot{color:var(--mut);font-size:12px;margin-top:30px;border-top:1px solid var(--line);padding-top:14px}
.scroll{overflow:auto;max-height:520px;border:1px solid var(--line);border-radius:10px}
.muted{color:var(--mut)}
.cell-amb{color:var(--mut);font-style:italic;font-weight:400}
.eval-banner{background:rgba(230,160,8,.10);border:1px solid #6b5410;color:#e3a008;
border-radius:8px;padding:6px 10px;font-size:12px;margin:0 0 8px}
tr.eval td{background:rgba(255,255,255,.02)}
.pill{display:inline-block;padding:0 6px;border-radius:999px;font-size:10px;border:1px solid var(--line);color:var(--mut);margin-left:6px}
.pop{cursor:help;border-bottom:1px dotted var(--mut)}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
button.f.tab{font-weight:600}
</style></head><body>
<div class="wrap">
  <h1>Mitochondrial SV caller comparison</h1>
  <p class="sub" id="subtitle"></p>
  <div class="cards" id="cards"></div>

  <h2>Runtime per caller</h2>
  <p class="legend">Wall-clock seconds per sample across the cohort (successful runs) — box = IQR,
    line = median, diamond = mean, whiskers = 1.5×IQR, dots = outliers. Hover the box for detail.</p>
  <div class="panel"><div id="chart-runtime"></div></div>
  <div class="panel" style="overflow:auto;margin-top:10px"><table id="runtime-table"></table></div>

  <h2>Accuracy per caller</h2>
  <p class="legend">Each caller as a binary deletion detector over the labelled cohort. Positives =
    detectable deletion events; negatives = wild-type / duplication / inversion samples (a deletion
    call there is a false positive). Forward-looking/ambiguous rows (wrap, dup-del, sub-size) are
    excluded from these numbers. Evaluation only — never gates the build.</p>
  <div class="panel" style="overflow:auto"><table id="accuracy"></table></div>

  <h2>Speed vs accuracy</h2>
  <p class="legend">Lower-left = fast; upper = more accurate (balanced accuracy = ½(sensitivity+specificity)).
    The reference caller is highlighted.</p>
  <div class="panel"><div id="chart-scatter"></div></div>

  <h2>Detection matrix — scenario &times; caller</h2>
  <p class="legend">Default = the <b>Deletions</b> suite (what we test today). Tabs add the
    forward-looking categories; <b>All</b> shows every scenario. Hover a sample name for what it tests.</p>
  <div class="controls" id="matrix-tabs"></div>
  <div id="eval-banner-slot"></div>
  <div class="panel scroll"><table class="matrix" id="matrix"></table></div>
  <div class="panel" style="overflow:auto;margin-top:10px"><table id="matrix-metrics"></table></div>
  <p class="legend"><span class="cell-yes" style="padding:1px 6px;border-radius:4px">detected (TP)</span>
    &nbsp; <span class="cell-fp" style="padding:1px 6px;border-radius:4px">false-positive deletion call</span>
    &nbsp; <span class="cell-amb">amb/gap/sub = evaluation-only row (not scored)</span></p>

  <h2>All calls</h2>
  <div class="controls">
    <input id="callfilter" type="text" placeholder="filter (sample / caller / type)…">
  </div>
  <div class="panel scroll"><table id="calls"></table></div>

  <div class="foot" id="foot"></div>
</div>
<div id="tip"></div>
<script>
const PAYLOAD = /*__PAYLOAD__*/;
const D = PAYLOAD.data, M = PAYLOAD.meta;
const CAL = D.callers;
const COLOR = {mitohpc:'#b692f6',eklipse:'#6ea8fe',mitosalt:'#2fbf71',
  splicebreak2:'#e3a008',mitomut:'#f178b6',mitoseek:'#4dd0e1'};
const col = c => COLOR[c] || '#9aa4b2';
const tip = document.getElementById('tip');
function showTip(e,html){tip.innerHTML=html;tip.style.opacity=1;
  tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY+12)+'px';}
function hideTip(){tip.style.opacity=0;}

// ---- subtitle + cards ----
document.getElementById('subtitle').textContent =
  `${M.n_samples} samples · ${CAL.length} callers · scope ${M.scope} · generated ${M.generated}`
  + (M.image ? ` · ${M.image}` : '');
const fastest = CAL.filter(c=>D.runtime[c].mean!=null)
  .sort((a,b)=>D.runtime[a].mean-D.runtime[b].mean)[0];
const pct1 = v => v==null ? '–' : Math.round(100*v)+'%';
const M_ALL = c => (D.metrics[c]||{}).all || {};
const M_DEL = c => (D.metrics[c]||{}).del || {};
const bySens = CAL.filter(c=>M_DEL(c).sensitivity!=null)
  .sort((a,b)=>M_DEL(b).sensitivity-M_DEL(a).sensitivity);
const bySpec = CAL.filter(c=>M_ALL(c).specificity!=null)
  .sort((a,b)=>M_ALL(b).specificity-M_ALL(a).specificity);
const nPos = (M_DEL(CAL[0]).n_pos)||0, nNeg = (M_ALL(CAL[0]).n_neg)||0;
const cards=[
  ['Callers', CAL.length],
  ['Scenarios', M.n_samples],
  ['Pos / neg', `${nPos} / ${nNeg}`],
  ['Fastest', fastest ? `${fastest} · ${D.runtime[fastest].mean}s` : '–'],
  ['Most sensitive', bySens[0] ? `${bySens[0]} · ${pct1(M_DEL(bySens[0]).sensitivity)}` : '–'],
  ['Best specificity', bySpec[0] ? `${bySpec[0]} · ${pct1(M_ALL(bySpec[0]).specificity)}` : '–'],
];
document.getElementById('cards').innerHTML = cards.map(([k,v])=>
  `<div class="card"><div class="k">${k}</div><div class="v">${v}</div></div>`).join('');

// ---- horizontal bar chart ----
function barChart(id, items, unit){
  const W=860, rowH=30, padL=120, padR=60, padT=10, padB=24;
  const H=padT+padB+items.length*rowH;
  const max=Math.max(1,...items.map(d=>d.v||0));
  const x=v=>padL+(W-padL-padR)*(v/max);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  // gridlines
  for(let i=0;i<=4;i++){const gx=padL+(W-padL-padR)*i/4, gv=(max*i/4);
    s+=`<line class="gl" x1="${gx}" y1="${padT}" x2="${gx}" y2="${H-padB}"/>`;
    s+=`<text x="${gx}" y="${H-8}" text-anchor="middle" fill="#9aa4b2">${gv.toFixed(gv<10?1:0)}</text>`;}
  items.forEach((d,i)=>{const y=padT+i*rowH+4, bw=x(d.v||0)-padL;
    s+=`<text x="${padL-8}" y="${y+13}" text-anchor="end">${d.label}</text>`;
    s+=`<rect class="bar" x="${padL}" y="${y}" width="${Math.max(0,bw)}" height="${rowH-12}" rx="4"
        fill="${d.color}" data-t="${encodeURIComponent(d.tip)}"/>`;
    s+=`<text x="${padL+Math.max(0,bw)+6}" y="${y+13}">${d.vlabel ?? d.v}${unit||''}</text>`;});
  s+=`</svg>`;
  const el=document.getElementById(id); el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{
    b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t)); b.onmouseleave=hideTip;});
}

// ---- horizontal box-and-whisker chart ----
function boxChart(id, items){
  items=items.filter(d=>d.b);
  const el=document.getElementById(id);
  if(!items.length){el.innerHTML='<span class="legend">no runtime data</span>';return;}
  const W=860,rowH=34,padL=120,padR=70,padT=10,padB=26,hb=9;
  const H=padT+padB+items.length*rowH;
  let mx=0; items.forEach(d=>{mx=Math.max(mx,d.b.whi,...(d.b.outliers||[]));}); mx=Math.max(1,mx);
  const X=v=>padL+(W-padL-padR)*(v/mx);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gx=X(mx*i/4);s+=`<line class="gl" x1="${gx}" y1="${padT}" x2="${gx}" y2="${H-padB}"/>`;
    s+=`<text x="${gx}" y="${H-9}" text-anchor="middle" fill="#9aa4b2">${(mx*i/4).toFixed(mx<10?1:0)}</text>`;}
  s+=`<text x="${(padL+W-padR)/2}" y="${H-0}" text-anchor="middle" fill="#9aa4b2">seconds / sample</text>`;
  items.forEach((d,i)=>{const b=d.b,y=padT+i*rowH+rowH/2-2,c=d.color;
    s+=`<text x="${padL-8}" y="${y+3}" text-anchor="end">${d.label}</text>`;
    s+=`<line x1="${X(b.wlo)}" y1="${y}" x2="${X(b.whi)}" y2="${y}" stroke="${c}" stroke-width="1" opacity=".6"/>`;
    s+=`<line x1="${X(b.wlo)}" y1="${y-5}" x2="${X(b.wlo)}" y2="${y+5}" stroke="${c}"/>`;
    s+=`<line x1="${X(b.whi)}" y1="${y-5}" x2="${X(b.whi)}" y2="${y+5}" stroke="${c}"/>`;
    s+=`<rect class="bar" x="${X(b.q1)}" y="${y-hb}" width="${Math.max(1,X(b.q3)-X(b.q1))}" height="${2*hb}" rx="2" fill="${c}" fill-opacity=".28" stroke="${c}" data-t="${encodeURIComponent(`<b>${d.label}</b> (n=${b.n})<br>median ${b.med}s · mean ${b.mean}s<br>IQR ${b.q1}–${b.q3}s · range ${b.min}–${b.max}s`)}"/>`;
    s+=`<line x1="${X(b.med)}" y1="${y-hb}" x2="${X(b.med)}" y2="${y+hb}" stroke="${c}" stroke-width="2"/>`;
    const mp=X(b.mean);s+=`<path d="M${mp},${y-5} L${mp+5},${y} L${mp},${y+5} L${mp-5},${y} Z" fill="var(--bg)" stroke="${c}" stroke-width="1.5"/>`;
    (b.outliers||[]).forEach(o=>{s+=`<circle cx="${X(o)}" cy="${y}" r="2.3" fill="${c}" fill-opacity=".7"/>`;});
  });
  s+=`</svg>`;el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
}
boxChart('chart-runtime', CAL.map(c=>({label:c, color:col(c), b:D.runtime[c].box})));
// runtime summary table (fastest median first)
(function(){
  const rows=CAL.map(c=>({c,b:D.runtime[c].box})).filter(x=>x.b).sort((a,b)=>a.b.med-b.b.med);
  let h=`<thead><tr><th>caller</th><th class="num">n</th><th class="num">median</th><th class="num">mean</th><th class="num">p25</th><th class="num">p75</th><th class="num">min</th><th class="num">max</th><th class="num">total</th></tr></thead><tbody>`;
  rows.forEach(d=>{const b=d.b;h+=`<tr><td>${d.c}</td><td class="num">${b.n}</td><td class="num">${b.med}s</td><td class="num">${b.mean}s</td><td class="num">${b.q1}s</td><td class="num">${b.q3}s</td><td class="num">${b.min}s</td><td class="num">${b.max}s</td><td class="num">${D.runtime[d.c].total}s</td></tr>`;});
  h+=`</tbody>`;document.getElementById('runtime-table').innerHTML=h;
})();

// ---- accuracy table (overall, all scored scenarios) ----
const fmtPct = (v,n,d) => v==null ? '<span class="muted">–</span>'
  : `${Math.round(100*v)}%${(n!=null&&d!=null)?` <span class="muted">(${n}/${d})</span>`:''}`;
const fmtNum = v => v==null ? '<span class="muted">–</span>' : (Math.round(v*100)/100).toFixed(2);
(function(){
  const rows=CAL.slice().sort((a,b)=>((M_ALL(b).balanced_acc??-1)-(M_ALL(a).balanced_acc??-1)));
  let h=`<thead><tr><th>caller</th><th class="num">sensitivity</th><th class="num">specificity</th>`
    +`<th class="num">precision</th><th class="num">F1</th><th class="num">bal.acc</th><th class="num">MCC</th>`
    +`<th class="num">FP</th></tr></thead><tbody>`;
  rows.forEach(c=>{const a=M_ALL(c),dl=M_DEL(c);
    const lbl=c==='mitohpc'?`<span class="tag ref">${c}</span>`:c;
    h+=`<tr><td>${lbl}</td>`
      +`<td class="num">${fmtPct(dl.sensitivity,dl.tp,dl.n_pos)}</td>`
      +`<td class="num">${fmtPct(a.specificity,a.tn,a.n_neg)}</td>`
      +`<td class="num">${fmtPct(a.precision)}</td><td class="num">${fmtNum(a.f1)}</td>`
      +`<td class="num">${fmtPct(a.balanced_acc)}</td><td class="num">${fmtNum(a.mcc)}</td>`
      +`<td class="num">${a.fp||0}${(a.common_fp_samples&&a.common_fp_samples.length)?' <span class="muted">('+a.common_fp_samples.length+' common)</span>':''}</td></tr>`;});
  h+=`</tbody>`;document.getElementById('accuracy').innerHTML=h;
})();

// ---- scatter: mean runtime (x) vs balanced accuracy (y) ----
(function(){
  const W=860,H=320,padL=50,padR=20,padT=14,padB=40;
  const pts=CAL.map(c=>({c, x:D.runtime[c].mean, y:(M_ALL(c).balanced_acc!=null?100*M_ALL(c).balanced_acc:null)}))
    .filter(p=>p.x!=null && p.y!=null);
  const maxX=Math.max(1,...pts.map(p=>p.x)), maxY=100;
  const X=v=>padL+(W-padL-padR)*(v/maxX), Y=v=>H-padB-(H-padT-padB)*(v/maxY);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gy=Y(maxY*i/4);
    s+=`<line class="gl" x1="${padL}" y1="${gy}" x2="${W-padR}" y2="${gy}"/>`;
    s+=`<text x="${padL-6}" y="${gy+3}" text-anchor="end" fill="#9aa4b2">${(maxY*i/4)}%</text>`;}
  s+=`<text x="${(W)/2}" y="${H-6}" text-anchor="middle" fill="#9aa4b2">mean runtime (s)</text>`;
  s+=`<text transform="translate(14 ${H/2}) rotate(-90)" text-anchor="middle" fill="#9aa4b2">balanced accuracy</text>`;
  pts.forEach(p=>{const r=p.c==='mitohpc'?7:6;
    s+=`<circle class="bar" cx="${X(p.x)}" cy="${Y(p.y)}" r="${r}" fill="${col(p.c)}"
        stroke="${p.c==='mitohpc'?'#fff':'none'}" stroke-width="1.5"
        data-t="${encodeURIComponent(`<b>${p.c}</b><br>${p.x}s · ${Math.round(p.y)}% balanced acc`)}"/>`;
    s+=`<text x="${X(p.x)+9}" y="${Y(p.y)+3}">${p.c}</text>`;});
  s+=`</svg>`;
  const el=document.getElementById('chart-scatter'); el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
})();

// ---- detection matrix (category tabs) ----
const CATS = (D.categories||[]).concat([{key:'all',label:'All',blurb:'every scenario',eval_only:false}]);
let mtab = (CATS.find(c=>c.key==='del')||CATS[0]||{key:'all'}).key;
function tabCat(){ return CATS.find(c=>c.key===mtab) || {}; }
function renderMatrix(){
  const cat=tabCat();
  const rows=D.matrix.filter(m=> mtab==='all' ? true : m.category===mtab);
  // eval-only banner
  document.getElementById('eval-banner-slot').innerHTML = cat.eval_only
    ? `<div class="eval-banner">Evaluation-only — forward-looking / ambiguous scenarios; never gates the build. ${cat.blurb||''}</div>` : '';
  let h=`<thead><tr><th>sample</th><th>truth event</th>`+
    CAL.map(c=>`<th>${c==='mitohpc'?'<span class="tag ref">'+c+'</span>':c}</th>`).join('')+`</tr></thead><tbody>`;
  rows.forEach(m=>{
    const pop = m.popup ? ` data-pop="${encodeURIComponent(m.popup)}"` : '';
    const nm = m.popup ? `<span class="pop"${pop}>${m.sample}</span>` : m.sample;
    h+=`<tr class="${m.eval_only?'eval':''}"><td>${nm}</td><td>${m.label}${m.eval_only?'<span class="pill">eval</span>':''}</td>`;
    CAL.forEach(c=>{
      let cls='cell-no',txt='·';
      if(m.klass==='positive'){ if(m.detected&&m.detected[c]){cls='cell-yes';txt='✓';} }
      else if(m.klass==='negative'){ const f=m.fp&&m.fp[c]; if(f){cls='cell-fp';txt=f.common?'FP✱':'FP';} }
      else { const called=(m.detected&&m.detected[c])||(m.fp&&m.fp[c]);
        if(called){cls='cell-amb';txt=m.reason==='wrap'?'wrap':m.reason==='knownfp'?'gap':m.reason==='sub'?'sub':'amb';} }
      h+=`<td class="c ${cls}">${txt}</td>`;});
    h+=`</tr>`;});
  h+=`</tbody>`;
  document.getElementById('matrix').innerHTML=h;
  document.querySelectorAll('#matrix .pop').forEach(el=>{
    el.onmousemove=e=>showTip(e,decodeURIComponent(el.dataset.pop)); el.onmouseleave=hideTip;});
  renderMatrixMetrics();
}
function renderMatrixMetrics(){
  const key=mtab, get=c=>(D.metrics[c]||{})[key]||{};
  const any=get(CAL[0]), hasPos=(any.n_pos||0)>0, hasNeg=(any.n_neg||0)>0, full=key==='all';
  const slot=document.getElementById('matrix-metrics');
  if(!hasPos && !hasNeg){ slot.innerHTML='<p class="legend" style="padding:8px 4px">Evaluation-only category — forward-looking / ambiguous scenarios, not scored. The matrix above shows what each caller did.</p>'; return; }
  let cols=['caller']; if(hasPos)cols.push('sensitivity'); if(hasNeg)cols.push('specificity','FP'); if(full)cols.push('precision','F1','bal.acc','MCC');
  let h='<thead><tr>'+cols.map((c,i)=>`<th class="${i?'num':''}">${c}</th>`).join('')+'</tr></thead><tbody>';
  const order=CAL.slice().sort((a,b)=>{const x=hasPos?'sensitivity':'specificity';return (get(b)[x]??-1)-(get(a)[x]??-1);});
  order.forEach(c=>{const m=get(c),lbl=c==='mitohpc'?`<span class="tag ref">${c}</span>`:c;
    h+=`<tr><td>${lbl}</td>`;
    if(hasPos)h+=`<td class="num">${fmtPct(m.sensitivity,m.tp,m.n_pos)}</td>`;
    if(hasNeg){h+=`<td class="num">${fmtPct(m.specificity,m.tn,m.n_neg)}</td>`;
      h+=`<td class="num">${m.fp||0}${(m.fp_samples&&m.fp_samples.length)?' <span class="muted">('+m.fp_samples.join(', ')+')</span>':''}</td>`;}
    if(full)h+=`<td class="num">${fmtPct(m.precision)}</td><td class="num">${fmtNum(m.f1)}</td><td class="num">${fmtPct(m.balanced_acc)}</td><td class="num">${fmtNum(m.mcc)}</td>`;
    h+=`</tr>`;});
  slot.innerHTML=h+'</tbody>';
}
const mc=document.getElementById('matrix-tabs');
CATS.forEach(cat=>{
  const b=document.createElement('button');b.className='f tab'+(cat.key===mtab?' on':'');
  b.textContent=cat.label+(cat.eval_only?' ⚑':'');b.title=cat.blurb||'';
  b.onclick=()=>{mtab=cat.key;mc.querySelectorAll('button').forEach(x=>x.classList.remove('on'));
    b.classList.add('on');renderMatrix();};mc.appendChild(b);});
renderMatrix();

// ---- calls table ----
function renderCalls(q){
  q=(q||'').toLowerCase();
  const rows=D.calls.filter(c=>!q ||
    `${c.sample} ${c.caller} ${c.sv_type}`.toLowerCase().includes(q));
  let h=`<thead><tr><th>sample</th><th>caller</th><th>type</th><th>bp5</th><th>bp3</th>
    <th>svlen</th><th>support</th><th>het</th><th>common</th><th>extra</th></tr></thead><tbody>`;
  rows.forEach(c=>{h+=`<tr><td>${c.sample}</td><td style="color:${col(c.caller)}">${c.caller}</td>
    <td>${c.sv_type||''}</td><td>${c.bp5??''}</td><td>${c.bp3??''}</td><td>${c.svlen??''}</td>
    <td>${c.support??''}</td><td>${c.het==null?'':(+c.het).toFixed(3)}</td>
    <td>${c.common_deletion==1?'✓':''}</td><td class="muted">${c.extra||''}</td></tr>`;});
  h+=`</tbody>`;document.getElementById('calls').innerHTML=h;
}
renderCalls('');
document.getElementById('callfilter').oninput=e=>renderCalls(e.target.value);

document.getElementById('foot').innerHTML =
  `Generated by <code>pipeline/make_report.py</code>. Detection is an evaluation/comparison of `
  + `third-party callers (not a pass/fail gate). Runtimes are on small test inputs — relative, `
  + `not absolute. Reference caller: <span class="tag ref">mitohpc</span> (MitoHPC).`;
</script></body></html>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--runtime", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--scenarios", default="", help="scenarios.json (categories + BAM popups)")
    ap.add_argument("--samples", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--scope", default="full")
    ap.add_argument("--image", default="")
    ap.add_argument("--generated", default="")
    args = ap.parse_args(argv)

    truth = load_truth(args.truth)
    scen = args.scenarios or os.path.join(os.path.dirname(args.truth), "scenarios.json")
    categories, descriptions = load_scenarios(scen)
    calls = load_calls(args.calls)
    runtime = load_runtime(args.runtime)
    samples = [s for s in args.samples.replace(",", " ").split() if s]
    if not samples:
        samples = sorted({c["sample"] for c in calls} | {r["sample"] for r in runtime})

    data = build(calls, runtime, truth, samples, categories, descriptions)
    meta = {"scope": args.scope, "image": args.image,
            "generated": args.generated or "(unstamped)", "n_samples": len(samples)}
    html = render_html(data, meta)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(html)
    sys.stderr.write("make_report: wrote %s (%d calls, %d samples, %d callers)\n"
                     % (args.out, len(calls), len(samples), len(data["callers"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
