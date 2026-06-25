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

GEN_TOL = 250
REAL_TRUTH = {
    "spike_del4977_h20": [("del", 8469, 13447)],
    "NA12718": [("none", None, None)],
    "NA12748": [("none", None, None)],
    "NA12775": [("none", None, None)],
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
    truth = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                truth.setdefault(f[0], []).append((f[1], _num(f[2]), _num(f[3])))
    truth.update(REAL_TRUTH)
    return truth


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


def build(calls, runtime, truth, samples):
    callers = list(parsers.CALLERS)
    calls_by_sample = {}
    for c in calls:
        calls_by_sample.setdefault(c["sample"], []).append(c)

    # ---- scenario x caller detection matrix ----
    matrix = []
    for sample in samples:
        events = truth.get(base_sample(sample), truth.get(sample, []))
        carries_common = any(k in ("del", "delwrap")
                             and parsers.is_common_deletion(b5, b3)
                             for (k, b5, b3) in events)
        scalls = calls_by_sample.get(sample, [])
        common_fp_here = sorted({c["caller"] for c in scalls
                                 if parsers.is_common_deletion(c["bp5"], c["bp3"])}) \
            if not carries_common else []
        if not events or all(e[0] == "none" for e in events):
            matrix.append({"sample": sample, "label": "wild-type (no SV)",
                           "kind": "none", "common": False,
                           "detected": {}, "fp": common_fp_here})
        for (kind, e5, e3) in events:
            if kind not in ("del", "delwrap", "dup"):
                continue
            is_common = kind != "dup" and parsers.is_common_deletion(e5, e3)
            det = {}
            for c in scalls:
                if c["sv_type"] not in ("deletion", "duplication"):
                    continue
                hit = (parsers.is_common_deletion(c["bp5"], c["bp3"]) if is_common
                       else _match(c["bp5"], c["bp3"], e5, e3))
                if hit:
                    det[c["caller"]] = True
            label = "%s %s-%s%s" % ("dup" if kind == "dup" else "del",
                                    e5, e3, " · COMMON" if is_common else "")
            matrix.append({"sample": sample, "label": label, "kind": kind,
                           "common": is_common, "detected": det,
                           "fp": common_fp_here if kind != "dup" else []})

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

    # ---- detection summary per caller (from the matrix, deletion events) ----
    det_counts = {c: 0 for c in callers}
    n_truth_dels = 0
    for m in matrix:
        if m["kind"] in ("del", "delwrap"):
            n_truth_dels += 1
            for c in m["detected"]:
                det_counts[c] = det_counts.get(c, 0) + 1
    det_summary = {}
    for c in callers:
        n_calls = sum(1 for x in calls if x["caller"] == c)
        fp = sorted({m["sample"] for m in matrix if c in (m.get("fp") or [])})
        det_summary[c] = {
            "n_truth_dels": n_truth_dels,
            "n_detected": det_counts.get(c, 0),
            "sensitivity": round(100.0 * det_counts.get(c, 0) / n_truth_dels) if n_truth_dels else 0,
            "n_calls": n_calls,
            "common_fp_samples": fp,
            "ran": runtime_summary[c]["n"] > 0 or n_calls > 0,
        }

    return {
        "callers": callers,
        "samples": samples,
        "matrix": matrix,
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

  <h2>Sensitivity per caller</h2>
  <p class="legend">Share of truth deletions detected (any matching call; common deletion within
    &plusmn;80&nbsp;bp, others &plusmn;250&nbsp;bp). Evaluation only — callers are third-party.</p>
  <div class="panel"><div id="chart-sens"></div></div>

  <h2>Speed vs sensitivity</h2>
  <p class="legend">Lower-left = fast; upper = more sensitive. The reference caller is highlighted.</p>
  <div class="panel"><div id="chart-scatter"></div></div>

  <h2>Detection matrix — scenario &times; caller</h2>
  <div class="controls" id="matrix-controls"></div>
  <div class="panel scroll"><table class="matrix" id="matrix"></table></div>
  <p class="legend"><span class="cell-yes" style="padding:1px 6px;border-radius:4px">detected</span>
    &nbsp; <span class="cell-fp" style="padding:1px 6px;border-radius:4px">common-deletion call on a non-carrier (false positive)</span></p>

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
const sens = CAL.slice().sort((a,b)=>D.detection[b].n_detected-D.detection[a].n_detected);
const cards=[
  ['Callers', CAL.length],
  ['Samples', M.n_samples],
  ['Truth deletions', (D.detection[CAL[0]]||{}).n_truth_dels ?? 0],
  ['Fastest', fastest ? `${fastest} · ${D.runtime[fastest].mean}s` : '–'],
  ['Most sensitive', sens[0] ? `${sens[0]} · ${D.detection[sens[0]].n_detected}/${D.detection[sens[0]].n_truth_dels}` : '–'],
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

barChart('chart-sens', CAL.map(c=>{const d=D.detection[c];
  return {label:c, v:d.n_detected, color:col(c),
    vlabel:`${d.n_detected}/${d.n_truth_dels}`,
    tip:`<b>${c}</b><br>${d.n_detected}/${d.n_truth_dels} truth deletions (${d.sensitivity}%)<br>${d.n_calls} total calls`
      + (d.common_fp_samples.length?`<br><span style="color:#e5534b">common-del FP on: ${d.common_fp_samples.join(', ')}</span>`:'')};}));

// ---- scatter: mean runtime (x) vs sensitivity (y) ----
(function(){
  const W=860,H=320,padL=50,padR=20,padT=14,padB=40;
  const pts=CAL.map(c=>({c, x:D.runtime[c].mean, y:D.detection[c].sensitivity}))
    .filter(p=>p.x!=null);
  const maxX=Math.max(1,...pts.map(p=>p.x)), maxY=100;
  const X=v=>padL+(W-padL-padR)*(v/maxX), Y=v=>H-padB-(H-padT-padB)*(v/maxY);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gy=Y(maxY*i/4);
    s+=`<line class="gl" x1="${padL}" y1="${gy}" x2="${W-padR}" y2="${gy}"/>`;
    s+=`<text x="${padL-6}" y="${gy+3}" text-anchor="end" fill="#9aa4b2">${(maxY*i/4)}%</text>`;}
  s+=`<text x="${(W)/2}" y="${H-6}" text-anchor="middle" fill="#9aa4b2">mean runtime (s)</text>`;
  pts.forEach(p=>{const r=p.c==='mitohpc'?7:6;
    s+=`<circle class="bar" cx="${X(p.x)}" cy="${Y(p.y)}" r="${r}" fill="${col(p.c)}"
        stroke="${p.c==='mitohpc'?'#fff':'none'}" stroke-width="1.5"
        data-t="${encodeURIComponent(`<b>${p.c}</b><br>${p.x}s · ${p.y}% sensitivity`)}"/>`;
    s+=`<text x="${X(p.x)+9}" y="${Y(p.y)+3}">${p.c}</text>`;});
  s+=`</svg>`;
  const el=document.getElementById('chart-scatter'); el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
})();

// ---- detection matrix ----
let mfilter='all';
function renderMatrix(){
  const rows=D.matrix.filter(m=> mfilter==='all' ? true :
    mfilter==='common' ? m.common :
    mfilter==='del' ? (m.kind==='del'||m.kind==='delwrap') : true);
  let h=`<thead><tr><th>sample</th><th>truth event</th>`+
    CAL.map(c=>`<th>${c==='mitohpc'?'<span class="tag ref">'+c+'</span>':c}</th>`).join('')+`</tr></thead><tbody>`;
  rows.forEach(m=>{
    h+=`<tr><td>${m.sample}</td><td>${m.label}</td>`;
    CAL.forEach(c=>{
      let cls='cell-no',txt='·';
      if(m.detected && m.detected[c]){cls='cell-yes';txt='✓';}
      if((m.fp||[]).includes(c)){cls='cell-fp';txt='FP';}
      h+=`<td class="c ${cls}">${txt}</td>`;});
    h+=`</tr>`;});
  h+=`</tbody>`;
  document.getElementById('matrix').innerHTML=h;
}
const mc=document.getElementById('matrix-controls');
[['all','all events'],['del','deletions'],['common','common deletion only']].forEach(([k,lab])=>{
  const b=document.createElement('button');b.className='f'+(k==='all'?' on':'');b.textContent=lab;
  b.onclick=()=>{mfilter=k;mc.querySelectorAll('button').forEach(x=>x.classList.remove('on'));
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
    ap.add_argument("--samples", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--scope", default="full")
    ap.add_argument("--image", default="")
    ap.add_argument("--generated", default="")
    args = ap.parse_args(argv)

    truth = load_truth(args.truth)
    calls = load_calls(args.calls)
    runtime = load_runtime(args.runtime)
    samples = [s for s in args.samples.replace(",", " ").split() if s]
    if not samples:
        samples = sorted({c["sample"] for c in calls} | {r["sample"] for r in runtime})

    data = build(calls, runtime, truth, samples)
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
