#!/usr/bin/env python3
"""Build the interactive LOD report (self-contained HTML, no external deps).

Reads the LOD sweep + the aggregated lod_cells.tsv / lod_fits.tsv and renders
methods, per-caller detection-probability curves, a detection heatmap, a LOD
summary table, runtime, pipeline-vs-circular arm comparison, and an
interpretation guide. Interactive selectors (arm / deletion / depth) drive the
figures; pure inline SVG + vanilla JS.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lod_stats as S  # noqa: E402

CALLER_ORDER = ["mitohpc", "eklipse", "mitosalt", "splicebreak2", "mitomut", "mitoseek"]


def _f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def _round_box(b, nd=2):
    """Round a boxstats dict for compact JSON (cap outliers to keep size sane)."""
    r = {k: (round(v, nd) if isinstance(v, float) else v)
         for k, v in b.items() if k != "outliers"}
    outs = sorted(round(x, nd) for x in b["outliers"])
    r["outliers"] = outs[:40]
    return r


def load_tsv(path):
    if not path or not os.path.isfile(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(cells, fits, sweep, arms):
    callers = [c for c in CALLER_ORDER if any(r["caller"] == c for r in cells)]
    arms = [a for a in arms if any(r["arm"] == a for r in cells)] or \
        sorted({r["arm"] for r in cells})
    variants = sorted({r["variant"] for r in cells})
    depths = sorted({int(r["depth"]) for r in cells})
    vafs = sorted({_f(r["vaf"]) for r in cells if _f(r["vaf"]) is not None})

    cell_map = {}   # (arm,variant,depth,caller) -> [{vaf,rate,lo,hi,k,n,pass}]
    for r in cells:
        key = (r["arm"], r["variant"], int(r["depth"]), r["caller"])
        cell_map.setdefault(key, []).append({
            "vaf": _f(r["vaf"]), "rate": _f(r["det_rate"]),
            "lo": _f(r["det_lo"]), "hi": _f(r["det_hi"]),
            "k": int(r["det_k"]), "n": int(r["det_n"]),
            "passrate": _f(r.get("pass_rate")),
        })
    for v in cell_map.values():
        v.sort(key=lambda d: (d["vaf"] if d["vaf"] is not None else 0))

    fit_map = {}
    for r in fits:
        fit_map["%s|%s|%s|%s" % (r["arm"], r["variant"], r["depth"], r["caller"])] = {
            "emp_transition": _f(r["emp_transition"]), "emp_reliable": _f(r["emp_reliable"]),
            "near_separable": r["near_separable"] == "1",
            "lod50": _f(r["lod50"]), "lod95": _f(r["lod95"]),
            "lod95_lo": _f(r["lod95_lo"]), "lod95_hi": _f(r["lod95_hi"]),
        }

    # runtime DISTRIBUTIONS keyed "arm|depth|caller" (arm in all/pipeline/circular,
    # depth in all/<value>) from the raw per-cell runtime_s — boxplots, the summary
    # table, AND the runtime-vs-depth view (how depth warps caller runtime).
    rt = defaultdict(list)
    for r in sweep:
        s = _f(r.get("runtime_s"))
        if s is None or r.get("status") != "ok":
            continue
        c, arm = r["caller"], r["arm"]
        dp = _f(r.get("depth"))
        dp = int(dp) if dp is not None else None
        for a in ("all", arm):
            for d in ("all", dp):
                if d is not None:
                    rt[(a, d, c)].append(s)
    runtime_box = {}
    for (arm, depth, caller), vals in rt.items():
        b = S.boxstats(vals)
        if b:
            runtime_box["%s|%s|%s" % (arm, depth, caller)] = _round_box(b)

    return {
        "callers": callers, "arms": arms, "variants": variants,
        "depths": depths, "vafs": vafs,
        "cells": {"|".join(map(str, k)): v for k, v in cell_map.items()},
        "fits": fit_map, "runtime": runtime_box,
        "n_rows": len(sweep),
    }


def render(data, meta):
    payload = json.dumps({"data": data, "meta": meta}).replace("</", "<\\/")
    return TEMPLATE.replace("/*__PAYLOAD__*/", payload)


TEMPLATE = r"""<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>mtDNA SV caller — limit of detection</title>
<style>
:root{--bg:#0f1117;--panel:#171a23;--panel2:#1f2430;--ink:#e6e9ef;--mut:#9aa4b2;
--line:#2a2f3a;--accent:#6ea8fe;--good:#2fbf71;--bad:#e5534b;--ref:#b692f6;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1080px;margin:0 auto;padding:26px 20px 90px}
h1{font-size:24px;margin:0 0 4px}h2{font-size:18px;margin:32px 0 10px;border-bottom:1px solid var(--line);padding-bottom:6px}
h3{font-size:15px;margin:18px 0 6px;color:var(--mut)}
.sub{color:var(--mut);margin:0 0 16px}
p{color:#d2d7e0}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:12px 0}
.legend{color:var(--mut);font-size:12px;margin:6px 0}
.controls{display:flex;gap:14px;flex-wrap:wrap;align-items:center;margin:10px 0}
label{font-size:12px;color:var(--mut);margin-right:6px}
select{background:var(--panel2);color:var(--ink);border:1px solid var(--line);border-radius:8px;padding:5px 9px;font-size:13px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{color:var(--mut);font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.tag{display:inline-block;padding:1px 7px;border-radius:999px;font-size:11px;border:1px solid var(--line)}
.tag.ref{color:var(--ref);border-color:var(--ref)}
svg text{fill:var(--ink);font-size:11px}svg .gl{stroke:var(--line);stroke-dasharray:2 3;opacity:.5}
.dot{cursor:pointer}.swatch{width:11px;height:11px;border-radius:3px;display:inline-block;vertical-align:-1px;margin-right:4px}
#tip{position:fixed;pointer-events:none;background:#000d;border:1px solid var(--line);border-radius:8px;padding:6px 9px;font-size:12px;opacity:0;transition:opacity .1s;z-index:9}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:14px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:13px 15px}
.card .k{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em}.card .v{font-size:21px;font-weight:600;margin-top:3px}
.foot{color:var(--mut);font-size:12px;margin-top:30px;border-top:1px solid var(--line);padding-top:14px}
code{background:var(--panel2);padding:1px 5px;border-radius:5px;font-size:12px}
</style></head><body><div class="wrap">
<h1>Mitochondrial SV calling — limit of detection</h1>
<p class="sub" id="subtitle"></p>

<h2>What this measures</h2>
<p>For each deletion, we simulate reads from a mixture of wild-type and
deletion-bearing mtDNA molecules at a known <b>heteroplasmy</b> (the fraction of
molecules carrying the deletion) and <b>depth</b>, then ask each caller whether it
recovers the deletion. Repeating across replicates gives a <b>detection
probability</b> at every heteroplasmy, and the <b>limit of detection (LOD)</b> is
the heteroplasmy at which that probability crosses 50% / 95%. Lower LOD = more
sensitive. This is an evaluation/comparison of third-party callers — not a
pass/fail gate.</p>
<div class="cards" id="cards"></div>

<h2>Detection-probability curves</h2>
<div class="controls" id="controls"></div>
<div class="panel"><div id="curve"></div><div id="legend" style="margin-top:6px"></div></div>
<p class="legend">Detection rate (with 95% Wilson interval ticks) vs heteroplasmy, one line per caller,
for the selected input arm / deletion / depth. The 95% LOD is where a curve reaches 0.95.</p>

<h2>Detection heatmap</h2>
<div class="panel"><div id="heat"></div></div>
<p class="legend">Detection rate per caller × heteroplasmy (selected arm/deletion/depth). Greener = detected more often.</p>

<h2>LOD summary</h2>
<div class="panel" style="overflow:auto"><table id="summary"></table></div>
<p class="legend"><b>empirical transition</b> = highest heteroplasmy still detected &lt;50% of the time;
<b>reliable</b> = lowest at ≥90%. <b>LOD50/LOD95</b> = Firth-logistic fit (★ = near-separable, treat as approximate).
All values are heteroplasmy %.</p>

<h2>Runtime</h2>
<div class="panel"><div id="runtime"></div></div>
<p class="legend">Per-cell wall-clock seconds per caller for the selected arm (all depths) — box = IQR,
line = median, diamond = mean, whiskers = 1.5×IQR, dots = outliers. Relative on these inputs, not absolute.</p>

<h3 style="margin:18px 0 6px">Runtime vs depth — how sequencing depth warps cost</h3>
<div class="panel"><div id="runtime-depth"></div></div>
<p class="legend">Median seconds/cell per caller across simulated depth (selected arm). Real-world mtDNA depth
runs to thousands× — the alignment-heavy callers steepen sharply toward 2000×, while the lightweight ones
stay flat. This is the single biggest driver of LOD-sweep cost.</p>
<div class="panel" style="overflow:auto;margin-top:10px"><table id="runtime-table"></table></div>
<p class="legend">Median seconds/cell per caller × depth (selected arm).</p>

<h2>Pipeline vs circular-aware input</h2>
<p>Two input preparations are compared: <b>pipeline</b> re-normalises the reads with
<code>bwa mem</code> exactly as the production pipeline does; <b>circular</b> feeds
MitoHPC's circular-aware <code>minimap2 + circSam.pl</code> BAM directly. The table
shows each caller's 95% LOD under both (at the selected deletion / depth).</p>
<div class="panel" style="overflow:auto"><table id="armcmp"></table></div>

<h2>How to read this report</h2>
<ul>
<li><b>LOD curve</b> — a caller whose curve climbs to 1.0 at a lower heteroplasmy is more sensitive. A flat/low curve means it misses the deletion even at high heteroplasmy on this input.</li>
<li><b>LOD95</b> — the headline sensitivity number: the heteroplasmy you need before that caller catches the deletion ≥95% of the time. The empirical transition is the robust read-out when the fit is flagged near-separable (★).</li>
<li><b>Depth</b> — sensitivity usually improves with depth; compare the same caller across the depth selector.</li>
<li><b>Deletion</b> — <code>del4977</code> sits in a 13&nbsp;bp direct repeat (the common deletion); <code>del6000</code> is a non-repeat deletion. Repeat-mediated junctions shift the called breakpoint (absorbed by the 30&nbsp;bp detection tolerance).</li>
<li><b>Arm</b> — if a caller does much better under one input arm, that input prep matters for it; the pipeline arm reflects production behaviour.</li>
<li><b>Specificity</b> — the <code>vaf=0</code> column is the blank: detections there are false positives.</li>
</ul>

<div class="foot" id="foot"></div>
</div><div id="tip"></div>
<script>
const P=/*__PAYLOAD__*/;const D=P.data,M=P.meta;
const COLOR={mitohpc:'#b692f6',eklipse:'#6ea8fe',mitosalt:'#2fbf71',splicebreak2:'#e3a008',mitomut:'#f178b6',mitoseek:'#4dd0e1'};
const col=c=>COLOR[c]||'#9aa4b2';const pct=v=>v==null?'—':(100*v).toFixed(v<0.1?1:0);
const tip=document.getElementById('tip');
const showTip=(e,h)=>{tip.innerHTML=h;tip.style.opacity=1;tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY+12)+'px';};
const hideTip=()=>tip.style.opacity=0;
let sel={arm:D.arms[0],variant:D.variants[0],depth:D.depths[D.depths.length-1]};

document.getElementById('subtitle').textContent=
 `${D.callers.length} callers · arms: ${D.arms.join(', ')} · deletions: ${D.variants.join(', ')} · `+
 `depths: ${D.depths.join('/')}× · ${D.vafs.length} heteroplasmy levels · generated ${M.generated}`;

// cards: most sensitive (lowest LOD95) at production depth, fastest
function fitOf(arm,variant,depth,caller){return D.fits[[arm,variant,depth,caller].join('|')]||{};}
function cellsOf(arm,variant,depth,caller){return D.cells[[arm,variant,depth,caller].join('|')]||[];}
(function(){
  const pd=sel.depth, arm=D.arms[0], v=D.variants[0];
  const ranked=D.callers.map(c=>({c,lod:fitOf(arm,v,pd,c).lod95})).filter(x=>x.lod!=null).sort((a,b)=>a.lod-b.lod);
  const rt=D.callers.map(c=>({c,b:D.runtime['all|all|'+c]})).filter(x=>x.b!=null).sort((a,b)=>a.b.med-b.b.med);
  const cards=[['Callers',D.callers.length],['Deletions',D.variants.length],
    ['Most sensitive',ranked[0]?`${ranked[0].c} · ${pct(ranked[0].lod)}%`:'—'],
    ['Fastest (median)',rt[0]?`${rt[0].c} · ${rt[0].b.med}s`:'—']];
  document.getElementById('cards').innerHTML=cards.map(([k,v])=>`<div class="card"><div class="k">${k}</div><div class="v">${v}</div></div>`).join('');
})();

// controls
function mkSel(id,label,opts,cur,fmt){
  return `<span><label>${label}</label><select id="${id}">`+
    opts.map(o=>`<option value="${o}" ${o==cur?'selected':''}>${fmt?fmt(o):o}</option>`).join('')+`</select></span>`;
}
document.getElementById('controls').innerHTML=
  mkSel('selArm','input arm',D.arms,sel.arm)+
  mkSel('selVar','deletion',D.variants,sel.variant)+
  mkSel('selDep','depth',D.depths,sel.depth,d=>d+'×');
['selArm','selVar','selDep'].forEach(id=>document.getElementById(id).onchange=e=>{
  if(id=='selArm')sel.arm=e.target.value; if(id=='selVar')sel.variant=e.target.value; if(id=='selDep')sel.depth=+e.target.value; redraw();});

function redraw(){drawCurve();drawHeat();drawSummary();drawRuntime();drawRuntimeDepth();drawArmCmp();}
function rtbox(arm,depth,c){return D.runtime[arm+'|'+depth+'|'+c];}

function drawCurve(){
  const W=860,H=340,L=52,R=18,T=14,B=42;const vafs=D.vafs;
  const maxV=Math.max(...vafs);const X=v=>L+(W-L-R)*(maxV?v/maxV:0);const Y=p=>H-B-(H-T-B)*p;
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=5;i++){const gy=Y(i/5);s+=`<line class="gl" x1="${L}" y1="${gy}" x2="${W-R}" y2="${gy}"/>`;
    s+=`<text x="${L-6}" y="${gy+3}" text-anchor="end" fill="#9aa4b2">${(i*20)}%</text>`;}
  vafs.forEach(v=>{s+=`<text x="${X(v)}" y="${H-8}" text-anchor="middle" fill="#9aa4b2">${pct(v)}</text>`;});
  s+=`<text x="${(L+W-R)/2}" y="${H-1}" text-anchor="middle" fill="#9aa4b2">heteroplasmy</text>`;
  let leg='';
  D.callers.forEach(c=>{
    const cl=cellsOf(sel.arm,sel.variant,sel.depth,c);if(!cl.length)return;
    const pts=cl.filter(d=>d.rate!=null);
    let path=pts.map((d,i)=>`${i?'L':'M'}${X(d.vaf)},${Y(d.rate)}`).join(' ');
    s+=`<path d="${path}" fill="none" stroke="${col(c)}" stroke-width="2"/>`;
    pts.forEach(d=>{
      s+=`<line x1="${X(d.vaf)}" y1="${Y(d.lo)}" x2="${X(d.vaf)}" y2="${Y(d.hi)}" stroke="${col(c)}" stroke-width="1" opacity=".4"/>`;
      s+=`<circle class="dot" cx="${X(d.vaf)}" cy="${Y(d.rate)}" r="3.2" fill="${col(c)}" data-t="${encodeURIComponent(`<b>${c}</b> @ ${pct(d.vaf)}% het<br>detected ${d.k}/${d.n} = ${(100*d.rate).toFixed(0)}%<br>95% CI ${(100*d.lo).toFixed(0)}–${(100*d.hi).toFixed(0)}%`)}"/>`;
    });
    const f=fitOf(sel.arm,sel.variant,sel.depth,c);
    leg+=`<span style="margin-right:14px"><span class="swatch" style="background:${col(c)}"></span>${c} · LOD95 ${f.lod95!=null?pct(f.lod95)+'%':'—'}${f.near_separable?'★':''}</span>`;
  });
  s+=`</svg>`;const el=document.getElementById('curve');el.innerHTML=s;
  el.querySelectorAll('.dot').forEach(d=>{d.onmousemove=e=>showTip(e,decodeURIComponent(d.dataset.t));d.onmouseleave=hideTip;});
  document.getElementById('legend').innerHTML=leg;
}

function drawHeat(){
  const vafs=D.vafs;const cw=Math.max(34,Math.min(70,(820-150)/vafs.length));
  const rh=26,L=120,T=22;const W=L+vafs.length*cw+10,H=T+D.callers.length*rh+6;
  const grn=r=>r==null?'#222733':`rgba(47,191,113,${0.12+0.78*r})`;
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  vafs.forEach((v,i)=>s+=`<text x="${L+i*cw+cw/2}" y="${T-8}" text-anchor="middle" fill="#9aa4b2">${pct(v)}</text>`);
  D.callers.forEach((c,ri)=>{
    s+=`<text x="${L-8}" y="${T+ri*rh+rh/2+3}" text-anchor="end" fill="${col(c)}">${c}</text>`;
    const cl=cellsOf(sel.arm,sel.variant,sel.depth,c);const byv={};cl.forEach(d=>byv[d.vaf]=d);
    vafs.forEach((v,ci)=>{const d=byv[v];const x=L+ci*cw,y=T+ri*rh;
      s+=`<rect class="dot" x="${x}" y="${y}" width="${cw-2}" height="${rh-3}" rx="3" fill="${grn(d?d.rate:null)}" data-t="${encodeURIComponent(d?`<b>${c}</b> @ ${pct(v)}%<br>${(100*d.rate).toFixed(0)}% (${d.k}/${d.n})`:'no data')}"/>`;
      if(d&&d.rate!=null)s+=`<text x="${x+cw/2-1}" y="${y+rh/2+1}" text-anchor="middle" fill="#0c2a1c" font-size="10">${Math.round(100*d.rate)}</text>`;});
  });
  s+=`</svg>`;const el=document.getElementById('heat');el.innerHTML=s;
  el.querySelectorAll('.dot').forEach(d=>{d.onmousemove=e=>showTip(e,decodeURIComponent(d.dataset.t));d.onmouseleave=hideTip;});
}

function drawSummary(){
  let h=`<thead><tr><th>caller</th><th class="num">empirical transition</th><th class="num">reliable (≥90%)</th><th class="num">LOD50</th><th class="num">LOD95</th><th class="num">LOD95 95% CI</th></tr></thead><tbody>`;
  D.callers.forEach(c=>{const f=fitOf(sel.arm,sel.variant,sel.depth,c);
    const lbl=c=='mitohpc'?`<span class="tag ref">${c}</span>`:c;
    const ci=(f.lod95_lo!=null&&f.lod95_hi!=null)?`${pct(f.lod95_lo)}–${pct(f.lod95_hi)}%`:'—';
    h+=`<tr><td>${lbl}</td><td class="num">${pct(f.emp_transition)}%</td><td class="num">${pct(f.emp_reliable)}%</td>`+
       `<td class="num">${f.lod50!=null?pct(f.lod50)+'%':'—'}</td><td class="num">${f.lod95!=null?pct(f.lod95)+'%'+(f.near_separable?' ★':''):'—'}</td><td class="num">${ci}</td></tr>`;});
  h+=`</tbody>`;document.getElementById('summary').innerHTML=h;
}

function drawRuntime(){
  const items=D.callers.map(c=>({c,b:rtbox(sel.arm,'all',c)})).filter(x=>x.b!=null);
  const el=document.getElementById('runtime'), tbl=document.getElementById('runtime-table');
  if(!items.length){el.innerHTML='<span class="legend">no runtime data for this arm</span>';tbl.innerHTML='';return;}
  const W=820,rh=34,L=120,R=70,T=10,B=26,H=T+B+items.length*rh,hb=9;
  let mx=0; items.forEach(d=>{mx=Math.max(mx,d.b.whi,...(d.b.outliers||[]));}); mx=Math.max(1,mx);
  const X=v=>L+(W-L-R)*(v/mx);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gx=X(mx*i/4);s+=`<line class="gl" x1="${gx}" y1="${T}" x2="${gx}" y2="${H-B}"/>`;
    s+=`<text x="${gx}" y="${H-9}" text-anchor="middle" fill="#9aa4b2">${(mx*i/4).toFixed(mx<10?1:0)}</text>`;}
  s+=`<text x="${(L+W-R)/2}" y="${H-0}" text-anchor="middle" fill="#9aa4b2">seconds / cell</text>`;
  items.forEach((d,i)=>{const b=d.b,y=T+i*rh+rh/2-2,c=col(d.c);
    s+=`<text x="${L-8}" y="${y+3}" text-anchor="end" fill="${c}">${d.c}</text>`;
    s+=`<line x1="${X(b.wlo)}" y1="${y}" x2="${X(b.whi)}" y2="${y}" stroke="${c}" stroke-width="1" opacity=".6"/>`;
    s+=`<line x1="${X(b.wlo)}" y1="${y-5}" x2="${X(b.wlo)}" y2="${y+5}" stroke="${c}"/>`;
    s+=`<line x1="${X(b.whi)}" y1="${y-5}" x2="${X(b.whi)}" y2="${y+5}" stroke="${c}"/>`;
    s+=`<rect class="dot" x="${X(b.q1)}" y="${y-hb}" width="${Math.max(1,X(b.q3)-X(b.q1))}" height="${2*hb}" rx="2" fill="${c}" fill-opacity=".28" stroke="${c}" data-t="${encodeURIComponent(`<b>${d.c}</b> (n=${b.n})<br>median ${b.med}s · mean ${b.mean}s<br>IQR ${b.q1}–${b.q3}s · range ${b.min}–${b.max}s`)}"/>`;
    s+=`<line x1="${X(b.med)}" y1="${y-hb}" x2="${X(b.med)}" y2="${y+hb}" stroke="${c}" stroke-width="2"/>`;
    const mp=X(b.mean);s+=`<path d="M${mp},${y-5} L${mp+5},${y} L${mp},${y+5} L${mp-5},${y} Z" fill="#0f1117" stroke="${c}" stroke-width="1.5"/>`;
    (b.outliers||[]).forEach(o=>{s+=`<circle cx="${X(o)}" cy="${y}" r="2.3" fill="${c}" fill-opacity=".7"/>`;});
  });
  s+=`</svg>`;el.innerHTML=s;
  el.querySelectorAll('.dot').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
  // per-depth median table (selected arm): caller × depth, fastest overall first
  const order=items.slice().sort((a,b)=>a.b.med-b.b.med).map(x=>x.c);
  let h=`<thead><tr><th>caller</th>`+D.depths.map(d=>`<th class="num">${d}×</th>`).join('')+`<th class="num">all</th></tr></thead><tbody>`;
  order.forEach(c=>{const lbl=c=='mitohpc'?`<span class="tag ref">${c}</span>`:c;
    h+=`<tr><td>${lbl}</td>`+D.depths.map(d=>{const b=rtbox(sel.arm,d,c);
      return `<td class="num">${b?b.med+'s':'—'}</td>`;}).join('')+
      (()=>{const b=rtbox(sel.arm,'all',c);return `<td class="num">${b?b.med+'s':'—'}</td>`;})()+`</tr>`;});
  h+=`</tbody>`;tbl.innerHTML=h;
}

function drawRuntimeDepth(){
  const el=document.getElementById('runtime-depth');
  const dps=D.depths; if(!dps.length){el.innerHTML='';return;}
  const W=820,H=300,L=56,R=110,T=12,B=40;
  const xs=dps; const X=i=>L+(W-L-R)*(dps.length>1?i/(dps.length-1):0.5);
  let mx=0; D.callers.forEach(c=>dps.forEach(d=>{const b=rtbox(sel.arm,d,c);if(b)mx=Math.max(mx,b.med);})); mx=Math.max(1,mx);
  const Y=v=>H-B-(H-T-B)*(v/mx);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gy=Y(mx*i/4);s+=`<line class="gl" x1="${L}" y1="${gy}" x2="${W-R}" y2="${gy}"/>`;
    s+=`<text x="${L-6}" y="${gy+3}" text-anchor="end" fill="#9aa4b2">${Math.round(mx*i/4)}s</text>`;}
  dps.forEach((d,i)=>s+=`<text x="${X(i)}" y="${H-8}" text-anchor="middle" fill="#9aa4b2">${d}×</text>`);
  s+=`<text x="${(L+W-R)/2}" y="${H-0}" text-anchor="middle" fill="#9aa4b2">simulated depth</text>`;
  D.callers.forEach(c=>{
    const pts=dps.map((d,i)=>({i,b:rtbox(sel.arm,d,c)})).filter(p=>p.b);
    if(!pts.length)return; const cc=col(c);
    s+=`<path d="${pts.map((p,j)=>`${j?'L':'M'}${X(p.i)},${Y(p.b.med)}`).join(' ')}" fill="none" stroke="${cc}" stroke-width="2"/>`;
    pts.forEach(p=>s+=`<circle class="dot" cx="${X(p.i)}" cy="${Y(p.b.med)}" r="3" fill="${cc}" data-t="${encodeURIComponent(`<b>${c}</b> @ ${dps[p.i]}×<br>median ${p.b.med}s (n=${p.b.n})`)}"/>`);
    const last=pts[pts.length-1]; s+=`<text x="${X(last.i)+6}" y="${Y(last.b.med)+3}" fill="${cc}">${c}</text>`;
  });
  s+=`</svg>`;el.innerHTML=s;
  el.querySelectorAll('.dot').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
}

function drawArmCmp(){
  let h=`<thead><tr><th>caller</th>`+D.arms.map(a=>`<th class="num">${a} LOD95</th>`).join('')+`</tr></thead><tbody>`;
  D.callers.forEach(c=>{h+=`<tr><td>${c}</td>`+D.arms.map(a=>{const f=fitOf(a,sel.variant,sel.depth,c);
    return `<td class="num">${f.lod95!=null?pct(f.lod95)+'%'+(f.near_separable?'★':''):'—'}</td>`;}).join('')+`</tr>`;});
  h+=`</tbody>`;document.getElementById('armcmp').innerHTML=h;
}

redraw();
document.getElementById('foot').innerHTML=
  `Generated by <code>pipeline/lod/make_lod_report.py</code> from a ${P.data.n_rows}-row sweep. `+
  `Simulator: MitoHPC make_testdata; detection within 30&nbsp;bp summed breakpoint error; `+
  `LOD via Firth-penalized logistic + empirical transition (Wilson 95% CIs). `+
  `Reference caller: <span class="tag ref">mitohpc</span>.`;
</script></body></html>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--fits", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--arms", default="pipeline,circular")
    ap.add_argument("--generated", default="")
    ap.add_argument("--image", default="")
    args = ap.parse_args(argv)

    cells = load_tsv(args.cells)
    fits = load_tsv(args.fits)
    sweep = load_tsv(args.sweep)
    arms = [a for a in args.arms.split(",") if a]
    data = build(cells, fits, sweep, arms)
    meta = {"generated": args.generated or "(unstamped)", "image": args.image}
    html = render(data, meta)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(html)
    sys.stderr.write("[make_lod_report] wrote %s (%d cells, %d fits, %d callers)\n"
                     % (args.out, len(cells), len(fits), len(data["callers"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
