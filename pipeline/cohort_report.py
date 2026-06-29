#!/usr/bin/env python3
"""Build a self-contained, interactive cohort SV summary (cohort_sv_summary.html).

For REAL cohorts (no per-sample truth): global, grant-reviewer-oriented metrics
that give an intuitive sense of how the callers differ. Deliberately NOT bloated —
no per-call table, no samplot — just distributions and counts:

  * KNOWN vs NOVEL per caller — how many calls match a previously-reported
    MitoBreak breakpoint (a credibility/specificity proxy), incl. the del4977
    common deletion;
  * box plots of SVs called per individual, one box per caller (all / PASS);
  * called vs passed per caller;
  * mean calls per individual by SV type, per caller;
  * heteroplasmy (VAF) distribution per caller.

Driven straight off the per-sample caller outputs (parsers.py), annotated with
the MitoBreak catalogue (mitobreak.py). Pure stdlib; no CDN.

Usage:
  cohort_report.py --root <output-root> --out cohort_sv_summary.html \\
      [--mitobreak assets/mitobreak.tsv.gz] [--samples "s1 s2 ..."] \\
      [--generated DATE] [--image URI]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from statistics import median

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import parsers  # noqa: E402
import mitobreak as MB  # noqa: E402

# directory names under the output root that are NOT samples
_NON_SAMPLE = {"by_caller", "prepared", "logs", "shards", "work", "progress",
               "manifests", "lod_report", "refcache"}

# What counts as a structural-variant CALL for the cohort comparison: a typed SV
# carrying both breakpoints (a 5'-3' span). This deliberately excludes raw evidence
# rows like MitoSeek's discordant "breakpoint" records (bp3 absent), which are not
# curated calls and would otherwise dwarf every other caller by 1000x.
SV_CALL_TYPES = ("deletion", "duplication", "inversion")


def _is_call(r):
    return (r.get("sv_type") in SV_CALL_TYPES
            and r.get("bp5") is not None and r.get("bp3") is not None)


def _passed(rec):
    """A call counts as PASS unless the caller has an explicit FILTER that isn't
    PASS. Only MitoHPC is detect-and-flag (filter in `extra`); the other callers'
    output IS their final call set, so they pass by construction."""
    if rec["caller"] == "mitohpc":
        return "filter=PASS" in (rec.get("extra") or "")
    return True


def discover(root):
    """Return (pairs, present_callers) over BOTH layouts: the classic per-sample
    root/<sample>/ (all callers) and the per-caller cascade's
    root/by_caller/<caller>/<sample>/ (one caller each). present_callers = every
    caller that has an output directory somewhere (so a caller that RAN but called
    nothing — or failed — still shows in the report, the reference especially)."""
    pairs = []
    present = set()
    bc_samples = set()
    bycaller = os.path.join(root, "by_caller")
    if os.path.isdir(bycaller):
        for caller in sorted(os.listdir(bycaller)):
            cdir = os.path.join(bycaller, caller)
            if not os.path.isdir(cdir):
                continue
            if caller in parsers.CALLERS:
                present.add(caller)
            for s in sorted(os.listdir(cdir)):
                sd = os.path.join(cdir, s)
                if os.path.isdir(sd):
                    pairs.append((s, sd)); bc_samples.add(s)
    for name in sorted(os.listdir(root)):
        # skip non-sample dirs AND samples already covered by by_caller, else a
        # leftover classic dir from a prior run would be parsed too and double-count.
        if name in _NON_SAMPLE or name in bc_samples:
            continue
        sd = os.path.join(root, name)
        sub = [c for c in parsers.CALLERS if os.path.isdir(os.path.join(sd, c))]
        if os.path.isdir(sd) and (sub or os.path.isfile(os.path.join(sd, "status.tsv"))):
            pairs.append((name, sd))
            present.update(sub)
    return pairs, present


# ---- stats (five-number summary + Tukey whiskers) ------------------------- #
def _pctile(xs, q):
    if not xs:
        return None
    if len(xs) == 1:
        return xs[0]
    pos = (q / 100.0) * (len(xs) - 1)
    lo = int(pos); frac = pos - lo
    return xs[lo] * (1 - frac) + xs[min(lo + 1, len(xs) - 1)] * frac


def boxstats(vals, nd=2):
    xs = sorted(v for v in vals if v is not None)
    n = len(xs)
    if n == 0:
        return None
    q1, med, q3 = _pctile(xs, 25), _pctile(xs, 50), _pctile(xs, 75)
    iqr = q3 - q1
    lof, hif = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inside = [x for x in xs if lof <= x <= hif]
    r = lambda v: round(v, nd)
    return {"n": n, "min": r(xs[0]), "q1": r(q1), "med": r(med), "q3": r(q3),
            "max": r(xs[-1]), "mean": r(sum(xs) / n),
            "wlo": r(min(inside) if inside else xs[0]),
            "whi": r(max(inside) if inside else xs[-1]),
            "outliers": [r(x) for x in xs if x < lof or x > hif][:40]}


def build(records, samples, db, present=None):
    # count only curated SV calls (typed, both breakpoints); drop raw evidence rows
    calls = [r for r in records if _is_call(r)]
    types = [t for t in SV_CALL_TYPES if any(r["sv_type"] == t for r in calls)]
    # callers to show: those that ran (present) or produced a call, in catalogue order
    have = set(present or ()) | {r["caller"] for r in calls}
    callers = [c for c in parsers.CALLERS if c in have]

    # annotate each call with passed + known(MitoBreak) + common-deletion once
    for r in calls:
        r["_pass"] = _passed(r)
        hit = MB.match(db, r.get("sv_type"), r.get("bp5"), r.get("bp3"))
        r["_mid"] = hit[0] if hit else None
        r["_common"] = parsers.is_common_deletion(r.get("bp5"), r.get("bp3"))

    per = {}
    for c in callers:
        crecs = [r for r in calls if r["caller"] == c]
        passed = [r for r in crecs if r["_pass"]]
        # per-individual counts (0 for samples this caller called nothing on)
        n_all = {s: 0 for s in samples}
        n_pass = {s: 0 for s in samples}
        by_type = {t: {s: 0 for s in samples} for t in types}
        for r in crecs:
            s = r["sample"]
            n_all[s] = n_all.get(s, 0) + 1
            if r["_pass"]:
                n_pass[s] = n_pass.get(s, 0) + 1
                by_type[r.get("sv_type") or "other"][s] = \
                    by_type.setdefault(r.get("sv_type") or "other", {}).get(s, 0) + 1
        known = [r for r in passed if r["_mid"]]
        novel = [r for r in passed if not r["_mid"]]
        common_samples = sorted({r["sample"] for r in passed if r["_common"]})
        common_vaf = [r["het"] for r in passed if r["_common"] and r.get("het") is not None]
        vafs = [r["het"] for r in passed if r.get("het") is not None]
        nsmp = max(1, len(samples))
        per[c] = {
            "n_called": len(crecs),
            "n_passed": len(passed),
            "known": len(known),
            "novel": len(novel),
            "known_pct": round(100.0 * len(known) / len(passed), 1) if passed else None,
            "distinct_known": len({r["_mid"] for r in known}),
            "common_samples": len(common_samples),
            "common_med_vaf": round(median(common_vaf), 4) if common_vaf else None,
            "box_all": boxstats([n_all[s] for s in samples]),
            "box_pass": boxstats([n_pass[s] for s in samples]),
            "by_type_mean": {t: round(sum(by_type[t].values()) / nsmp, 3) for t in types},
            "vaf_box": boxstats(vafs, nd=4),
        }
    callers = [c for c in callers if c in per]
    return {"samples": samples, "callers": callers, "types": types,
            "per": per, "n_samples": len(samples)}


def render_html(data, meta):
    payload = json.dumps({"data": data, "meta": meta}).replace("</", "<\\/")
    return _HTML.replace("/*__PAYLOAD__*/", payload)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="output root with per-sample (or by_caller) dirs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mitobreak", default="", help="MitoBreak DB tsv(.gz); default = vendored")
    ap.add_argument("--samples", default="", help="explicit sample list (else discovered)")
    ap.add_argument("--generated", default="")
    ap.add_argument("--image", default="")
    args = ap.parse_args(argv)

    root = os.path.abspath(args.root)
    db = MB.load(args.mitobreak or None)
    pairs, present = discover(root)
    records = []
    seen = []
    for name, sd in pairs:
        try:
            records.extend(parsers.parse_sample_dir(sd, name, pass_only=False))
        except Exception as e:  # noqa: BLE001 — an in-progress/partial dir must not sink the report
            sys.stderr.write("[cohort_report] WARN %s skipped (%r)\n" % (name, e))
        if name not in seen:
            seen.append(name)
    samples = [s for s in args.samples.replace(",", " ").split() if s] or sorted(seen)

    data = build(records, samples, db, present)
    meta = {"generated": args.generated or "(unstamped)", "image": args.image,
            "n_samples": len(samples), "mitobreak_tol": MB.MITOBREAK_TOL,
            "mitobreak_loaded": bool(db),
            "common_bp5": parsers.COMMON_DEL_BP5, "common_bp3": parsers.COMMON_DEL_BP3}
    html = render_html(data, meta)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(html)
    sys.stderr.write("[cohort_report] wrote %s (%d samples, %d callers, %d calls, MitoBreak=%s)\n"
                     % (args.out, len(samples), len(data["callers"]),
                        len(records), "yes" if db else "MISSING"))
    return 0


# --------------------------------------------------------------------------- #
_HTML = r"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>mtDNA SV cohort summary</title>
<style>
:root{--bg:#0f1117;--panel:#171a23;--panel2:#1f2430;--ink:#e6e9ef;--mut:#9aa4b2;
--line:#2a2f3a;--accent:#6ea8fe;--good:#2fbf71;--known:#6ea8fe;--novel:#e3a008;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1100px;margin:0 auto;padding:28px 20px 80px}
h1{font-size:24px;margin:0 0 4px} h2{font-size:18px;margin:34px 0 8px;border-bottom:1px solid var(--line);padding-bottom:6px}
.sub{color:var(--mut);margin:0 0 18px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:18px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.card .k{color:var(--mut);font-size:12px;text-transform:uppercase;letter-spacing:.04em}
.card .v{font-size:22px;font-weight:600;margin-top:4px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:10px 0}
.legend{color:var(--mut);font-size:12px;margin:6px 0 2px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 8px;text-align:left;border-bottom:1px solid var(--line);white-space:nowrap}
th{color:var(--mut);font-weight:600}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.tag{display:inline-block;padding:1px 7px;border-radius:999px;font-size:11px;border:1px solid var(--line)}
.tag.ref{color:#b692f6;border-color:#b692f6}
.controls{margin:8px 0;display:flex;gap:8px;flex-wrap:wrap}
button.f{background:var(--panel2);color:var(--ink);border:1px solid var(--line);
border-radius:8px;padding:5px 10px;cursor:pointer;font-size:12px}
button.f.on{border-color:var(--accent);color:var(--accent)}
svg text{fill:var(--ink);font-size:11px}
svg .gl{stroke:var(--line);stroke-dasharray:2 3;opacity:.5}
.bar{cursor:pointer}.bar:hover{opacity:.85}
#tip{position:fixed;pointer-events:none;background:#000d;border:1px solid var(--line);
border-radius:8px;padding:6px 9px;font-size:12px;color:#fff;opacity:0;transition:opacity .1s;z-index:9}
.foot{color:var(--mut);font-size:12px;margin-top:30px;border-top:1px solid var(--line);padding-top:14px}
.warn{background:rgba(229,83,75,.12);border:1px solid #6b2b27;color:#e5534b;border-radius:8px;padding:6px 10px;font-size:12px;margin:8px 0}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:middle}
</style></head><body>
<div class="wrap">
  <h1>Mitochondrial SV cohort summary</h1>
  <p class="sub" id="subtitle"></p>
  <div id="warn-slot"></div>
  <div class="cards" id="cards"></div>

  <h2>Known-breakpoint fraction per caller</h2>
  <p class="legend">Of each caller's PASS calls, the % whose breakpoint is catalogued in
    <b>MitoBreak</b> (a previously-reported mtDNA breakpoint, matched within <b id="mbtol">20</b>&nbsp;bp
    on both ends). Scale-free, so it compares callers fairly regardless of how many calls each makes:
    a high fraction means calls that align with the literature; a low fraction flags a noisy call set.
    Absolute known / novel / called / passed counts are in the table.</p>
  <div class="panel"><div id="chart-known"></div></div>
  <div class="panel" style="overflow:auto"><table id="known-table"></table></div>

  <h2>Structural variants called per individual</h2>
  <p class="legend">Distribution across the cohort of the number of calls per individual, one box per
    caller (box = IQR, line = median, diamond = mean, whiskers = 1.5×IQR, dots = outliers). Log axis —
    callers differ by orders of magnitude. Toggle every call vs PASS-only; samples a caller called
    nothing on count as 0.</p>
  <div class="controls" id="box-toggle"></div>
  <div class="panel"><div id="chart-box"></div></div>

  <h2>Calls per individual, by SV type</h2>
  <p class="legend">Left: the SV-type <b>composition</b> of each caller's PASS calls (scale-free).
    Right table: the absolute mean number of PASS calls per individual, by type.</p>
  <div class="panel"><div id="chart-type"></div><div id="type-legend" class="legend"></div></div>
  <div class="panel" style="overflow:auto"><table id="type-table"></table></div>

  <h2>Common deletion (del4977) detection</h2>
  <p class="legend">Individuals in which each caller makes a PASS call at the canonical
    <b id="cdcoord"></b> common deletion, with the median reported heteroplasmy.</p>
  <div class="panel"><div id="chart-common"></div></div>

  <h2>Reported heteroplasmy (VAF) per caller</h2>
  <p class="legend">Distribution of the heteroplasmy each caller reports for its PASS calls — callers
    differ in how they estimate VAF (coverage dosage vs junction fraction).</p>
  <div class="panel"><div id="chart-vaf"></div></div>

  <div class="foot" id="foot"></div>
</div>
<div id="tip"></div>
<script>
const P=/*__PAYLOAD__*/;const D=P.data,M=P.meta;
const CAL=D.callers, PER=D.per, S=D.samples, TYPES=D.types;
const COLOR={mitohpc:'#b692f6',eklipse:'#6ea8fe',mitosalt:'#2fbf71',
  splicebreak2:'#e3a008',mitomut:'#f178b6',mitoseek:'#4dd0e1'};
const TCOL={deletion:'#6ea8fe',duplication:'#2fbf71',inversion:'#e3a008',other:'#9aa4b2'};
const col=c=>COLOR[c]||'#9aa4b2';
const tip=document.getElementById('tip');
function showTip(e,h){tip.innerHTML=h;tip.style.opacity=1;tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY+12)+'px';}
function hideTip(){tip.style.opacity=0;}
function lbl(c){return c==='mitohpc'?`<span class="tag ref">${c}</span>`:c;}

document.getElementById('subtitle').textContent =
  `${M.n_samples} individuals · ${CAL.length} callers · generated ${M.generated}`+(M.image?` · ${M.image}`:'');
document.getElementById('mbtol').textContent=M.mitobreak_tol;
document.getElementById('cdcoord').textContent=`m.${M.common_bp5}–${M.common_bp3}`;
if(!M.mitobreak_loaded) document.getElementById('warn-slot').innerHTML=
  '<div class="warn">MitoBreak catalogue not found — known/novel split unavailable (all calls shown as novel).</div>';

// ---- cards ----
const totPass=CAL.reduce((a,c)=>a+PER[c].n_passed,0);
const totKnown=CAL.reduce((a,c)=>a+PER[c].known,0);
const commonCallers=CAL.filter(c=>PER[c].common_samples>0).length;
const mostCommon=CAL.slice().sort((a,b)=>PER[b].common_samples-PER[a].common_samples)[0];
const cards=[
  ['Individuals',M.n_samples],['Callers',CAL.length],
  ['PASS calls (cohort)',totPass],
  ['Known (MitoBreak)',totPass?`${totKnown} · ${Math.round(100*totKnown/totPass)}%`:'–'],
  ['Detect common del',`${commonCallers}/${CAL.length} callers`],
  ['Top common-del caller',mostCommon?`${mostCommon} · ${PER[mostCommon].common_samples}`:'–'],
];
document.getElementById('cards').innerHTML=cards.map(([k,v])=>
  `<div class="card"><div class="k">${k}</div><div class="v">${v}</div></div>`).join('');

// ---- generic stacked horizontal bar (segments=[{v,color,name}]) ----
function stackChart(id, rows, unit, xlabel){
  const W=860,rowH=30,padL=120,padR=70,padT=10,padB=28;
  const H=padT+padB+rows.length*rowH;
  const max=Math.max(1,...rows.map(r=>r.segs.reduce((a,s)=>a+s.v,0)));
  const X=v=>padL+(W-padL-padR)*(v/max);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gx=X(max*i/4);s+=`<line class="gl" x1="${gx}" y1="${padT}" x2="${gx}" y2="${H-padB}"/>`;
    s+=`<text x="${gx}" y="${H-9}" text-anchor="middle" fill="#9aa4b2">${(max*i/4).toFixed(max<10?1:0)}</text>`;}
  if(xlabel)s+=`<text x="${(padL+W-padR)/2}" y="${H-0}" text-anchor="middle" fill="#9aa4b2">${xlabel}</text>`;
  rows.forEach((r,i)=>{let x=padL,y=padT+i*rowH+4;
    s+=`<text x="${padL-8}" y="${y+13}" text-anchor="end">${r.label}</text>`;
    r.segs.forEach(seg=>{const w=X(x===padL?seg.v:seg.v)-padL+0;const ww=(W-padL-padR)*(seg.v/max);
      if(seg.v>0){s+=`<rect class="bar" x="${x}" y="${y}" width="${Math.max(0,ww)}" height="${rowH-12}" fill="${seg.color}" data-t="${encodeURIComponent(`<b>${r.label}</b><br>${seg.name}: ${seg.v}${unit||''}`)}"/>`;}
      x+=ww;});
    const tot=r.segs.reduce((a,sg)=>a+sg.v,0);
    s+=`<text x="${x+6}" y="${y+13}">${r.total!=null?r.total:tot}${unit||''}</text>`;});
  s+=`</svg>`;const el=document.getElementById(id);el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
}

// ---- horizontal percent bar (0-100%); rows=[{label,v,color,tip}] ----
function pctBar(id, rows){
  const W=860,rowH=28,padL=120,padR=60,padT=8,padB=24,H=padT+padB+rows.length*rowH;
  const X=v=>padL+(W-padL-padR)*(Math.max(0,Math.min(100,v))/100);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  for(let i=0;i<=4;i++){const gx=X(100*i/4);s+=`<line class="gl" x1="${gx}" y1="${padT}" x2="${gx}" y2="${H-padB}"/>`;
    s+=`<text x="${gx}" y="${H-8}" text-anchor="middle" fill="#9aa4b2">${100*i/4}%</text>`;}
  rows.forEach((r,i)=>{const y=padT+i*rowH+4;
    s+=`<text x="${padL-8}" y="${y+12}" text-anchor="end">${r.label}</text>`;
    if(r.v!=null){s+=`<rect class="bar" x="${padL}" y="${y}" width="${Math.max(0,X(r.v)-padL)}" height="${rowH-12}" rx="3" fill="${r.color}" data-t="${encodeURIComponent(r.tip)}"/>`;
      s+=`<text x="${X(r.v)+6}" y="${y+12}">${r.v}%</text>`;}
    else s+=`<text x="${padL+4}" y="${y+12}" fill="#9aa4b2">no PASS calls</text>`;});
  s+=`</svg>`;const el=document.getElementById(id);el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
}

// ---- horizontal box-and-whisker (items=[{label,color,b}]); log: log1p X-axis ----
function boxChart(id, items, xlabel, log){
  items=items.filter(d=>d.b);const el=document.getElementById(id);
  if(!items.length){el.innerHTML='<span class="legend">no data</span>';return;}
  const W=860,rowH=34,padL=120,padR=70,padT=10,padB=26,hb=9;
  const H=padT+padB+items.length*rowH;
  let mx=0;items.forEach(d=>{mx=Math.max(mx,d.b.whi,...(d.b.outliers||[]));});mx=Math.max(mx||1,0.0001);
  const tx=v=>log?Math.log1p(Math.max(0,v)):v, tmx=tx(mx);
  const X=v=>padL+(W-padL-padR)*(tmx>0?tx(v)/tmx:0);
  // gridline tick values: log -> nice 1/3/10/... ladder; linear -> quarters
  let ticks=[];
  if(log){const cand=[0,1,3,10,30,100,300,1000,3000,10000,30000];
    ticks=cand.filter(v=>v<=mx); if(ticks[ticks.length-1]!==mx)ticks.push(Math.round(mx));}
  else for(let i=0;i<=4;i++)ticks.push(mx*i/4);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  ticks.forEach(t=>{const gx=X(t);s+=`<line class="gl" x1="${gx}" y1="${padT}" x2="${gx}" y2="${H-padB}"/>`;
    s+=`<text x="${gx}" y="${H-9}" text-anchor="middle" fill="#9aa4b2">${t<10?(+t).toFixed(t%1?1:0):Math.round(t)}</text>`;});
  s+=`<text x="${(padL+W-padR)/2}" y="${H-0}" text-anchor="middle" fill="#9aa4b2">${xlabel||''}${log?' (log)':''}</text>`;
  items.forEach((d,i)=>{const b=d.b,y=padT+i*rowH+rowH/2-2,c=d.color;
    s+=`<text x="${padL-8}" y="${y+3}" text-anchor="end">${d.label}</text>`;
    s+=`<line x1="${X(b.wlo)}" y1="${y}" x2="${X(b.whi)}" y2="${y}" stroke="${c}" opacity=".6"/>`;
    s+=`<line x1="${X(b.wlo)}" y1="${y-5}" x2="${X(b.wlo)}" y2="${y+5}" stroke="${c}"/>`;
    s+=`<line x1="${X(b.whi)}" y1="${y-5}" x2="${X(b.whi)}" y2="${y+5}" stroke="${c}"/>`;
    s+=`<rect class="bar" x="${X(b.q1)}" y="${y-hb}" width="${Math.max(1,X(b.q3)-X(b.q1))}" height="${2*hb}" rx="2" fill="${c}" fill-opacity=".28" stroke="${c}" data-t="${encodeURIComponent(`<b>${d.label}</b> (n=${b.n})<br>median ${b.med} · mean ${b.mean}<br>IQR ${b.q1}–${b.q3} · range ${b.min}–${b.max}`)}"/>`;
    s+=`<line x1="${X(b.med)}" y1="${y-hb}" x2="${X(b.med)}" y2="${y+hb}" stroke="${c}" stroke-width="2"/>`;
    const mp=X(b.mean);s+=`<path d="M${mp},${y-5} L${mp+5},${y} L${mp},${y+5} L${mp-5},${y} Z" fill="var(--bg)" stroke="${c}" stroke-width="1.5"/>`;
    (b.outliers||[]).forEach(o=>{s+=`<circle cx="${X(o)}" cy="${y}" r="2.3" fill="${c}" fill-opacity=".7"/>`;});});
  s+=`</svg>`;el.innerHTML=s;
  el.querySelectorAll('.bar').forEach(b=>{b.onmousemove=e=>showTip(e,decodeURIComponent(b.dataset.t));b.onmouseleave=hideTip;});
}

// ---- known-breakpoint fraction (scale-free %; sorted desc) ----
(function(){
  const order=CAL.slice().sort((a,b)=>(PER[b].known_pct??-1)-(PER[a].known_pct??-1));
  pctBar('chart-known', order.map(c=>({label:c,v:PER[c].known_pct,color:col(c),
    tip:`<b>${c}</b><br>${PER[c].known}/${PER[c].n_passed} PASS calls known<br>${PER[c].novel} novel · ${PER[c].distinct_known} distinct catalogue entries`})));
  let h=`<thead><tr><th>caller</th><th class="num">called</th><th class="num">passed</th>`
    +`<th class="num">known</th><th class="num">novel</th><th class="num">known%</th>`
    +`<th class="num">distinct known</th><th class="num">common-del indiv.</th></tr></thead><tbody>`;
  order.forEach(c=>{const p=PER[c];h+=`<tr><td>${lbl(c)}</td><td class="num">${p.n_called}</td>`
    +`<td class="num">${p.n_passed}</td><td class="num">${p.known}</td><td class="num">${p.novel}</td>`
    +`<td class="num">${p.known_pct==null?'–':p.known_pct+'%'}</td><td class="num">${p.distinct_known}</td>`
    +`<td class="num">${p.common_samples}${p.common_med_vaf!=null?` <span style="color:var(--mut)">(VAF ${p.common_med_vaf})</span>`:''}</td></tr>`;});
  document.getElementById('known-table').innerHTML=h+'</tbody>';
})();

// ---- box: calls per individual (all <-> pass) ----
let boxMode='pass';
function renderBox(){
  const key=boxMode==='pass'?'box_pass':'box_all';
  const order=CAL.slice().sort((a,b)=>((PER[b][key]||{}).med||0)-((PER[a][key]||{}).med||0));
  boxChart('chart-box', order.map(c=>({label:c,color:col(c),b:PER[c][key]})), 'calls / individual', true);
}
const bt=document.getElementById('box-toggle');
[['pass','PASS only'],['all','all calls']].forEach(([k,t])=>{const b=document.createElement('button');
  b.className='f'+(k===boxMode?' on':'');b.textContent=t;
  b.onclick=()=>{boxMode=k;bt.querySelectorAll('button').forEach(x=>x.classList.remove('on'));b.classList.add('on');renderBox();};
  bt.appendChild(b);});
renderBox();

// ---- type composition (% mix, scale-free) + absolute mean/indiv table ----
(function(){
  const sum=c=>TYPES.reduce((x,t)=>x+(PER[c].by_type_mean[t]||0),0);
  const order=CAL.slice().sort((a,b)=>sum(b)-sum(a));
  stackChart('chart-type', order.map(c=>{const s=sum(c);
    return {label:c, total:'',
      segs:TYPES.map(t=>({v:s>0?Math.round(1000*(PER[c].by_type_mean[t]||0)/s)/10:0,
        color:TCOL[t]||'#9aa4b2',name:t}))};}), '%', 'type composition of PASS calls');
  document.getElementById('type-legend').innerHTML=
    TYPES.map(t=>`<span class="sw" style="background:${TCOL[t]||'#9aa4b2'}"></span>${t}`).join(' &nbsp; ');
  let h=`<thead><tr><th>caller</th>`+TYPES.map(t=>`<th class="num">${t}</th>`).join('')
    +`<th class="num">total / indiv.</th></tr></thead><tbody>`;
  order.forEach(c=>{h+=`<tr><td>${lbl(c)}</td>`
    +TYPES.map(t=>`<td class="num">${(PER[c].by_type_mean[t]||0).toFixed(2)}</td>`).join('')
    +`<td class="num">${sum(c).toFixed(2)}</td></tr>`;});
  document.getElementById('type-table').innerHTML=h+'</tbody>';
})();

// ---- common deletion detection (bar of # individuals) ----
(function(){
  const order=CAL.slice().sort((a,b)=>PER[b].common_samples-PER[a].common_samples);
  stackChart('chart-common', order.map(c=>({label:c,total:PER[c].common_samples,
    segs:[{v:PER[c].common_samples,color:col(c),name:'individuals with del4977'}]})),
    '', 'individuals (PASS) with the common deletion');
})();

// ---- VAF distribution per caller ----
boxChart('chart-vaf', CAL.map(c=>({label:c,color:col(c),b:PER[c].vaf_box})), 'reported heteroplasmy (VAF)');

document.getElementById('foot').innerHTML=
  `Generated by <code>pipeline/cohort_report.py</code>. Global cohort metrics only — no per-call table, `
  +`no per-call visualisation. A <b>call</b> is a typed SV (deletion/duplication/inversion) carrying both `
  +`breakpoints; raw discordant-read evidence is not counted. PASS = the caller's confident set (MitoHPC `
  +`FILTER==PASS; the other callers' output is already their final call set). "Known" = breakpoint within `
  +`${M.mitobreak_tol} bp (both ends) of a MitoBreak catalogue entry. Reference: <span class="tag ref">mitohpc</span> (MitoHPC).`;
</script></body></html>
"""


if __name__ == "__main__":
    raise SystemExit(main())
