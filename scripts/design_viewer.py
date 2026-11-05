"""Inspect a design's specificity: IGV-loadable tracks plus a self-contained HTML viewer.

    python scripts/design_viewer.py --prefix design/output/panel --reference ref.fa --outdir panel_view

Writes to --outdir:
  amplicons.bed          BED12: amplicon span with both primers as thick blocks, coloured by pool
  mispriming_sites.bed   3'-anchored near-match priming sites of every primer (excluding its own site)
  unintended_products.bed  products predicted from any two facing sites (own on-target pair excluded)
  session.xml            IGV session loading the tracks on the reference (untested in IGV desktop)
  specificity.html       self-contained viewer: overview, per-target table, locus view

"Near-match" = BLAST alignment of >= --min-len bases that includes the primer 3' end, with at most
--max-mismatch mismatches and an exact last 4 bases (a 3' mismatch blocks extension). This is
informational: it flags sites that exact-hit screening in design.py cannot see.
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
import design  # noqa: E402

POOL_COLORS = ["#2563eb", "#d97706", "#059669", "#9333ea", "#dc2626", "#0891b2"]


def hex_to_rgb(h):
    return ",".join(str(int(h[i:i + 2], 16)) for i in (1, 3, 5))


def load_design(prefix):
    """All pool CSV rows as dicts with integer coordinates and primer footprints (0-based half-open)."""
    rows = []
    for path in sorted(glob.glob(f"{prefix}_pool_*.csv"), key=lambda p: int(p.rsplit("_pool_", 1)[1][:-4])):
        for r in csv.DictReader(open(path)):
            r['pool'] = int(r['pool'])
            r['start'], r['end'] = int(r['start']), int(r['end'])
            fl, rl = len(r['fwd_primer']), len(r['rev_primer'])
            if r['strand'] == '+':
                r['f'], r['r'] = (r['start'], r['start'] + fl), (r['end'] - rl, r['end'])
            else:   # sense primer is at the right-hand end of a minus-strand amplicon
                r['f'], r['r'] = (r['end'] - fl, r['end']), (r['start'], r['start'] + rl)
            rows.append(r)
    if not rows:
        raise SystemExit(f"No pool CSVs found for prefix '{prefix}'.")
    return rows


def near_matches(primers, db_prefix, min_len, max_mismatch, threads=4):
    """3'-anchored near-match sites. primers: {name: seq}. Returns list of site dicts."""
    query = "".join(f">{n}\n{s}\n" for n, s in primers.items())
    cmd = ['blastn', '-query', '-', '-db', db_prefix, '-task', 'blastn-short', '-word_size', '7',
           '-evalue', '1000', '-perc_identity', '70', '-dust', 'no', '-soft_masking', 'false',
           '-max_target_seqs', '5000', '-num_threads', str(threads),
           '-outfmt', '6 qseqid sseqid length mismatch gaps qstart qend sstart send qseq sseq']
    out = subprocess.run(cmd, input=query, capture_output=True, text=True, check=True).stdout
    sites = []
    for line in out.splitlines():
        q, contig, length, mm, gaps, qs, qe, ss, se, qseq, sseq = line.split('\t')
        length, mm, gaps, qe, ss, se = int(length), int(mm), int(gaps), int(qe), int(ss), int(se)
        if gaps or qe != len(primers[q]) or length < min_len or mm > max_mismatch:
            continue
        if qseq[-4:].upper() != sseq[-4:].upper():
            continue                                    # 3' end must match to extend
        plus = ss <= se
        sites.append({'primer': q, 'contig': contig, 'lo': min(ss, se) - 1, 'hi': max(ss, se),
                      'strand': '+' if plus else '-', 'len': length, 'mm': mm})
    return sites


def overlaps(site, contig, lo, hi):
    return site['contig'] == contig and site['lo'] < hi and site['hi'] > lo


def build(rows, sites, max_product):
    """Annotate sites as own/other; predict unintended products between facing sites in the same pool."""
    own = {}
    for r in rows:
        own[f"{r['target_id']}_F"] = (r['contig'], *r['f'])
        own[f"{r['target_id']}_R"] = (r['contig'], *r['r'])
    pool_of = {f"{r['target_id']}_{x}": r['pool'] for r in rows for x in 'FR'}
    mis = [s for s in sites if not overlaps(s, *own[s['primer']])]
    for s in mis:
        s['pool'] = pool_of[s['primer']]
        s['target'] = s['primer'].rsplit('_', 1)[0]

    # every priming site in a pool: on-target footprints + near-match sites
    allsites = []
    for r in rows:
        for x, fp in (('F', r['f']), ('R', r['r'])):
            strand = (r['strand'] if x == 'F' else ('-' if r['strand'] == '+' else '+'))
            allsites.append({'primer': f"{r['target_id']}_{x}", 'target': r['target_id'], 'pool': r['pool'],
                             'contig': r['contig'], 'lo': fp[0], 'hi': fp[1], 'strand': strand, 'own': True})
    allsites += [dict(s, own=False) for s in mis]
    products = []
    for pool in sorted({s['pool'] for s in allsites}):
        ps = [s for s in allsites if s['pool'] == pool]
        for a in (s for s in ps if s['strand'] == '+'):
            for b in (s for s in ps if s['strand'] == '-'):
                if a['contig'] != b['contig'] or not a['lo'] < b['hi']:
                    continue
                size = b['hi'] - a['lo']
                if size > max_product or (a['own'] and b['own'] and a['target'] == b['target']):
                    continue
                products.append({'pool': pool, 'contig': a['contig'], 'lo': a['lo'], 'hi': b['hi'], 'size': size,
                                 'a': a['primer'], 'b': b['primer'], 'a_own': a['own'], 'b_own': b['own']})
    return mis, sorted(products, key=lambda p: p['size'])


def write_tracks(outdir, rows, mis, products, reference_fasta):
    os.makedirs(outdir, exist_ok=True)
    with open(f"{outdir}/amplicons.bed", "w") as f:
        f.write('track name="amplicons" description="amplicon with primers as blocks" itemRgb="On"\n')
        for r in rows:
            (fl, fh), (rl, rh) = sorted([r['f'], r['r']])
            rgb = hex_to_rgb(POOL_COLORS[(r['pool'] - 1) % len(POOL_COLORS)])
            f.write(f"{r['contig']}\t{r['start']}\t{r['end']}\t{r['name']}\t0\t{r['strand']}\t{r['start']}\t{r['end']}"
                    f"\t{rgb}\t2\t{fh - fl},{rh - rl}\t0,{rl - r['start']}\n")
    with open(f"{outdir}/mispriming_sites.bed", "w") as f:
        f.write('track name="mispriming_sites" description="3-prime anchored near-match priming sites" itemRgb="On"\n')
        for s in mis:
            f.write(f"{s['contig']}\t{s['lo']}\t{s['hi']}\t{s['primer']}_{s['len']}bp_{s['mm']}mm\t{s['len'] * 50}"
                    f"\t{s['strand']}\t{s['lo']}\t{s['hi']}\t220,38,38\n")
    with open(f"{outdir}/unintended_products.bed", "w") as f:
        f.write('track name="unintended_products" description="predicted from facing priming sites" itemRgb="On"\n')
        for p in products:
            f.write(f"{p['contig']}\t{p['lo']}\t{p['hi']}\t{p['a']}__{p['b']}_{p['size']}bp\t0\t.\t{p['lo']}\t{p['hi']}\t217,119,6\n")
    first = rows[0]
    with open(f"{outdir}/session.xml", "w") as f:
        f.write(f'''<?xml version="1.0" encoding="UTF-8"?>
<Session genome="{os.path.abspath(reference_fasta)}" locus="{first['contig']}:{max(1, first['start'] - 500)}-{first['end'] + 500}" version="8">
    <Resources>
        <Resource path="amplicons.bed" type="bed"/>
        <Resource path="mispriming_sites.bed" type="bed"/>
        <Resource path="unintended_products.bed" type="bed"/>
    </Resources>
</Session>
''')


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Design specificity viewer</title>
<style>
:root{--bg:#fafafa;--fg:#1f2937;--mut:#6b7280;--card:#fff;--line:#d1d5db;--bad:#dc2626;--warn:#d97706;--ok:#059669}
@media (prefers-color-scheme:dark){:root{--bg:#111827;--fg:#e5e7eb;--mut:#9ca3af;--card:#1f2937;--line:#374151}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,sans-serif}
main{max-width:1180px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:22px 0 8px}
.sub{color:var(--mut);margin:0 0 14px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.card b{display:block;font-size:22px}.card span{color:var(--mut);font-size:12px}
.panel{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px;overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{color:var(--mut);font-weight:600;cursor:pointer;position:sticky;top:0;background:var(--card)}
tr.row{cursor:pointer}tr.row:hover{background:rgba(127,127,127,.12)}tr.sel{background:rgba(37,99,235,.15)}
.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px}
.flag{color:var(--bad);font-weight:600}.good{color:var(--ok)}
.legend span{margin-right:14px;font-size:12px;color:var(--mut)}
.sw{display:inline-block;width:14px;height:8px;margin-right:4px;vertical-align:middle}
svg text{fill:var(--fg);font:11px system-ui,sans-serif}
.note{color:var(--mut);font-size:12px;margin-top:6px}
</style></head><body><main>
<h1>Design specificity viewer</h1>
<p class="sub" id="sub"></p>
<div class="cards" id="cards"></div>

<h2>Genome overview</h2>
<div class="panel"><div class="legend" id="legend"></div><div id="overview"></div>
<div class="note">Amplicons are drawn at least 3 px wide so they stay visible. Click a row below to zoom.</div></div>

<h2>Targets <span class="note">(click a column to sort, a row to open the locus view)</span></h2>
<div class="panel" style="max-height:340px;overflow:auto"><table id="tbl"></table></div>

<h2 id="locusTitle">Locus view</h2>
<div class="panel"><div id="locus"><span class="note">Select a target.</span></div></div>
<p class="note">Near-match = BLAST alignment ≥ <span id="p_len"></span> bases including the primer 3′ end, ≤ <span id="p_mm"></span> mismatches, exact last 4 bases.
Informational: exact-hit screening in <code>design.py</code> cannot see these. Unintended products need two facing sites in the <em>same pool</em> within <span id="p_prod"></span> bp.</p>
</main>
<script>
const D = /*DATA*/;
const C = D.colors, el = id => document.getElementById(id);
const fmt = n => n.toLocaleString();
el('sub').textContent = `${D.targets.length} targets, ${D.n_pools} pool(s), reference ${D.contigs.map(c=>c.name+' ('+fmt(c.len)+' bp)').join(', ')}`;
el('p_len').textContent = D.params.min_len; el('p_mm').textContent = D.params.max_mismatch; el('p_prod').textContent = fmt(D.params.max_product);
const stats = [[D.targets.length,'targets'],[D.n_pools,'pools'],[D.targets.length*2,'primers'],
  [D.sites.length,'near-match sites'],[D.products.length,'unintended products'],
  [D.targets.filter(t=>t.nsites+t.nprod===0).length,'targets with no flags']];
el('cards').innerHTML = stats.map(([v,l])=>`<div class="card"><b>${v}</b><span>${l}</span></div>`).join('');
el('legend').innerHTML = D.pools.map(p=>`<span><i class="sw" style="background:${C[(p-1)%C.length]}"></i>pool ${p} amplicon</span>`).join('')
  + '<span><i class="sw" style="background:#dc2626"></i>near-match site</span><span><i class="sw" style="background:#d97706"></i>unintended product</span>';

function overview(){
  const W=1100, out=[];
  D.contigs.forEach(c=>{
    const x=p=>20+p/c.len*(W-40), tracks=D.pools.length+2, H=56+tracks*22;
    let s=`<svg width="${W}" height="${H}" role="img" aria-label="overview of ${c.name}"><text x="20" y="14">${c.name}</text>`;
    s+=`<line x1="20" x2="${W-20}" y1="22" y2="22" stroke="var(--line)"/>`;
    for(let i=0;i<=10;i++){const p=c.len*i/10;s+=`<line x1="${x(p)}" x2="${x(p)}" y1="19" y2="25" stroke="var(--mut)"/><text x="${x(p)-14}" y="36">${(p/1e6).toFixed(2)}M</text>`}
    D.pools.forEach((pl,i)=>{const y=46+i*22;s+=`<text x="2" y="${y+9}" style="font-size:10px">P${pl}</text>`;
      D.targets.filter(t=>t.pool===pl&&t.contig===c.name).forEach(t=>{s+=`<rect data-t="${t.name}" x="${x(t.start)}" y="${y}" width="${Math.max(3,x(t.end)-x(t.start))}" height="12" fill="${C[(pl-1)%C.length]}" style="cursor:pointer"><title>${t.target} ${fmt(t.start)}-${fmt(t.end)} (${t.size} bp)</title></rect>`})});
    let y=46+D.pools.length*22; s+=`<text x="2" y="${y+9}" style="font-size:10px">sites</text>`;
    D.sites.filter(q=>q.contig===c.name).forEach(q=>{s+=`<line x1="${x(q.lo)}" x2="${x(q.lo)}" y1="${y}" y2="${y+12}" stroke="#dc2626" stroke-width="1.5"><title>${q.primer} near-match ${q.len} bp, ${q.mm} mm at ${fmt(q.lo)}</title></line>`});
    y+=22; s+=`<text x="2" y="${y+9}" style="font-size:10px">prod.</text>`;
    D.products.filter(q=>q.contig===c.name).forEach(q=>{s+=`<rect x="${x(q.lo)}" y="${y+3}" width="${Math.max(2,x(q.hi)-x(q.lo))}" height="6" fill="#d97706"><title>${q.a} + ${q.b}: ${q.size} bp (pool ${q.pool})</title></rect>`});
    out.push(s+'</svg>');
  });
  el('overview').innerHTML=out.join('');
  el('overview').querySelectorAll('rect[data-t]').forEach(r=>r.addEventListener('click',()=>select(r.dataset.t)));
}

let sortKey='start', asc=true, current=null;
const cols=[['name','target'],['pool','pool'],['contig','contig'],['start','start'],['size','size'],['strand','strand'],
  ['tm','Tm F/R'],['hits','exact hits'],['nsites','near-match sites'],['nprod','unintended products']];
function table(){
  const rows=[...D.targets].sort((a,b)=>{const x=a[sortKey],y=b[sortKey];return (x<y?-1:x>y?1:0)*(asc?1:-1)});
  el('tbl').innerHTML='<tr>'+cols.map(([k,l])=>`<th data-k="${k}">${l}${k===sortKey?(asc?' ▲':' ▼'):''}</th>`).join('')+'</tr>'+
   rows.map(t=>`<tr class="row ${t.name===current?'sel':''}" data-t="${t.name}"><td><i class="dot" style="background:${C[(t.pool-1)%C.length]}"></i>${t.target}</td><td>${t.pool}</td><td>${t.contig}</td><td>${fmt(t.start)}</td><td>${t.size}</td><td>${t.strand}</td><td>${t.tm}</td><td>${t.hits}</td>
   <td class="${t.nsites?'flag':'good'}">${t.nsites}</td><td class="${t.nprod?'flag':'good'}">${t.nprod}</td></tr>`).join('');
  el('tbl').querySelectorAll('th').forEach(h=>h.onclick=()=>{const k=h.dataset.k;asc=(k===sortKey)?!asc:true;sortKey=k;table()});
  el('tbl').querySelectorAll('tr.row').forEach(r=>r.onclick=()=>select(r.dataset.t));
}
function select(name){current=name;table();locus(D.targets.find(t=>t.name===name))}
function arrow(x1,x2,y,col,plus,label,tip){const h=10,hd=Math.min(7,Math.abs(x2-x1));
  const pts=plus?`${x1},${y} ${x2-hd},${y} ${x2},${y+h/2} ${x2-hd},${y+h} ${x1},${y+h}`:`${x2},${y} ${x1+hd},${y} ${x1},${y+h/2} ${x1+hd},${y+h} ${x2},${y+h}`;
  return `<polygon points="${pts}" fill="${col}"><title>${tip}</title></polygon>`+(label?`<text x="${x1}" y="${y-3}" style="font-size:10px">${label}</text>`:'')}
function locus(t){
  const pad=600, lo=Math.max(0,t.start-pad), hi=t.end+pad, W=1100, x=p=>20+(p-lo)/(hi-lo)*(W-40);
  const near=D.sites.filter(q=>q.contig===t.contig&&q.hi>lo&&q.lo<hi);
  const prods=D.products.filter(q=>q.contig===t.contig&&q.hi>lo&&q.lo<hi);
  const mine=D.sites.filter(q=>q.target===t.target);
  let s=`<svg width="${W}" height="190"><text x="20" y="12">${t.contig}:${fmt(lo)}-${fmt(hi)}</text><line x1="20" x2="${W-20}" y1="30" y2="30" stroke="var(--line)"/>`;
  for(let i=0;i<=6;i++){const p=lo+(hi-lo)*i/6;s+=`<line x1="${x(p)}" x2="${x(p)}" y1="27" y2="33" stroke="var(--mut)"/><text x="${x(p)-18}" y="44">${fmt(Math.round(p))}</text>`}
  s+=`<text x="2" y="66" style="font-size:10px">amplicon</text><rect x="${x(t.start)}" y="62" width="${x(t.end)-x(t.start)}" height="4" fill="${C[(t.pool-1)%C.length]}"/>`;
  s+=arrow(x(t.f[0]),x(t.f[1]),58,'#2563eb',t.strand==='+',t.target+'_F','forward (sense) primer '+t.fwd)+arrow(x(t.r[0]),x(t.r[1]),58,'#2563eb',t.strand!=='+',t.target+'_R','reverse (antisense) primer '+t.rev);
  s+=`<text x="2" y="100" style="font-size:10px">near-match</text>`;
  near.forEach((q,i)=>{const own=q.target===t.target;s+=arrow(x(q.lo),x(q.hi),90+(i%3)*14,own?'#d97706':'#dc2626',q.strand==='+','',`${q.primer}: ${q.len} bp, ${q.mm} mismatch(es), ${q.strand} strand at ${fmt(q.lo)}`)});
  s+=`<text x="2" y="150" style="font-size:10px">products</text>`;
  prods.forEach((q,i)=>{const y=140+(i%3)*10;s+=`<rect x="${x(Math.max(q.lo,lo))}" y="${y}" width="${Math.max(2,x(Math.min(q.hi,hi))-x(Math.max(q.lo,lo)))}" height="5" fill="#d97706"><title>${q.a} + ${q.b}: ${q.size} bp</title></rect>`});
  el('locus').innerHTML=s+'</svg>'+`<div class="note">Blue arrows = designed primers (arrowhead = extension direction). Red/amber arrows = near-match priming sites in this window (amber = this target's own primers). ${near.length} site(s), ${prods.length} predicted product(s) in view.</div>`
   +(mine.length?'<h2 style="margin-top:14px">Near-match sites of this target\'s primers (anywhere in the reference)</h2><table><tr><th>primer</th><th>location</th><th>strand</th><th>aligned bp</th><th>mismatches</th></tr>'+mine.map(q=>`<tr><td>${q.primer}</td><td>${q.contig}:${fmt(q.lo+1)}-${fmt(q.hi)}</td><td>${q.strand}</td><td>${q.len}</td><td>${q.mm}</td></tr>`).join('')+'</table>':'')
   +(prods.length?'<table>'+'<tr><th>pool</th><th>primer A</th><th>primer B</th><th>size</th></tr>'+prods.slice(0,15).map(q=>`<tr><td>${q.pool}</td><td>${q.a}${q.a_own?'':' (near-match)'}</td><td>${q.b}${q.b_own?'':' (near-match)'}</td><td>${q.size} bp</td></tr>`).join('')+'</table>':'');
  el('locusTitle').textContent='Locus view: '+t.target;
}
overview();table();if(D.targets.length)select(D.targets[0].name);
</script></body></html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--prefix", required=True, help="design --output-prefix")
    ap.add_argument("--reference", required=True, nargs="+", help="reference FASTA(s) to search for near-matches")
    ap.add_argument("--blast-dir", default="ncbi_data")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-len", type=int, default=16, help="min aligned bases (default 16)")
    ap.add_argument("--max-mismatch", type=int, default=1)
    ap.add_argument("--max-product", type=int, default=1000, help="products up to this size are flagged (match design.py --max-cross-product)")
    a = ap.parse_args()

    rows = load_design(a.prefix)
    primers = {}
    for r in rows:
        primers[f"{r['target_id']}_F"], primers[f"{r['target_id']}_R"] = r['fwd_primer'], r['rev_primer']
    sites = []
    for fa in a.reference:
        prefix = os.path.join(a.blast_dir, os.path.splitext(os.path.basename(fa))[0])
        design.ensure_blast_db(fa, prefix)
        sites += near_matches(primers, prefix, a.min_len, a.max_mismatch)
    mis, products = build(rows, sites, a.max_product)
    write_tracks(a.outdir, rows, mis, products, a.reference[0])

    contigs = {}
    for fa in a.reference:
        for line in open(fa):
            if line.startswith(">"):
                cur = line[1:].split()[0]
                contigs[cur] = 0
            else:
                contigs[cur] += len(line.strip())
    used = sorted({r['contig'] for r in rows} | {s['contig'] for s in mis})
    targets = []
    for r in rows:
        tid = r['target_id']
        ns = sum(1 for s in mis if s['target'] == tid)
        npd = sum(1 for p in products if tid in (p['a'].rsplit('_', 1)[0], p['b'].rsplit('_', 1)[0]))
        targets.append({'name': r['name'], 'target': tid, 'pool': r['pool'], 'contig': r['contig'], 'start': r['start'],
                        'end': r['end'], 'size': int(r['amplicon_size']), 'strand': r['strand'],
                        'tm': f"{float(r['fwd_tm']):.1f}/{float(r['rev_tm']):.1f}", 'hits': r['specificity_hits'],
                        'fwd': r['fwd_primer'], 'rev': r['rev_primer'], 'f': list(r['f']), 'r': list(r['r']),
                        'nsites': ns, 'nprod': npd})
    data = {'targets': targets, 'sites': mis, 'products': products, 'colors': POOL_COLORS,
            'pools': sorted({t['pool'] for t in targets}), 'n_pools': len({t['pool'] for t in targets}),
            'contigs': [{'name': c, 'len': contigs[c]} for c in used if c in contigs],
            'params': {'min_len': a.min_len, 'max_mismatch': a.max_mismatch, 'max_product': a.max_product}}
    with open(os.path.join(a.outdir, "specificity.html"), "w") as f:
        f.write(HTML.replace("/*DATA*/", json.dumps(data).replace("</", "<\\/")))
    flagged = sum(1 for t in targets if t['nsites'] or t['nprod'])
    print(f"{len(targets)} targets in {data['n_pools']} pool(s): {len(mis)} near-match sites, "
          f"{len(products)} unintended predicted products, {flagged} targets flagged")
    print(f"Wrote {a.outdir}/specificity.html, amplicons.bed, mispriming_sites.bed, unintended_products.bed, session.xml")


if __name__ == "__main__":
    main()
