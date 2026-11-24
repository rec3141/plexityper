// Headless check of the Plexi Pond simulation core (extracted from index.html).
// Usage: node headless_test.js [ticks] [seeds]
'use strict';
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8');
const m = html.match(/<script id="core">([\s\S]*?)<\/script>/);
if (!m) throw new Error('core script not found');
const C = new Function('module', m[1] + '\nreturn PlexiCore;')(undefined);
const T = +(process.argv[2] || 3000), nSeeds = +(process.argv[3] || 2);
let fails = 0;
function ok(cond, msg){ console.log((cond ? 'PASS ' : 'FAIL ') + msg); if (!cond) fails++; }

// 1. unit checks of the pure graph logic
ok(C.minCover([[1],[2]]).d === 2, 'minCover AND of two sole donors -> d=2');
ok(C.minCover([[1,2]]).d === 1, 'minCover OR (substitutable donors) -> d=1');
ok(C.minCover([[1,2],[2,3]]).d === 1, 'minCover stacking donor covers both -> d=1');
// Astra's counterexample (relay 248981991110131712): 15 alternatives for req 1, donor 14 serves both.
ok(C.minCover([[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],[14]]).d === 1, 'minCover no donor cap: donor 14 covers both -> d=1');
{ const big = []; for (let j = 0; j < 6; j++){ const r = []; for (let d = 0; d < 40; d++) if ((d + j) % 7 === 0) r.push(d); big.push(r); }
  // brute-force check against exhaustive search over donor subsets up to size 6
  const donors = [...new Set(big.flat())]; let best = 99;
  const rec = (i, chosen) => { if (chosen.length >= best) return; if (big.every(r => r.some(d => chosen.includes(d)))){ best = chosen.length; return; }
    if (i >= donors.length) return; rec(i+1, chosen.concat([donors[i]])); rec(i+1, chosen); };
  rec(0, []); ok(C.minCover(big).d === best, 'minCover DP equals exhaustive search on a 40-donor instance (d=' + best + ')'); }
ok(C.minCover([]).d === 0 && C.minCover([]).feasible, 'minCover: no requirement -> d=0 (autonomous)');
ok(!C.minCover([[1],[]]).feasible && isNaN(C.minCover([[1],[]]).d), 'minCover: uncoverable requirement -> infeasible, not d=0');
ok(C.classify(3, [[0,1],[1,2]]) === 'chain', 'classify path -> chain');
ok(C.classify(3, [[0,2],[1,2]]) === 'fan-in', 'classify collider -> fan-in');
ok(C.classify(3, [[0,1],[1,2],[2,0]]) === 'mesh', 'classify 3-cycle -> mesh');
ok(C.classify(4, [[0,1],[0,2],[0,3]]) === 'hub', 'classify star -> hub');
ok(C.findCycles(3, [[1],[0],[]], 8).length === 0, 'a 2-cycle is not a cycle >= 3');
const plan2 = {suppliers:[5,6], tests:[{name:'control'},{name:'remove',target:5},{name:'remove',target:6},{name:'none'}]};
const plan3 = {suppliers:[5,6,7], tests:[{name:'control'},{name:'remove',target:5},{name:'remove',target:6},{name:'remove',target:7},
  {name:'only',target:5},{name:'only',target:6},{name:'only',target:7},{name:'none'}]};
let r;
r = C.koInterpret(plan2, [10,9,9,0]); ok(r.status === 'one-suffices' && r.dLower === 1, 'KO k=2: either alone suffices -> one-suffices (redundant use)');
r = C.koInterpret(plan2, [10,1,1,0]); ok(r.status === 'at-least' && r.dLower === 2, 'KO k=2: both needed -> at least 2 (lower bound)');
r = C.koInterpret(plan2, [10,9,9,9]); ok(r.status === 'not-needed' && r.dLower === 0, 'KO: survives without tested suppliers -> not needed');
r = C.koInterpret(plan2, [0.2,9,9,0]); ok(r.status === 'inconclusive' && r.dLower === null, 'KO: extinct control -> inconclusive');
r = C.koInterpret(plan3, [10,9,9,9,1,1,1,0]); ok(r.status === 'at-least' && r.dLower === 2, 'KO k=3: any two of three -> at least 2, not exactly');
r = C.koInterpret(plan3, [10,1,9,9,8,1,1,0]); ok(r.status === 'one-suffices', 'KO k=3: supplier 5 alone suffices -> one-suffices');
r = C.koInterpret(plan3, [10,1,9,9,1,1,1,9]); ok(r.status === 'ambiguous', 'KO: single removal fatal but removing all harmless -> ambiguous');
ok(C.koPlan({nodes:[{g:1},{g:2},{g:3},{g:4},{g:5},{g:6},{g:7},{g:8}], used:[1,2,3,4,5,6,7].map(f => ({from:f, to:0}))}, 1).suppliers.length === 7, 'koPlan tests all 7 suppliers (no cap)');

// 2. determinism
const a1 = C.create({}, 7), a2 = C.create({}, 7);
for (let t = 0; t < 300; t++){ C.step(a1); C.step(a2); }
ok(a1.n === a2.n && a1.f[100] === a2.f[100], 'same seed -> identical state after 300 ticks');

// 3. runs under contrasting settings
const settings = {
  'default': {},
  'well-mixed': {mixing: 0.9},
  'no-leak, cheap genes': {leak: 0.05, geneCost: 0.1},
  'leaky, costly, particles': {leak: 0.8, geneCost: 0.9, particles: 6},
  'no incompatibility': {incompat: 0}
};
const hdr = ['setting','seed','cells','genotypes','auxGenoFrac','cov%','resolved%','plexi%','use>=2%','redund%','fanIn%','reqCyc','topoReq','topoUsed','algae','ms/tick'];
console.log(hdr.join('\t'));
let maxMs = 0;
for (const [name, p] of Object.entries(settings)){
  for (let sd = 1; sd <= nSeeds; sd++){
    const S = C.create(p, sd);
    const t0 = process.hrtime.bigint();
    for (let t = 0; t < T; t++) C.step(S);
    const ms = Number(process.hrtime.bigint() - t0)/1e6/T;
    maxMs = Math.max(maxMs, ms);
    const A = C.analyze(S);
    let algae = 0; for (let i = 0; i < S.n; i++) if (S.gen[i] & C.B_EUK) algae++;
    const auxG = A.nodes.filter(n => (n.g & 63) !== 63).length / Math.max(1, A.nodes.length);
    console.log([name, sd, S.n, A.nGenotypes, auxG.toFixed(2), (100*A.coverage).toFixed(0), (100*A.resolvedFrac).toFixed(0), (100*A.plexiFrac).toFixed(0), (100*A.multiUseFrac).toFixed(0),
      (100*A.redundantFrac).toFixed(0), (100*A.fanInFrac).toFixed(0), A.reqCycles.length, A.topoReq, A.topoUsed, algae, ms.toFixed(2)].join('\t'));
    if (name === 'default' && sd === 1){
      // sanity of invariants
      let bad = 0; for (let i = 0; i < S.f.length; i++) if (!(S.f[i] >= -1e-6) || !isFinite(S.f[i])) bad++;
      ok(bad === 0, 'fields finite and non-negative');
      let cnt = 0; for (let g = 0; g < C.NGT; g++) cnt += S.gcount[g];
      ok(cnt === S.n, 'genotype counts sum to population');
      // observer from 12 random samples
      const sm = []; const r = new C.RNG(99);
      for (let k = 0; k < 12; k++) sm.push(C.sample(S, 4 + r.next()*88, 4 + r.next()*56, 6, 'EV' + k));
      const O = C.observe(sm, {compl: 0.85});
      const cmp = C.compare(A, O);
      console.log('observer: nodes', O.nodes.length, 'edges', O.edges.length, 'topo', O.topo, '| eligible true edges', cmp.truthElig, 'of', cmp.truthAll,
        '(coverage', cmp.coverage.toFixed(2) + ') | conditional precision', cmp.precision.toFixed(2), 'recall', cmp.recall.toFixed(2),
        '| end-to-end recall', cmp.recallE2E.toFixed(2), '| unscorable predicted', cmp.unscorable, '| plexi agree', cmp.dAgree.toFixed(2), 'of', cmp.dN);
      ok(cmp.tp + cmp.fn === cmp.truthElig && cmp.truthElig <= cmp.truthAll, 'compare: precision/recall share one eligible pair set');
      ok(A.nodes.every(n => n.dUsed === 0 ? n.dEff === 0 : (n.dEff >= 1 - 1e-9 && n.dEff <= n.dUsed + 1e-9)), 'd_eff: 0 when unused, else 1 <= d_eff <= used');
      ok(A.nodes.every(n => n.dUsed !== 1 || Math.abs(n.dEff - 1) < 1e-9), 'd_eff = 1 when exactly one supplier is used');
      ok(A.nodes.every(n => n.reqs.every(r => r.ch !== 8)), 'H2 removal is never a requirement (fermentation is facultative in step())');
      const protoUsed = A.nodes.filter(n => (n.g & 63) === 63 && (n.g & (1<<6))).filter(n => n.dUsed > 0).length;
      ok(protoUsed === 0, 'phototrophic prototrophs register no required-channel suppliers (' + protoUsed + ')');
      const tsv = C.exportTSV(sm, 0.85, 17); fs.writeFileSync(process.env.PLEXI_EXPORT || '/tmp/plexi_export_test.tsv', tsv);
      console.log('export: wrote', tsv.split('\n').length, 'lines to', process.env.PLEXI_EXPORT || '/tmp/plexi_export_test.tsv');
      const Op = C.observe(sm, {compl: 0.85, mode: 'pcor'});
      console.log('observer (partial corr): edges', Op.edges.length, 'topo', Op.topo);
      // top recipients table
      A.nodes.slice(0, 8).forEach(n => console.log('  node', n.g.toString(2).padStart(11,'0'), 'n=' + n.count, 'used', n.dUsed, 'dEff', n.dEff.toFixed(2), 'dA', n.dA, 'dC', n.dC, 'dT', n.dT, n.inCycle ? 'cycle' : '', n.plexi ? 'PLEXI' : (n.redundantOnly ? 'redundant' : '')));
      // a knock-out on the recipient using the most suppliers
      const rec = A.nodes.filter(n => n.dUsed >= 1).sort((x, y) => y.dUsed - x.dUsed)[0];
      if (rec){
        const pl = C.koPlan(A, rec.g), H = 600, acc = pl.tests.map(() => 0);
        const cl = pl.tests.map(t => { const c = C.clone(S); C.knockOut(c, t.remove); return c; });
        for (let t = 0; t < H; t++) cl.forEach((c, i) => { C.step(c); if (t >= H*0.75) acc[i] += c.gcount[rec.g]; });
        const out = C.koInterpret(pl, acc.map(v => v/(H*0.25)));
        console.log('KO recipient', rec.g.toString(2), 'tested', out.usedN, 'min-cover dT', rec.dT, 'status', out.status, 'dLower', out.dLower,
          out.rows.map(r => r.name + ':' + (r.target!=null ? r.target : '') + '=' + r.ratio.toFixed(2)).join(' '));
      }
    }
  }
}
ok(maxMs < 20, 'mean step time < 20 ms (max over runs: ' + maxMs.toFixed(2) + ' ms)');
console.log(fails ? fails + ' FAILED' : 'all checks passed');
process.exit(fails ? 1 : 0);
