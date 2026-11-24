const fs=require('fs');const html=fs.readFileSync(require('path').join(__dirname,'index.html'),'utf8');
const C=new Function('module',html.match(/<script id="core">([\s\S]*?)<\/script>/)[1]+'\nreturn PlexiCore;')(undefined);
console.log(['leak','geneCost','mixing','incompat','particles','cells','algae','genotypes','cov%','resolved%','plexi%','use>=2%','fanIn%','reqCycles','usedCycles','topoReq'].join('\t'));
for (const leak of [0.3,0.8]) for (const geneCost of [0.3,0.9]) for (const mixing of [0,0.3]) for (const incompat of [0,1]) for (const particles of [0,6]){
  const S=C.create({leak,geneCost,mixing,incompat,particles},1); for(let t=0;t<3000;t++) C.step(S);
  const A=C.analyze(S); let alg=0; for(let i=0;i<S.n;i++) if(S.gen[i]&C.B_EUK) alg++;
  console.log([leak,geneCost,mixing,incompat,particles,S.n,alg,A.nGenotypes,(100*A.coverage).toFixed(0),(100*A.resolvedFrac).toFixed(0),(100*A.plexiFrac).toFixed(0),(100*A.multiUseFrac).toFixed(0),(100*A.fanInFrac).toFixed(0),A.reqCycles.length,A.usedCycles.length,A.topoReq].join('\t'));
}
