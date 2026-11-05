const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const out = process.argv[2];
  for (const [cfg, name, h] of [["cfg_overview.json","overview",900],["cfg_gene.json","gene",1100],["cfg_amplicon.json","amplicon",1100]]) {
    const p = await b.newPage({ viewport: { width: 1400, height: h }, deviceScaleFactor: 2 });
    p.on('console', m => { if (m.type()==='error') console.log(name, m.text().slice(0,200)); });
    await p.goto(`http://127.0.0.1:8765/index.html?cfg=${cfg}`);
    await p.waitForFunction(() => window.igvReady === true, null, { timeout: 60000 });
    await p.waitForTimeout(7000);
    await p.screenshot({ path: `${out}/igv_${name}.png`, fullPage: true });
    await p.close();
  }
  await b.close();
})();
