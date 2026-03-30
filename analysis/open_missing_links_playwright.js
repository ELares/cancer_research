const { chromium } = require('playwright');

const urls = [
  'https://jamanetwork.com/journals/jama/articlepdf/2475463/jpc150007.pdf',
  'https://jamanetwork.com/journals/jama/articlepdf/2666504/jama_stupp_2017_oi_170143.pdf',
  'https://acsjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.3322/caac.21613',
  'http://www.cell.com/article/S0092867421002233/pdf',
  'http://www.clinical-breast-cancer.com/article/S1526820921000598/pdf',
  'https://doi.org/10.1016/j.biopha.2021.112512',
  'https://dr.ntu.edu.sg/bitstream/10356/171012/2/Polymeric%20STING%20Pro-agonists%20for%20Tumor-Specific%20Sonodynamic%20Immunotherapy.pdf',
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC12522170/pdf/nihms-2110220.pdf',
  'https://infoscience.epfl.ch/handle/20.500.14299/205010',
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC12704925/pdf/nihms-2122299.pdf',
  'https://doi.org/10.1016/s1470-2045(24)00508-4',
  'https://doi.org/10.1016/j.cell.2024.12.010',
  'https://doi.org/10.1016/j.celrep.2025.116738',
  'https://doi.org/10.1016/j.pdpdt.2026.105442',
];

(async () => {
  const context = await chromium.launchPersistentContext('', {
    headless: false,
  });

  const page = context.pages()[0] || await context.newPage();
  await page.goto(urls[0], { waitUntil: 'domcontentloaded' });

  for (const url of urls.slice(1)) {
    const tab = await context.newPage();
    await tab.goto(url, { waitUntil: 'domcontentloaded' });
  }

  console.log(`Opened ${urls.length} tabs. Browser will remain open until you close it.`);

  await new Promise(() => {});
})();
