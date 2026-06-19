#!/usr/bin/env node
// Build one self-contained HTML for the whole course, ready for Chrome -> PDF.
// Handles: LaTeX math (MathJax), Mermaid diagrams, tables, emoji, local figures.
const fs = require('fs'), path = require('path');
const MarkdownIt = require('markdown-it');

const ROOT = path.resolve(__dirname, '..');
const COURSE = path.join(ROOT, 'course');
const FIG = path.join(ROOT, 'figures');

const files = ['README.md','00_mental_map.md','01_biology_and_breeding.md','02_the_data.md',
  '03_phenotyping_spatial_BLUPs.md','04_genotyping_GBS_SNPs.md','05_quant_genetics_foundations.md',
  '06_genomic_relationship_matrix.md','07_GBLUP.md','08_RKHS_kernels.md','09_GWAS_FarmCPU.md',
  '10_GWAS_assisted_GP.md','11_NIRS_RSI.md','12_multitrait_models.md','13_cross_validation_accuracy.md',
  '14_across_cycles_updating.md','15_results_and_takehome.md','16_reproduce_it_yourself.md','GLOSSARY.md'];

const md = new MarkdownIt({ html:true, linkify:true, typographer:true });

function render(srcRaw){
  const store = [];
  const keep = s => '@@KEEP' + (store.push(s)-1) + 'KEEP@@';   // collision-proof token
  let src = srcRaw
    .replace(/```mermaid\n([\s\S]*?)```/g, (_,code)=> keep('<div class="mermaid">'+code.trim()+'</div>'))
    .replace(/\$\$([\s\S]*?)\$\$/g, (_,m)=> keep('<div class="math-display">$$'+m+'$$</div>'))
    .replace(/\$([^\$\n]+?)\$/g, (_,m)=> keep('<span class="math-inline">$'+m+'$</span>'));
  let html = md.render(src);
  html = html.replace(/@@KEEP(\d+)KEEP@@/g, (_,i)=> store[+i]);          // restore
  html = html.replace(/src="(?:\.\.\/)?figures\/([^"]+)"/g, (_,f)=> 'src="file://'+FIG+'/'+f+'"'); // figures
  return html;
}

let sections = '', toc = '';
files.forEach((f,i)=>{
  const raw = fs.readFileSync(path.join(COURSE, f), 'utf8');
  const m = raw.match(/^#\s+(.+)$/m);
  const title = m ? m[1].replace(/[`*]/g,'') : f;
  const id = 'sec'+i;
  toc += '<li><a href="#'+id+'">'+md.renderInline(title)+'</a></li>';
  sections += '<section class="lesson" id="'+id+'">'+render(raw)+'</section>';
});

const head = `<!doctype html><html><head><meta charset="utf-8">
<title>Genomic Prediction in Black Beans — Course</title>
<script>
window.MathJax = { tex:{ inlineMath:[['$','$']], displayMath:[['$$','$$']] },
  options:{ skipHtmlTags:['script','noscript','style','textarea','pre','code'] }, svg:{ fontCache:'global' } };
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js" async></script>
<script src="https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.min.js"></script>
<script>mermaid.initialize({ startOnLoad:true, theme:'neutral', flowchart:{useMaxWidth:true} });</script>
<style>
  @page { size:A4; margin:16mm 15mm; }
  html { font-size:11pt; }
  body { font-family:"Georgia","Times New Roman",serif; color:#1a1a1a; line-height:1.5; }
  .cover { text-align:center; page-break-after:always; padding-top:28vh; }
  .cover h1 { font-size:30pt; color:#14532d; border:none; margin-bottom:6pt; }
  .cover .sub { font-size:13pt; color:#444; }
  .cover .meta { margin-top:40px; font-size:10pt; color:#666; }
  .toc { page-break-after:always; } .toc h2 { color:#14532d; } .toc ol { line-height:1.9; }
  .lesson { page-break-before:always; }
  h1 { font-size:20pt; color:#14532d; border-bottom:3px solid #14532d; padding-bottom:4px; }
  h2 { font-size:15pt; color:#15803d; margin-top:1.2em; border-bottom:1px solid #ccc; }
  h3 { font-size:12.5pt; color:#166534; }
  a { color:#1d4ed8; text-decoration:none; }
  code { font-family:"Menlo","Consolas",monospace; font-size:9.5pt; background:#f3f4f6; padding:1px 4px; border-radius:3px; }
  pre { background:#f6f8fa; border:1px solid #e5e7eb; border-radius:6px; padding:10px 12px;
        white-space:pre-wrap; word-wrap:break-word; font-size:8.5pt; line-height:1.35; page-break-inside:avoid; }
  pre code { background:none; padding:0; font-size:8.5pt; }
  blockquote { border-left:4px solid #86efac; background:#f0fdf4; margin:0.6em 0; padding:6px 14px; }
  table { border-collapse:collapse; width:100%; margin:0.8em 0; font-size:9.5pt; page-break-inside:avoid; }
  th,td { border:1px solid #cbd5e1; padding:5px 8px; text-align:left; vertical-align:top; }
  th { background:#dcfce7; } tr:nth-child(even){ background:#f8fafc; }
  img { max-width:96%; height:auto; display:block; margin:10px auto; page-break-inside:avoid; }
  .mermaid { text-align:center; margin:14px 0; page-break-inside:avoid; }
  .math-display { overflow-x:auto; } hr { border:none; border-top:1px solid #e5e7eb; margin:1.4em 0; }
</style></head><body>
<div class="cover">
  <h1>Genomic Prediction in Black Beans</h1>
  <div class="sub">A beginner-friendly course that reproduces a real GWAS &amp; genomic-selection study</div>
  <div class="sub" style="font-size:11pt;margin-top:14px;">Built around Izquierdo, Wright &amp; Cichy (2025), <i>G3</i> jkaf007</div>
  <div class="meta">19 documents · toy-first worked examples · reproduced figures · glossary</div>
</div>
<nav class="toc"><h2>Contents</h2><ol>`;

const tail = '</ol></nav>' + sections + '</body></html>';
fs.writeFileSync(path.join(ROOT,'course.html'), head + toc + tail);
console.log('Wrote course.html ('+((head+toc+tail).length/1024).toFixed(0)+' KB, '+files.length+' docs)');
