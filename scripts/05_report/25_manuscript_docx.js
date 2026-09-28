// Manuscript md → docx (Times New Roman 12pt, double spacing, continuous line numbers, A4, page numbers; references from docs/references.md, figures from letter_fig1/2.png)
// Run: after npm i docx (in a temporary folder),  NODE_PATH=<that folder>/node_modules node scripts/05_report/25_manuscript_docx.js "$(pwd)"
const fs = require("fs"), path = require("path");
const { Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType, Footer, PageNumber,
        LineNumberRestartFormat, PageBreak } = require("docx");
const ROOT = process.argv[2];
const MDNAME = process.argv[3] || "docs/manuscript_v3_GB.md";
const MD = fs.readFileSync(path.join(ROOT, MDNAME), "utf8");
const REF = fs.readFileSync(path.join(ROOT, "docs/references.md"), "utf8");
const OUT = path.join(ROOT, MDNAME.replace(/\.md$/, ".docx"));
const FONT = "Times New Roman", SIZE = 24, LINE = 480;   // 12 pt, double spacing

// inline markdown → runs (**bold**, *italic*, `code`)
function runs(text, base = {}) {
  const out = []; const re = /(\*\*[^*]+\*\*|\*[^*\s][^*]*\*|`[^`]+`)/g; let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index), ...base }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ text: t.slice(2, -2), bold: true, ...base }));
    else if (t.startsWith("`")) out.push(new TextRun({ text: t.slice(1, -1), font: "Courier New", ...base }));
    else out.push(new TextRun({ text: t.slice(1, -1), italics: true, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last), ...base }));
  return out;
}
const P = (text, o = {}) => new Paragraph({ children: runs(text, o.run || {}), spacing: { line: LINE, after: o.after ?? 0 },
  indent: o.indent, alignment: o.align, heading: o.heading, keepNext: o.keepNext });
const H = (text, level) => new Paragraph({ children: [new TextRun({ text, bold: true, size: level === 1 ? 28 : 24, italics: level === 2 })],
  heading: level === 1 ? HeadingLevel.HEADING_1 : HeadingLevel.HEADING_2, spacing: { line: LINE, before: level === 1 ? 240 : 120 }, keepNext: true });
const PB = () => new Paragraph({ children: [new PageBreak()] });

// ---- split markdown into sections
const lines = MD.split("\n");
let title = lines[0].replace(/^#\s+/, "");
const body = [], titlePage = []; let stop = false;
for (let i = 1; i < lines.length && !stop; i++) {
  const l = lines[i];
  if (/^## Notes for the authors/.test(l)) { stop = true; break; }
  if (/^> /.test(l)) { titlePage.push(l.replace(/^>\s+/, "")); continue; }          // title page: authors, affiliations, correspondence
  if (/^\*Short version/.test(l) || /^\*Version/.test(l) || /^---\s*$/.test(l) || !l.trim()) continue;   // working notes are not part of the document
  body.push(l);
}
// ^x^ → superscript run
function supRuns(text) { const out = []; const re = /\^([^^]+)\^/g; let last = 0, m;
  while ((m = re.exec(text))) { if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index) })); out.push(new TextRun({ text: m[1].replace(/\\\*/g, "*"), superScript: true })); last = m.index + m[0].length; }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last).replace(/\\\*/g, "*") })); return out; }
const children = [];
// title page
children.push(new Paragraph({ children: [new TextRun({ text: title, bold: true, size: 32 })], spacing: { line: LINE, after: 240 } }));
if (titlePage.length) titlePage.forEach(t => children.push(new Paragraph({ children: supRuns(t), spacing: { line: LINE } })));
else ["[Author names]", "[Affiliations]", "Correspondence: [corresponding author, e-mail]"].forEach(t =>
  children.push(new Paragraph({ children: [new TextRun({ text: t, color: "666666" })], spacing: { line: LINE } })));
children.push(PB());

let figuresInserted = false, sect = "";
const legends = [];
for (const l of body) {
  if (/^## /.test(l)) {
    sect = l.replace(/^##\s+/, "");
    if (sect === "Figure legends") {                         // references go before legends
      children.push(PB(), H("References", 1));
      REF.split("\n").filter(x => /^\d+\.\s/.test(x)).map(x => [parseInt(x), x.replace(/^\d+\.\s+/, "").replace(/\s*\*\([^]*?\)\*/g, "").trim()])
        .sort((a, b) => a[0] - b[0]).forEach(([n, t]) =>
          children.push(new Paragraph({ children: runs(`${n}. ${t}`), spacing: { line: LINE }, indent: { left: 360, hanging: 360 } })));
      children.push(PB());
    }
    if (sect === "Abstract" || sect === "Background") { if (sect === "Background") children.push(PB()); }
    children.push(H(sect, 1)); continue;
  }
  if (/^### /.test(l)) { children.push(H(l.replace(/^###\s+/, ""), 2)); continue; }
  if (sect === "Figure legends") { legends.push(l); children.push(P(l, { after: 240 })); continue; }
  children.push(P(l, { indent: (sect === "Abstract" || sect === "Declarations" || sect === "Methods" || sect === "Additional files") ? undefined : { firstLine: 360 } }));
}
// figures at the end, one per page with a short caption — skipped when NO_FIGS=1 (journal submission: figures are uploaded as separate files)
if (!process.env.NO_FIGS) {
const V4 = fs.existsSync(path.join(ROOT, "results/figures/draft_v4/panels/Fig3_composite.png")) && fs.existsSync(path.join(ROOT, "results/figures/draft_v4/panels/Fig4_composite.png"));
const figs = [["draft_v3/panels/Figure1_keynote.png", null, null, "Figure 1"], ["draft_v3/panels/Fig2_composite.png", null, null, "Figure 2"],
  [V4 ? "draft_v4/panels/Fig3_composite.png" : "draft_v3/panels/Fig3_composite.png", null, null, "Figure 3"], [V4 ? "draft_v4/panels/Fig4_composite.png" : "draft_v3/panels/Fig4_composite.png", null, null, "Figure 4"]];
for (const [f, w, h, cap] of figs) {
  const buf = fs.readFileSync(path.join(ROOT, "results/figures", f));
  let ww = w, hh = h; if (!ww) { ww = buf.readUInt32BE(16); hh = buf.readUInt32BE(20); }   // PNG header
  const W = 600, Hh = Math.round(W * hh / ww);
  children.push(PB(), new Paragraph({ children: [new TextRun({ text: cap, bold: true })], spacing: { after: 120 }, suppressLineNumbers: true }),
    new Paragraph({ suppressLineNumbers: true, children: [new ImageRun({ type: "png", data: buf, transformation: { width: W, height: Hh } })], alignment: AlignmentType.CENTER }));
}
}

const doc = new Document({
  creator: "LSG", title,
  styles: { default: { document: { run: { font: FONT, size: SIZE } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: FONT, size: 28, bold: true }, paragraph: { outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: FONT, size: 24, bold: true, italics: true }, paragraph: { outlineLevel: 1 } } ] },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } },
                  lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, suppressLineNumbers: true, children: [new TextRun({ children: [PageNumber.CURRENT] })] })] }) },
    children }]
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(OUT, b); console.log("WROTE", OUT, b.length); });
