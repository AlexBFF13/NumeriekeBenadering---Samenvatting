/**
 * Converts Claude-style HTML snippets (no DOCTYPE, uses CSS vars) to
 * standalone HTML files that work in any browser.
 *
 * Usage:
 *   node tools/convert.js              — converts all HTML in Algoritmes/ once
 *   node tools/convert.js --watch      — watches for changes and auto-converts
 */

const fs   = require('fs');
const path = require('path');

const WATCH_DIR = path.join(__dirname, '..', 'Algoritmes');

// Map every Claude CSS variable to a real value
const CSS_VARS = {
  '--font-sans':                  '-apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif',
  '--font-mono':                  '"Menlo", "Monaco", "Courier New", monospace',
  '--color-background-primary':   '#ffffff',
  '--color-background-secondary': '#f8f9fb',
  '--color-background-success':   '#eefaf3',
  '--color-background-danger':    '#fef2f2',
  '--color-background-warning':   '#fef5e7',
  '--color-text-primary':         '#111111',
  '--color-text-secondary':       '#6b7280',
  '--color-text-tertiary':        '#adb5c0',
  '--color-text-success':         '#166534',
  '--color-text-danger':          '#991b1b',
  '--color-border-primary':       '#c0c4cc',
  '--color-border-secondary':     '#c0c4cc',
  '--color-border-tertiary':      '#e2e5ec',
  '--border-radius-sm':           '4px',
  '--border-radius-md':           '6px',
  '--border-radius-lg':           '10px',
};

function isClaudeSnippet(content) {
  // Claude snippets start with <style> or whitespace + <style>, no DOCTYPE
  return !content.trimStart().startsWith('<!') &&
         (content.includes('var(--') || content.trimStart().startsWith('<style'));
}

function replaceCSSVars(content) {
  let out = content;
  for (const [varName, value] of Object.entries(CSS_VARS)) {
    // Replace var(--name) globally
    out = out.replaceAll(`var(${varName})`, value);
  }
  return out;
}

function deriveTitle(filename) {
  return path.basename(filename, '.html')
    .replace(/_/g, ' ')
    .replace(/-/g, ' ')
    .replace(/\b\w/g, c => c.toUpperCase());
}

function wrap(content, title) {
  return `<!DOCTYPE html>
<html lang="nl">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>${title}</title>
  <style>
    *, *::before, *::after { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 24px 20px 48px;
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      background: #f5f5f5;
      color: #111;
      min-height: 100vh;
    }
    .sr-only { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0,0,0,0); }
  </style>
</head>
<body>
${content}
</body>
</html>`;
}

function convertFile(filePath) {
  const raw = fs.readFileSync(filePath, 'utf8');

  if (!isClaudeSnippet(raw)) return false; // already a full HTML file, skip

  const noVars  = replaceCSSVars(raw);
  const title   = deriveTitle(filePath);
  const full    = wrap(noVars, title);

  fs.writeFileSync(filePath, full, 'utf8');
  console.log(`[converted] ${path.basename(filePath)}`);
  return true;
}

function convertAll() {
  const files = fs.readdirSync(WATCH_DIR).filter(f => f.endsWith('.html'));
  let count = 0;
  for (const f of files) {
    if (convertFile(path.join(WATCH_DIR, f))) count++;
  }
  console.log(`Done. ${count} file(s) converted.`);
}

// ── Entry point ──────────────────────────────────────────────────────────────

if (process.argv.includes('--watch')) {
  convertAll();
  console.log(`\nWatching ${WATCH_DIR} for changes…`);

  fs.watch(WATCH_DIR, (event, filename) => {
    if (!filename || !filename.endsWith('.html')) return;
    const filePath = path.join(WATCH_DIR, filename);
    if (!fs.existsSync(filePath)) return;

    // Small delay so VS Code finishes writing before we read
    setTimeout(() => {
      try { convertFile(filePath); }
      catch (e) { console.error(`Error processing ${filename}:`, e.message); }
    }, 150);
  });
} else {
  convertAll();
}
