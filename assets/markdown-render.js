/* =======================================================================
   UMtts markdown renderer
   Goal: render the source .md VERBATIM (no truncation, no paraphrase,
   no collapsing, no "read more"). Every line ends up in the DOM.
   Also builds a table of contents and tracks the current section.
   ======================================================================= */

(function () {
  'use strict';

  function escHTML(s) {
    return s
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
  }

  // Inline transforms - italic, bold, inline code - but do not drop any chars.
  function inline(s) {
    // `code`
    s = s.replace(/`([^`]+)`/g, (_, t) => '<code>' + t + '</code>');
    // **bold**
    s = s.replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>');
    // *italic*  (but not inside already-escaped markers)
    s = s.replace(/(^|[\s(])\*([^*\n]+)\*(?=[\s.,;:)]|$)/g, '$1<em>$2</em>');
    return s;
  }

  // Detect heading level from markdown line
  function headingLevel(line) {
    const m = line.match(/^(#{1,6})\s+(.*)$/);
    return m ? { level: m[1].length, text: m[2] } : null;
  }

  // Render a markdown string into HTML that preserves every source line.
  function render(md) {
    // normalize line endings - CRLF → LF
    md = md.replace(/\r\n?/g, '\n');
    const lines = md.split('\n');
    let out = [];
    let i = 0;
    let tocEntries = [];
    let slugCount = {};

    function slug(text) {
      let s = text.toLowerCase()
        .replace(/[^\w\s-]/g, '')
        .replace(/\s+/g, '-')
        .replace(/-+/g, '-')
        .replace(/^-|-$/g, '') || 'section';
      if (slugCount[s]) { slugCount[s]++; s = s + '-' + slugCount[s]; }
      else { slugCount[s] = 1; }
      return s;
    }

    // Look for table: a sequence of lines with | and the second being the separator
    function isTableRow(line) { return /\|/.test(line); }
    function isTableSep(line) { return /^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)+\|?\s*$/.test(line); }

    // Split a row by |
    function splitRow(line) {
      // strip leading/trailing |
      let s = line.trim();
      if (s.startsWith('|')) s = s.slice(1);
      if (s.endsWith('|')) s = s.slice(0, -1);
      return s.split('|').map(c => c.trim());
    }

    // equation detector - lines that look like math (contain specific symbols)
    // We treat them as a display equation if:
    //   - they aren't a heading or list item
    //   - they contain ψ, ∇, ∂, ∫, ∑, =, →, ², ³, ⁿ, φ, √, etc., AND no natural-language run
    // To be safe we use a light heuristic and fall back to paragraphs.
    const MATH_CHARS = /[ψφΦκΚβ∇∂∫∑πωαρℏ·×→⇒≥≤≠≈∞]/;
    const MATH_HEAVY = /[ψφΦ∇∂∫∑→⇒∞αωπ]|[⁰¹²³⁴⁵⁶⁷⁸⁹]/;
    function looksLikeEquation(line) {
      const t = line.trim();
      if (!t) return false;
      if (t.length > 200) return false;
      if (/^(?:[*\-•]\s+|\d+\.\s+)/.test(t)) return false; // actual list item only
      if (/^#/.test(t)) return false;
      const hasMath = MATH_HEAVY.test(t);
      const hasEq = /=|⇒|→/.test(t);
      return hasMath && hasEq && t.length < 180;
    }

    while (i < lines.length) {
      const line = lines[i];
      const raw = line;
      const trimmed = line.trim();

      // skip no-op empty lines (they already separate blocks)
      if (trimmed === '') { i++; continue; }

      // headings - explicit markdown (#, ##, ###)
      const h = headingLevel(line);
      if (h) {
        const id = slug(h.text);
        tocEntries.push({ level: h.level, text: h.text, id });
        out.push(`<h${h.level} id="${id}">${inline(escHTML(h.text))}</h${h.level}>`);
        i++;
        continue;
      }

      // Part/APPENDIX heading (ALL CAPS "PART n:" or "APPENDIX X:") - treat as h2
      const partMatch = trimmed.match(/^(PART\s+\d+[A-Z]?:\s+.*|APPENDIX\s+[A-Z]:\s+.*|DECLARATION OF ZERO FREE PARAMETERS|THE ZERO-INPUT FRAMEWORK|THE TEN COMMANDMENTS OF MASS HARMONICS|ABSTRACT|HISTORICAL NOTATIONS AND SYMBOL EQUIVALENCY)$/);
      if (partMatch) {
        const id = slug(trimmed);
        tocEntries.push({ level: 2, text: trimmed, id });
        out.push(`<h2 id="${id}">${inline(escHTML(trimmed))}</h2>`);
        i++;
        continue;
      }

      // Subsection like "3.1 Expansion of the Lagrangian" or "4.2 DERIVATION: φ as..."
      const subMatch = trimmed.match(/^(\d+[A-Z]?\.\d+(?:\.\d+)?(?::|\s+).+|Axiom\s+\d+:.*|Step\s+\d+:.*|Reason\s+\d+\s*[---:].*|Commandment\s+[IVX]+.*)$/);
      if (subMatch && !looksLikeEquation(trimmed)) {
        // 3-part numbers (like "3.1.2") count as level 4; 2-part as level 3
        const dotCount = (trimmed.match(/\./g) || []).length;
        const level = dotCount >= 2 ? 4 : 3;
        const id = slug(trimmed);
        tocEntries.push({ level, text: trimmed, id });
        out.push(`<h${level} id="${id}">${inline(escHTML(trimmed))}</h${level}>`);
        i++;
        continue;
      }

      // STOP / CRITICAL callout
      if (trimmed === 'STOP' || trimmed === 'CRITICAL' || trimmed === 'STOP. FULL READ REQUIRED.') {
        // collect consecutive callout lines (until a blank or a heading)
        let buf = [trimmed];
        let j = i + 1;
        while (j < lines.length) {
          const nxt = lines[j].trim();
          if (nxt === '') break;
          if (headingLevel(lines[j])) break;
          if (/^(PART|APPENDIX)\s/.test(nxt)) break;
          buf.push(nxt);
          j++;
        }
        out.push('<div class="callout callout--stop"><h4>⚠ Critical</h4>' +
                 buf.map(b => '<p>' + inline(escHTML(b)) + '</p>').join('') +
                 '</div>');
        i = j;
        continue;
      }

      // Table - current line has |, next line is separator
      if (isTableRow(line) && i + 1 < lines.length && isTableSep(lines[i + 1])) {
        const header = splitRow(line);
        i += 2;
        let rows = [];
        while (i < lines.length && isTableRow(lines[i])) {
          rows.push(splitRow(lines[i]));
          i++;
        }
        let html = '<table><thead><tr>' +
          header.map(c => '<th>' + inline(escHTML(c)) + '</th>').join('') +
          '</tr></thead><tbody>' +
          rows.map(r => '<tr>' + r.map(c => '<td>' + inline(escHTML(c)) + '</td>').join('') + '</tr>').join('') +
          '</tbody></table>';
        out.push(html);
        continue;
      }

      // tab-separated table - like "Old\tNew\tWhat it names" then rows
      if (/\t/.test(line) && i + 1 < lines.length && /\t/.test(lines[i + 1])) {
        // Collect contiguous tab-bearing lines
        let buf = [line];
        let j = i + 1;
        while (j < lines.length && /\t/.test(lines[j]) && lines[j].trim() !== '') {
          buf.push(lines[j]);
          j++;
        }
        if (buf.length >= 2) {
          const header = buf[0].split('\t').map(c => c.trim());
          const rows = buf.slice(1).map(r => r.split('\t').map(c => c.trim()));
          // pad short rows
          rows.forEach(r => { while (r.length < header.length) r.push(''); });
          out.push('<table><thead><tr>' +
            header.map(c => '<th>' + inline(escHTML(c)) + '</th>').join('') +
            '</tr></thead><tbody>' +
            rows.map(r => '<tr>' + r.map(c => '<td>' + inline(escHTML(c)) + '</td>').join('') + '</tr>').join('') +
            '</tbody></table>');
          i = j;
          continue;
        }
      }

      // bullet list
      if (/^\s*[•\-*]\s+/.test(line)) {
        let items = [];
        while (i < lines.length && /^\s*[•\-*]\s+/.test(lines[i])) {
          items.push(lines[i].replace(/^\s*[•\-*]\s+/, ''));
          i++;
        }
        out.push('<ul>' + items.map(it => '<li>' + inline(escHTML(it)) + '</li>').join('') + '</ul>');
        continue;
      }

      // numbered list (1. 2. ...)
      if (/^\s*\d+\.\s+/.test(line) && !/^\d+[A-Z]?\.\d/.test(trimmed)) {
        let items = [];
        while (i < lines.length && /^\s*\d+\.\s+/.test(lines[i])) {
          items.push(lines[i].replace(/^\s*\d+\.\s+/, ''));
          i++;
        }
        out.push('<ol>' + items.map(it => '<li>' + inline(escHTML(it)) + '</li>').join('') + '</ol>');
        continue;
      }

      // equation block (single-line, math-heavy)
      if (looksLikeEquation(trimmed)) {
        out.push('<div class="equation">' + inline(escHTML(trimmed)) + '</div>');
        i++;
        continue;
      }
 
      // Raw HTML passthrough
      if (trimmed === '<!-- RAWHTML -->') {
        i++;
        let buf = [];
        while (i < lines.length && lines[i].trim() !== '<!-- /RAWHTML -->') {
          buf.push(lines[i]);
          i++;
        }
        if (i < lines.length) i++; // consume closing marker
        out.push(buf.join('\n'));
        continue;
      }

      // quote (start with > or a leading quote char)
      if (/^>\s*/.test(trimmed) || /^"[^"]*"$/.test(trimmed)) {
        let buf = [trimmed.replace(/^>\s*/, '')];
        let j = i + 1;
        while (j < lines.length && /^>\s*/.test(lines[j])) {
          buf.push(lines[j].replace(/^>\s*/, ''));
          j++;
        }
        out.push('<blockquote>' + buf.map(b => '<p>' + inline(escHTML(b)) + '</p>').join('') + '</blockquote>');
        i = j;
        continue;
      }

      // Plain paragraph - collect consecutive non-blank, non-structured lines
      let pbuf = [line];
      i++;
      while (i < lines.length) {
        const nxt = lines[i];
        if (nxt.trim() === '') break;
        if (headingLevel(nxt)) break;
        if (/^\s*[•\-*]\s+/.test(nxt)) break;
        if (/^\s*\d+\.\s+/.test(nxt) && !/^\d+[A-Z]?\.\d/.test(nxt.trim())) break;
        if (isTableRow(nxt) && isTableSep(lines[i + 1] || '')) break;
        if (/\t/.test(nxt) && i + 1 < lines.length && /\t/.test(lines[i + 1])) break;
        if (looksLikeEquation(nxt.trim())) break;
        if (/^(PART|APPENDIX)\s/.test(nxt.trim())) break;
        const subM = nxt.trim().match(/^(\d+[A-Z]?\.\d+(?:\.\d+)?(?::|\s+).+|Axiom\s+\d+:.*|Step\s+\d+:.*|Reason\s+\d+\s*[---:].*)$/);
        if (subM && !looksLikeEquation(nxt.trim())) break;
        pbuf.push(nxt);
        i++;
      }
      // Join with spaces but preserve explicit line breaks with <br>
      const para = pbuf.map(l => inline(escHTML(l))).join('<br>');
      out.push('<p>' + para + '</p>');
    }

    return { html: out.join('\n'), toc: tocEntries };
  }

  // Build a TOC element
  function buildTOC(entries, targetEl) {
    const ol = document.createElement('ol');
    entries.forEach(e => {
      if (e.level > 3) return;
      const li = document.createElement('li');
      li.className = 'lvl-' + e.level;
      const a = document.createElement('a');
      a.href = '#' + e.id;
      a.textContent = e.text.length > 72 ? e.text.slice(0, 70) + '…' : e.text;
      a.dataset.target = e.id;
      li.appendChild(a);
      ol.appendChild(li);
    });
    targetEl.appendChild(ol);
  }

  // Scrollspy - mark current section without moving the page.
  // The old build used page-level TOC auto-scrolling. On phones, that can scroll
  // the whole document back toward the top because the TOC lives above the body.
  // This version only scrolls the TOC container on desktop, and never the window.
  function scrollspy(tocEl, bodyEl) {
    const headings = Array.from(bodyEl.querySelectorAll('h1, h2, h3'));
    const links = Array.from(tocEl.querySelectorAll('a'));
    if (!headings.length) return;

    const linkById = {};
    links.forEach(a => { linkById[a.dataset.target] = a; });

    let ticking = false;

    function update() {
      ticking = false;
      const scrollY = window.scrollY + 140;
      let current = headings[0];

      for (const h of headings) {
        if (h.offsetTop <= scrollY) current = h;
        else break;
      }

      links.forEach(a => a.classList.remove('is-current'));

      if (current && linkById[current.id]) {
        const link = linkById[current.id];
        link.classList.add('is-current');

        const isMobile = window.matchMedia('(max-width: 980px)').matches;
        if (!isMobile && tocEl.scrollHeight > tocEl.clientHeight) {
          const linkTop = link.offsetTop;
          const linkBottom = linkTop + link.offsetHeight;
          const viewTop = tocEl.scrollTop;
          const viewBottom = viewTop + tocEl.clientHeight;

          if (linkTop < viewTop || linkBottom > viewBottom) {
            tocEl.scrollTop = Math.max(0, linkTop - Math.floor(tocEl.clientHeight * 0.45));
          }
        }
      }
    }

    function requestUpdate() {
      if (!ticking) {
        ticking = true;
        window.requestAnimationFrame(update);
      }
    }

    window.addEventListener('scroll', requestUpdate, { passive: true });
    window.addEventListener('resize', requestUpdate, { passive: true });
    update();
  }

  // Public API
  window.UMtts = {
    renderInto: function (opts) {
      const source = document.getElementById(opts.sourceId).textContent;
      // Preserve the source byte-count for verification
      const expectedLen = source.length;
      const { html, toc } = render(source);
      const bodyEl = document.getElementById(opts.bodyId);
      bodyEl.insertAdjacentHTML('beforeend', html);
      if (opts.tocId) {
        const tocEl = document.getElementById(opts.tocId);
        buildTOC(toc, tocEl);
        scrollspy(tocEl, bodyEl);
      }

      
    // Integrity tag - count characters with rigor: content vs formatting
const bodyText = bodyEl.innerText.replace(/\s+/g, '');

// Strip RAWHTML comment blocks from source before counting
const sourceWithoutRawhtml = source.replace(/<!--\s*RAWHTML\s*-->[\s\S]*?<!--\s*\/RAWHTML\s*-->/g, '');
const sourceCondensed = sourceWithoutRawhtml.replace(/\s+/g, '');

// Calculate formatting overhead (markdown syntax characters)
const contentChars = bodyText.length;
const totalChars = sourceCondensed.length;
const formattingChars = totalChars - contentChars;

console.log('[UMtts] source bytes:', expectedLen,
            'rendered content:', contentChars,
            'total source:', totalChars,
            'formatting symbols:', formattingChars);

if (opts.integrityEl) {
  // Format: "580 / 595 • Formatting: 15"
  // This keeps "Source characters rendered:" as the label context
  const integrityText = contentChars.toLocaleString() + ' / ' + totalChars.toLocaleString() + 
    ' • Formatting: ' + formattingChars.toLocaleString();
  
  document.getElementById(opts.integrityEl).textContent = integrityText;
}
}
      
};
  
})();
