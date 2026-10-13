// dsXNAbuilder front end. All computation runs in worker.js (Pyodide); this file is UI only.
// Issue messages come in English from the Python package, each with a stable code.

const $ = (s) => document.querySelector(s);
const $$ = (s) => [...document.querySelectorAll(s)];
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

const EXAMPLES = [
  { id: 'dickerson', label: 'Dickerson dodecamer (B-DNA)', top: 'CGCGAATTCGCG', bottom: 'GCGCTTAAGCGC', model: 4 },
  { id: 'arna', label: 'RNA duplex (A-RNA)', top: 'GGCGCUUGCGUC', bottom: 'CCGCGAACGCAG', model: 20 },
  { id: 'hybrid', label: 'RNA:DNA hybrid', top: 'r(GACUAGCUGA)', bottom: 'd(CTGATCGACT)', model: 22 },
  { id: 'adna', label: 'A-DNA', top: 'GGTATACC', bottom: 'CCATATGG', model: 54 },
  { id: 'zdna', label: 'Z-DNA (CG)₃', top: 'CGCGCG', bottom: 'GCGCGC', model: 15 },
  { id: 'mismatch', label: 'Mismatch + overhang', top: 'CGCGAGTTCGCGAA', bottom: 'GCGCTAAAGCGC--', model: 55 },
  { id: 'abasic', label: 'Abasic site', top: 'CGCGA-TTCGCG', bottom: 'GCGCTTAAGCGC', model: 4 },
  { id: 'sirna', label: 'siRNA 2′-OMe/2′-F + PS', top: 'mG*mA*fC*mU*fA*mG*fC*mU*mG*mA', bottom: 'C*U*GAUCGAC*U', model: 20 },
  { id: 'lna', label: 'LNA gapmer : RNA', top: 'r(GCAAUGUACUGAAUGC)', bottom: '+C*+G*+T*dT*dA*dC*dA*dT*dG*dA*dC*dT*dT*+A*+C*+G', model: 22 },
  { id: 'moe', label: 'MOE gapmer : RNA', top: 'r(GACUAGCUGA)', bottom: 'eC*eT*eG*dA*dT*dC*dG*eA*eC*eT', model: 22 },
  { id: 'mc', label: '5mC + PS on B-DNA', top: 'CGC^5mCGAA*2TTCGCG', bottom: 'GCGCTT*1AAGCGC', model: 4 },
  // explicit modXNA triplet {BACKBONE/SUGAR/BASE:pairing letter}: phosphorothioate,
  // 2'-amino-2'-deoxyribose and 5,6-dihydrocytosine, with an abasic site and overhangs
  { id: 'modxna', label: 'modXNA code', top: 'CG{PS2/2AD/RYC:C}GAA*2TTC-CG--', bottom: 'GCGCTT*1AAGCGCAA', model: 4 },
];

const STATUS = {
  loading_python: 'Loading Python runtime…', loading_numpy: 'Loading numpy…', loading_app: 'Loading dsXNAbuilder…',
};
const POLICY = {
  any: () => 'generic: any sequence',
  swap: (m) => `derived for ${m.native[0]}; other bases are swapped onto its backbone`,
  pattern: (m) => `repeat ${m.native[0]}; keep its purine/pyrimidine pattern`,
};
// marks the generic fiber models (built for any sequence) in the menu and the table
const GENERIC = '★';

const state = {
  ready: false, models: [], byId: new Map(), last: null, stage: null, comp: null,
  repr: 'sticks', color: 'base', spinning: false, validateTimer: null, validSeq: 0,
  residues: new Map(), names: {}, buildSeq: 0,
};

// ------------------------------------------------------------------ worker
const worker = new Worker(new URL('worker.js?v=0.2.0-dev', import.meta.url), { type: 'module' });
let seq = 0;
const pending = new Map();
worker.onmessage = (ev) => {
  const d = ev.data;
  if (d.status) return onStatus(d);
  const p = pending.get(d.id);
  if (!p) return;
  pending.delete(d.id);
  d.ok ? p.resolve(d.result) : p.reject(new Error(d.error));
};
const call = (name, payload = {}) => new Promise((resolve, reject) => {
  const id = ++seq;
  pending.set(id, { resolve, reject });
  worker.postMessage({ id, name, payload });
});

function onStatus(d) {
  if (d.status === 'ready') {
    state.ready = true;
    $('#ver').textContent = 'v' + d.version;
    $('#status').textContent = '';
    call('models').then((info) => { setModels(info); applyURL(); });
    if (!$('#view-library').classList.contains('hidden')) loadLibrary();
  } else if (d.status === 'error') {
    $('#status').innerHTML = `<span class="error">Could not start the engine: ${esc(d.error)}</span>`;
  } else {
    $('#status').innerHTML = `<span class="spinner"></span>${esc(STATUS[d.status] || d.status)}`;
  }
}

// ------------------------------------------------------------------ models
function setModels(info) {
  state.models = info.models;
  state.excluded = info.excluded;
  state.byId = new Map(info.models.map((m) => [m.id, m]));
  fillModelSelect();
  renderModelTables();
  showModelInfo();
  validateSoon(0);
}
function modelLabel(m) {
  const pol = m.policy === 'any' ? 'generic' : m.policy === 'swap' ? 'any (base swap)' : 'pattern';
  return `${m.policy === 'any' ? GENERIC + ' ' : ''}#${m.id} · ${m.family} · ${m.twist.toFixed(1)}°/${m.rise.toFixed(2)} Å · ${m.native[0]} · ${pol}`;
}
function fillModelSelect() {
  const sel = $('#model');
  const cur = +sel.value || 4;
  const fam = new Map();
  for (const m of state.models) {
    if (!fam.has(m.family)) fam.set(m.family, []);
    fam.get(m.family).push(m);
  }
  const order = [...fam.keys()].sort((a, b) => famRank(a) - famRank(b) || a.localeCompare(b));
  sel.innerHTML = order.map((f) => `<optgroup label="${esc(f)}">${fam.get(f)
    .sort((a, b) => (a.policy === 'any' ? 0 : 1) - (b.policy === 'any' ? 0 : 1) || a.id - b.id)
    .map((m) => `<option value="${m.id}"${m.policy === 'any' ? ' class="generic"' : ''}>${esc(modelLabel(m))}</option>`).join('')}</optgroup>`).join('');
  sel.value = state.byId.has(cur) ? cur : 4;
}
function famRank(f) {
  const r = ['B-DNA', 'A-DNA', 'A-RNA', 'DNA:RNA hybrid', "A'-RNA", 'C-DNA', 'Z-DNA', 'D-DNA'];
  const i = r.indexOf(f);
  return i < 0 ? 99 : i;
}
function showModelInfo() {
  const m = state.byId.get(+$('#model').value);
  if (!m) return;
  $('#model').classList.toggle('generic', m.policy === 'any');
  const info = `${m.twist.toFixed(1)}° and ${m.rise.toFixed(2)} Å per bp, ${m.handedness}-handed · ` +
    `strands ${m.strands[0]}/${m.strands[1]} · sugars ${m.pucker[0] || '?'}/${m.pucker[1] || '?'}`;
  const badge = m.policy === 'any' ? `<span class="badge generic">${GENERIC} generic</span> ` : '';
  $('#modelinfo').innerHTML = badge + esc(`${info} · ${POLICY[m.policy](m)} · ${m.description}`);
}
function renderModelTables() {
  const head = '<thead><tr><th>#</th><th>family</th><th>twist (°/bp)</th><th>rise (Å/bp)</th><th>hand</th><th>repeat</th><th>strands</th><th>sequence</th><th>3DNA description</th><th></th></tr></thead>';
  const rows = state.models.map((m) => {
    const cls = [m.id === +$('#model').value ? 'sel' : '', m.policy === 'any' ? 'generic' : ''].filter(Boolean).join(' ');
    const badge = m.policy === 'any' ? ` <span class="badge generic" title="generic model: built for any sequence">${GENERIC} generic</span>` : '';
    return `<tr${cls ? ` class="${cls}"` : ''}>
    <td class="num">${m.id}</td><td>${esc(m.family)}${badge}</td><td class="num">${m.twist.toFixed(2)}</td><td class="num">${m.rise.toFixed(3)}</td>
    <td>${m.handedness}</td><td class="num">${m.k}</td><td>${m.strands.join('/')}</td>
    <td class="seq">${esc(m.native[0])}/${esc(m.native[1])}</td><td class="desc">${esc(m.description)}</td>
    <td><button class="btn ghost small" data-use="${m.id}">use</button></td></tr>`;
  }).join('');
  $('#modeltable').innerHTML = head + `<tbody>${rows}</tbody>`;
  $('#excludedtable').innerHTML = '<thead><tr><th>#</th><th>family</th><th>3DNA description</th><th>reason</th></tr></thead><tbody>' +
    (state.excluded || []).map((m) => `<tr><td class="num">${m.id}</td><td>${esc(m.family)}</td><td class="desc">${esc(m.description)}</td><td class="desc">${esc((m.issues || []).join('; '))}</td></tr>`).join('') + '</tbody>';
  for (const b of $$('[data-use]')) b.onclick = () => { $('#model').value = b.dataset.use; showModelInfo(); selectTab('build'); validateSoon(0); };
}

// ------------------------------------------------------------------ modXNA library
// Every fragment of the modXNA release, by section; "prebuilt" ones are in the residue
// library (every combination of them is available) and come with their notation.
async function loadLibrary() {
  if (state.fragments || !state.ready) return;
  try {
    state.fragments = await call('fragments');
  } catch (e) {
    $('#libtables').innerHTML = `<p class="error">${esc(e.message)}</p>`;
    return;
  }
  const t = state.fragments;
  $('#libsource').innerHTML = `${t.n_fragments} fragments of modXNA ${esc(t.modxna.version)} ` +
    `(@${esc((t.modxna.commit || '').slice(0, 7))}). Names and notes from the modXNA catalogue page, read on ` +
    `${esc(t.names_source.read)}. Flagged entries: ${esc(t.review_source)} ` +
    `${t.prebuilt.length} residues prebuilt; unavailable ones are listed with the reason ` +
    '(missing modXNA parameters, duplicate atom names, a sugar or base that a fiber scaffold cannot carry).';
  renderLibrary();
  renderComposer();
}
function writeCell(w) {
  // "* or *1" -> code spans; descriptive values ("default linkage") stay text
  return w.split(' or ').map((x) => (x.includes(' ') ? esc(x) : `<code>${esc(x)}</code>`)).join(' or ');
}
function noteHTML(n) {
  return /^https?:\/\//.test(n) ? `<a href="${esc(n)}" target="_blank" rel="noopener">${esc(n)}</a>` : esc(n);
}
const FLAG = { name: '⚠ catalogue error: ', duplicate: '⚠ ', protected: 'protected: ', stereo: 'ⓘ ' };
function renderLibrary() {
  const t = state.fragments;
  if (!t) return;
  const q = $('#libfilter').value.trim().toLowerCase();
  const only = $('#libavail').checked;
  const match = (e) => (!only || e.available) &&
    (!q || [e.codes.join('/'), e.name, ...e.notes].join(' ').toLowerCase().includes(q));
  let html = '';
  const jump = [];
  t.sections.forEach((s, i) => {
    const rows = s.entries.filter(match);
    if (!rows.length) return;
    const n = rows.reduce((k, e) => k + e.codes.length, 0);
    const avail = rows.filter((e) => e.available).reduce((k, e) => k + e.codes.length, 0);
    const short = s.moiety === 'base' ? (s.parent === 'other' ? 'other bases' : `${s.parent}-based`) : s.title.split(' ')[0].toLowerCase();
    jump.push(`<button type="button" data-libjump="${i}">${esc(short)} <span>${n}</span></button>`);
    html += `<h3 id="libsec-${i}">${esc(s.title)} <span class="count">${n} fragment${n === 1 ? '' : 's'}` +
      `${avail ? ` · ${avail} prebuilt` : ''}</span></h3>` +
      '<div class="flowtable"><table class="libtable"><thead><tr><th>code</th><th>name</th><th>notes</th><th>in dsXNAbuilder</th></tr></thead><tbody>' +
      rows.map((e) => `<tr${e.available ? ' class="avail"' : ''}>
        <td class="code">${e.codes.map(esc).join('<br>')}</td>
        <td class="name">${esc(e.name)}${e.flags.map((f) => `<div class="flag ${esc(f.kind)}">${FLAG[f.kind] || ''}${esc(f.text)}</div>`).join('')}</td>
        <td class="desc">${e.notes.map(noteHTML).join('<br>')}</td>
        <td class="use">${e.available ? `<span class="badge avail">✓ prebuilt</span> ${writeCell(e.write || '')}` : ''}</td></tr>`).join('') +
      '</tbody></table></div>';
  });
  $('#libtables').innerHTML = html || '<p class="hint">No fragment matches.</p>';
  $('#libjump').innerHTML = jump.length > 1 ? 'Jump to ' + jump.join('') : '';
  for (const b of $$('[data-libjump]')) b.onclick = () => $(`#libsec-${b.dataset.libjump}`).scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// ------------------------------------------------------------------ residue names
// Modified residues get X01, X02, ... by default; the user may rename them. The names go
// to Python with the request (validated there) and into the PDB, the package and the URL.
const CAPS = { int: 'internal', '5cap': '5′ end', '3cap': '3′ end' };
function renderNames(out) {
  const named = (out.summary && out.summary.named) || [];
  $('#names').classList.toggle('hidden', !named.length);
  $('#nametable').innerHTML = '<thead><tr><th>name</th><th>modXNA residue</th><th>position</th><th>what it is</th></tr></thead><tbody>' +
    named.map((n) => `<tr><td><input type="text" class="resname" maxlength="3" spellcheck="false"
      data-key="${esc(n.key)}" data-auto="${n.custom ? '' : esc(n.name)}" value="${esc(n.name)}" aria-label="name of ${esc(n.key)}"></td>
      <td><code>${esc(n.triplet.join('/'))}</code></td><td>${CAPS[n.cap] || esc(n.cap)}</td><td class="desc">${esc(n.description)}</td></tr>`).join('') +
    '</tbody>';
}
function applyNames() {
  const names = {};
  for (const i of $$('#nametable input.resname')) {
    const v = i.value.trim().toUpperCase();
    // keep names the user typed or chose before; automatic ones stay automatic
    if (v && v !== i.dataset.auto) names[i.dataset.key] = v;
  }
  state.names = names;
  build();
}

// ------------------------------------------------------------------ residue composer
// Click a backbone, a sugar and a base to write {BACKBONE/SUGAR/BASE:LETTER} into a strand.
// Only prebuilt residues can be built here; the status line says, for each position in a
// strand, whether this exact combination is prebuilt, unavailable (and why) or not built.
const comp = { backbone: 'DPO', sugar: 'DC2', base: 'DAA', family: 'A', letter: 'A' };
const caret = { top: null, bottom: null };      // last caret position in each strand
const LETTERS = ['A', 'C', 'G', 'T', 'U'];
const FAMILIES = ['A', 'C', 'G', 'T', 'U', 'other'];
const CANONICAL_LETTER = { DAA: 'A', RAA: 'A', DCC: 'C', RCC: 'C', DGG: 'G', RGG: 'G', DTT: 'T', RUU: 'U', M5C: 'C' };
function fragIndex() {
  if (!state.fragIdx) {
    const idx = new Map();
    for (const s of state.fragments.sections) {
      for (const e of s.entries) for (const c of e.codes) idx.set(c, { moiety: s.moiety, family: s.parent, e });
    }
    state.fragIdx = idx;
    state.prebuilt = new Set(state.fragments.prebuilt);
  }
  return state.fragIdx;
}
function chip(code, moiety) {
  const f = fragIndex().get(code);
  const title = `${code} — ${f.e.name}` + (f.e.available ? '' : ` (not available here: ${f.e.reason})`);
  return `<button type="button" class="chip-btn${f.e.available ? '' : ' unavail'}" data-${moiety}="${esc(code)}"
    aria-pressed="${comp[moiety] === code}" title="${esc(title)}">${esc(code)}</button>`;
}
function fragInfo(code) {
  const f = fragIndex().get(code);
  if (!f) return '';
  const flags = f.e.flags.map((x) => `<div class="flag ${esc(x.kind)}">${FLAG[x.kind] || ''}${esc(x.text)}</div>`).join('');
  const notes = f.e.notes.length ? ` · <span class="muted">${f.e.notes.map(esc).join(' · ')}</span>` : '';
  const why = f.e.available ? '' : `<div class="flag name">✗ not available here: ${esc(f.e.reason)}</div>`;
  return `<b>${esc(code)}</b> ${esc(f.e.name)}${notes}${flags}${why}`;
}
function renderComposer() {
  if (!state.fragments) return;
  const idx = fragIndex();
  const all = $('#c-all').checked;
  const codes = (moiety, pred = () => true) => [...idx.entries()]
    .filter(([c, f]) => f.moiety === moiety && c !== '5PO' && pred(f) && (all || f.e.available || comp[moiety] === c))
    .map(([c]) => c);
  $('#c-backbone').innerHTML = codes('backbone').map((c) => chip(c, 'backbone')).join('');
  $('#c-sugar').innerHTML = codes('sugar').map((c) => chip(c, 'sugar')).join('');
  $('#c-family').innerHTML = FAMILIES.map((f) => {
    const n = codes('base', (x) => x.family === f).length;
    return `<button type="button" data-family="${f}" aria-pressed="${comp.family === f}"${n ? '' : ' disabled'}>${f === 'other' ? 'other' : f + '*'} <span>${n}</span></button>`;
  }).join('');
  $('#c-base').innerHTML = codes('base', (x) => x.family === comp.family).map((c) => chip(c, 'base')).join('');
  $('#c-letter').innerHTML = LETTERS.map((l) => `<button type="button" data-letter="${l}" aria-pressed="${comp.letter === l}">${l}</button>`).join('');
  for (const m of ['backbone', 'sugar', 'base']) $(`#i-${m}`).innerHTML = fragInfo(comp[m]);
  const fam = idx.get(comp.base) ? idx.get(comp.base).family : null;
  $('#c-letterhint').textContent = CANONICAL_LETTER[comp.base] ? '' : fam === 'other'
    ? 'choose how this base pairs: there is no canonical parent'
    : `taken from its modXNA family (${fam}-based); change it if the base pairs otherwise`;
  $('#c-code').textContent = composedCode();
  $('#c-status').innerHTML = composedStatus();
}
function composedCode() {
  return `{${comp.backbone}/${comp.sugar}/${comp.base}:${comp.letter}}`;
}
function composedStatus() {
  const un = state.fragments.unavailable || {};
  const canonical = ['DPO', 'RPO'].includes(comp.backbone) &&
    ((comp.sugar === 'DC2' && ['DAA', 'DCC', 'DGG', 'DTT'].includes(comp.base)) ||
     (comp.sugar === 'RC3' && ['RAA', 'RCC', 'RGG', 'RUU'].includes(comp.base)));
  const pos = [['internal', `${comp.backbone}-${comp.sugar}-${comp.base}-int`],
    ['3′ end', `${comp.backbone}-${comp.sugar}-${comp.base}-3cap`],
    ['5′ end (5′-OH, no phosphate)', `5PO-${comp.sugar}-${comp.base}-5cap`]];
  return pos.map(([label, key]) => {
    let s;
    if (canonical && label === 'internal') s = '<span class="ok">✓ canonical (Amber)</span>';
    else if (state.prebuilt.has(key)) s = '<span class="ok">✓ prebuilt</span>';
    else if (un[key] || un[`${comp.backbone}-*-*-*`]) s = `<span class="bad">✗ ${esc(un[key] || un[`${comp.backbone}-*-*-*`])}</span>`;
    else s = `<span class="muted">– not prebuilt: this combination is not in the library (<code>tools/build_residue_library.py --only ${esc(key)}</code>)</span>`;
    return `<div><span class="pos">${label}</span> ${s}</div>`;
  }).join('');
}
function insertComposed(which) {
  const input = $(which === 'top' ? '#top' : '#bottom');
  const code = composedCode();
  const at = caret[which] ?? input.value.length;
  input.value = input.value.slice(0, at) + code + input.value.slice(at);
  caret[which] = at + code.length;
  selectTab('build');
  input.focus();
  input.setSelectionRange(caret[which], caret[which]);
  validateSoon(0);
}
function composerClick(ev) {
  const b = ev.target.closest('button');
  if (!b) return;
  if (b.dataset.backbone) comp.backbone = b.dataset.backbone;
  else if (b.dataset.sugar) comp.sugar = b.dataset.sugar;
  else if (b.dataset.family) comp.family = b.dataset.family;
  else if (b.dataset.letter) comp.letter = b.dataset.letter;
  else if (b.dataset.base) {
    comp.base = b.dataset.base;
    const fam = fragIndex().get(comp.base).family;
    comp.letter = CANONICAL_LETTER[comp.base] || (LETTERS.includes(fam) ? fam : comp.letter);
  } else return;
  renderComposer();
}

// ------------------------------------------------------------------ request
function request() {
  return {
    top: $('#top').value, bottom: $('#bottom').value, model: +$('#model').value || 4,
    names: { ...state.names },
    amber: { dna_ff: $('#dnaff').value, water: $('#water').value, box: $('#box').value, buffer: +$('#buffer').value || 10 },
  };
}
function validateSoon(delay = 250) {
  clearTimeout(state.validateTimer);
  state.validateTimer = setTimeout(validate, delay);
}
async function validate() {
  if (!state.ready) return;
  const req = request();
  $('#boxrow').classList.toggle('hidden', req.amber.water === 'none');
  if (!req.top.trim() || !req.bottom.trim()) {
    $('#grid').className = 'grid empty-grid';
    $('#grid').innerHTML = '<span class="hint">Type both strands to see the pairing.</span>';
    $('#issues').innerHTML = '';
    $('#build').disabled = true;
    return;
  }
  const mySeq = ++state.validSeq;
  try {
    const v = await call('validate', req);
    if (mySeq !== state.validSeq) return;
    renderValidation(v);
  } catch (e) {
    $('#issues').innerHTML = `<li class="error">${esc(e.message)}</li>`;
  }
}
function renderValidation(v) {
  $('#build').disabled = !v.ok;
  if (v.columns) {
    const sym = { wc: '|', wobble: ':', mismatch: '×' };
    const cols = v.columns.map((c) => `<div class="col">
      <span class="tok${tokClass(c.top, c.top_mod)}">${esc(c.top || '–')}</span>
      <span class="pair ${c.pair || ''}">${c.pair ? sym[c.pair] : ''}</span>
      <span class="tok${tokClass(c.bottom, c.bottom_mod)}">${esc(c.bottom || '–')}</span>
      <span class="num">${(c.i + 1) % 5 === 0 || c.i === 0 ? c.i + 1 : ''}</span></div>`).join('');
    $('#grid').className = 'grid';
    $('#grid').innerHTML = `<div class="end"><div>5′</div><div></div><div>3′</div></div>${cols}<div class="end"><div>3′</div><div></div><div>5′</div></div>`;
  }
  renderIssues(v.issues || []);
}
function tokClass(tok, mod) {
  if (!tok) return ' gap';
  if (tok === '-') return ' abasic';
  return mod ? ' mod' : '';
}
function renderIssues(list) {
  const order = { error: 0, warning: 1, info: 2 };
  $('#issues').innerHTML = [...list].sort((a, b) => order[a.level] - order[b.level]).map((i) =>
    `<li class="${i.level}">${esc(i.message)}<span class="code">${esc(i.code)}</span></li>`).join('');
}

// ------------------------------------------------------------------ build
async function build() {
  if (!state.ready) return;
  const mine = ++state.buildSeq;      // a reset meanwhile discards this build
  const req = request();
  const t0 = performance.now();
  $('#build').disabled = true;
  $('#status').innerHTML = '<span class="spinner"></span>Building…';
  try {
    const v = await call('validate', req);
    if (mine !== state.buildSeq) return;
    renderValidation(v);
    if (!v.ok) return;
    req.libs = await fetchLibs(v.needed_libs || []);
    if ((v.needed_libs || []).length) req.frcmod = await fetchText('lib/frcmod.modxna');
    const out = await call('build', req);
    if (mine !== state.buildSeq) return;
    renderIssues(out.issues || []);
    if (!out.ok) return;
    state.last = out;
    renderSummary(out);
    renderNames(out);
    // wall time from click to result (validation, library fetch, build); read by
    // scripts/bench_browser.py and handy when reporting performance
    document.body.dataset.buildMs = String(Math.round(performance.now() - t0));
    state.lastReq = req;
    try {
      await showStructure(out.pdb, out.summary);
    } catch (e) {
      $('#viewer').innerHTML = '<div class="empty">3D view unavailable (this browser has no WebGL). The structure and the package can still be downloaded.</div>';
      state.stage = null;
    }
    for (const id of ['#dlzip', '#dlpdb', '#link']) $(id).disabled = false;
    writeURL(req);
  } catch (e) {
    $('#issues').innerHTML = `<li class="error">${esc(e.message)}</li>`;
  } finally {
    $('#status').textContent = '';
    $('#build').disabled = false;
  }
}
async function fetchText(path) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
  return r.text();
}
async function fetchLibs(keys) {
  const texts = await Promise.all(keys.map((k) => fetchText(`lib/residues/${encodeURIComponent(k)}.lib`)));
  return Object.fromEntries(keys.map((k, i) => [k, texts[i]]));
}
function renderSummary(out) {
  const s = out.summary;
  const stat = (v, l) => `<div class="stat"><div class="v">${esc(v)}</div><div class="l">${esc(l)}</div></div>`;
  $('#summary').classList.remove('hidden');
  $('#summary').innerHTML = stat(s.n_bp, 'base pairs') + stat(s.n_atoms, `atoms (${s.n_heavy} heavy)`) +
    stat(`#${s.model}`, `model · ${s.family}`) +
    stat(s.o3p ? `${s.o3p[0].toFixed(2)}–${s.o3p[1].toFixed(2)}` : '–', 'O3′–P (Å)') +
    stat(s.clashes, 'heavy-atom contacts < 2.2 Å') + (s.n_modified ? stat(s.n_modified, 'modXNA residues') : '');
}

// ------------------------------------------------------------------ viewer
// Colours. "base" paints the base atoms (and their hydrogens) by base, so that
// Watson-Crick pairs read red-blue (A-T/U) and green-amber (G-C) and a mismatch
// stands out; the sugar-phosphate stays grey. "modified" paints whole residues:
// canonical grey, modified (modXNA chemistry) red, abasic amber. "element" is
// CPK with carbon and hydrogen darkened on the light background.
const BASE_COLORS = { A: 0xe5484d, T: 0x3e7bfa, U: 0x36b3d9, G: 0x30a46c, C: 0xf5a524 };
const OTHER_BASE = 0xa855f7;
const BACKBONE = 0x8c8c8c;
const KIND_COLORS = { canonical: 0x8c8c8c, modified: 0xdc2626, abasic: 0xeab308 };
// atoms without a prime are base atoms, except the phosphate, the
// phosphorothioate sulfurs and the 2'-O-methyl carbon of modXNA (CM2)
const NOT_BASE = new Set(['P', 'OP1', 'OP2', 'OP3', 'HOP3', 'S1', 'S2', 'CM2']);
const isBase = (name) => !name.includes("'") && !NOT_BASE.has(name);
const ELEMENT_COLORS = {
  dark: { C: 0x909090, H: 0xf0f0f0, N: 0x3050f8, O: 0xff0d0d, P: 0xff8000, S: 0xe6c619, F: 0x90e050 },
  light: { C: 0x6e6e6e, H: 0xb4b4b4, N: 0x3050f8, O: 0xe00000, P: 0xff8000, S: 0xc9a800, F: 0x6cc040 },
};
const SCHEMES = {};
const theme = () => (document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark');
const residueOf = (atom) => state.residues.get(`${atom.chainname}:${atom.resno}`);
function registerSchemes() {
  if (SCHEMES.base) return;         // once per page: the viewer may be rebuilt after a reset
  const reg = (label, color) => NGL.ColormakerRegistry.addScheme(function () { this.atomColor = color; }, label);
  SCHEMES.base = reg('dsxna-base', (atom) => {
    const r = residueOf(atom);
    if (!r || !isBase(atom.atomname)) return BACKBONE;
    return BASE_COLORS[r.letter] ?? OTHER_BASE;
  });
  SCHEMES.modified = reg('dsxna-modified', (atom) => {
    const r = residueOf(atom);
    return KIND_COLORS[r ? r.kind : 'canonical'];
  });
  SCHEMES.element = reg('dsxna-element', (atom) => {
    const c = ELEMENT_COLORS[theme()];
    return c[atom.element] ?? c.C;
  });
}
function legend() {
  const sw = (hex, label) => `<span class="sw"><i style="background:#${hex.toString(16).padStart(6, '0')}"></i>${esc(label)}</span>`;
  const kinds = new Set([...state.residues.values()].map((r) => r.kind));
  const letters = new Set([...state.residues.values()].map((r) => r.letter));
  let html = '';
  if (state.color === 'base') {
    html = ['A', 'C', 'G', 'T', 'U'].filter((l) => letters.has(l)).map((l) => sw(BASE_COLORS[l], l)).join('') +
      ([...letters].some((l) => l !== '-' && !(l in BASE_COLORS)) ? sw(OTHER_BASE, 'other base') : '') +
      sw(BACKBONE, 'sugar–phosphate');
  } else if (state.color === 'modified') {
    html = ['canonical', 'modified', 'abasic'].filter((k) => kinds.has(k)).map((k) => sw(KIND_COLORS[k], k)).join('');
  } else {
    const c = ELEMENT_COLORS[theme()];
    html = ['C', 'H', 'N', 'O', 'P', 'S', 'F'].filter((e) => e === 'C' || e === 'H' || hasElement(e)).map((e) => sw(c[e], e)).join('');
  }
  $('#legend').innerHTML = html;
}
function hasElement(e) {
  return !!state.last && state.last.pdb.split('\n').some((l) => l.startsWith('ATOM') && l.slice(76, 78).trim() === e);
}
function bg() { return theme() === 'light' ? '#f2f2f2' : '#0b0b0b'; }
function hasWebGL() {
  try {
    const c = document.createElement('canvas');
    return !!(c.getContext('webgl2') || c.getContext('webgl'));
  } catch (e) { return false; }
}
async function showStructure(pdb, summary) {
  if (!window.NGL || !hasWebGL()) throw new Error('no WebGL');
  if (!state.stage) {
    $('#viewer').innerHTML = '';
    registerSchemes();
    state.stage = new NGL.Stage('viewer', { backgroundColor: bg(), tooltip: true });
    window.addEventListener('resize', () => state.stage && state.stage.handleResize());
  }
  state.residues = new Map((summary.residues || []).map(([chain, resseq, letter, kind]) => [`${chain}:${resseq}`, { letter, kind }]));
  state.stage.removeAllComponents();
  state.comp = await state.stage.loadFile(new Blob([pdb], { type: 'text/plain' }), { ext: 'pdb', name: 'dsxnabuilder' });
  state.modSele = Object.values(summary.name_map || {}).map((n) => `[${n}]`).join(' or ');
  applyRepr();
  state.comp.autoView();
}
function applyRepr() {
  legend();
  const c = state.comp;
  if (!c) return;
  c.removeAllRepresentations();
  const colorScheme = SCHEMES[state.color];
  if (state.repr === 'cartoon') {
    c.addRepresentation('cartoon', { colorScheme: 'chainid', aspectRatio: 4, radiusScale: 1.2 });
    c.addRepresentation('base', { colorScheme });
    if (state.modSele) c.addRepresentation('ball+stick', { sele: state.modSele, colorScheme });
  } else {
    c.addRepresentation('licorice', { colorScheme });
  }
}

// ------------------------------------------------------------------ downloads, URL
function download(name, data, type) {
  const url = URL.createObjectURL(new Blob([data], { type }));
  const a = Object.assign(document.createElement('a'), { href: url, download: name });
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
function b64bytes(b64) {
  const s = atob(b64);
  const u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
  return u;
}
function writeURL(req) {
  const p = new URLSearchParams({ top: req.top, bottom: req.bottom, model: req.model,
    ff: req.amber.dna_ff, water: req.amber.water, box: req.amber.box, buffer: req.amber.buffer });
  const names = Object.entries(req.names || {});
  if (names.length) p.set('names', names.map(([k, v]) => `${k}=${v}`).join(';'));
  history.replaceState(null, '', '?' + p.toString());
}
function applyURL() {
  const p = new URLSearchParams(location.search);
  if (p.get('tab')) selectTab(p.get('tab'));
  const set = (id, k) => { if (p.get(k) !== null) $(id).value = p.get(k); };
  set('#top', 'top'); set('#bottom', 'bottom'); set('#model', 'model'); set('#dnaff', 'ff');
  set('#water', 'water'); set('#box', 'box'); set('#buffer', 'buffer');
  state.names = Object.fromEntries((p.get('names') || '').split(';').filter((x) => x.includes('='))
    .map((x) => [x.slice(0, x.lastIndexOf('=')), x.slice(x.lastIndexOf('=') + 1)]));
  showModelInfo();
  if (p.get('top') && p.get('bottom')) build(); else validateSoon(0);
}

// ------------------------------------------------------------------ tabs, examples
function selectTab(name) {
  if (name === 'library') loadLibrary();
  for (const tab of ['build', 'models', 'library', 'docs']) {
    $(`#tab-${tab}`).setAttribute('aria-selected', String(tab === name));
    $(`#view-${tab}`).classList.toggle('hidden', tab !== name);
  }
  if (name === 'build' && state.stage) setTimeout(() => state.stage.handleResize(), 0);
}
function renderExamples() {
  $('#examples').innerHTML = EXAMPLES.map((e) => `<button class="btn ghost small" type="button" data-ex="${e.id}">${esc(e.label)}</button>`).join('');
  for (const b of $$('[data-ex]')) b.onclick = () => {
    const e = EXAMPLES.find((x) => x.id === b.dataset.ex);
    $('#top').value = e.top; $('#bottom').value = e.bottom; $('#model').value = e.model;
    state.names = {};
    showModelInfo(); build();
  };
}

// ------------------------------------------------------------------ reset
// Back to the starting state, as on first load: empty strands, B-DNA #4, default Amber
// options and viewer settings, no structure, no names, a clean URL.
function resetAll() {
  clearTimeout(state.validateTimer);
  state.validSeq++;                 // results still on their way are dropped
  state.buildSeq++;
  $('#top').value = '';
  $('#bottom').value = '';
  $('#model').value = 4;
  showModelInfo();
  $('#dnaff').value = 'OL15';
  $('#water').value = 'none';
  $('#box').value = 'octahedral';
  $('#buffer').value = 10;
  $('#advanced').open = false;
  state.names = {};
  state.last = state.lastReq = state.comp = null;
  state.residues = new Map();
  state.spinning = false;
  if (state.stage) {
    state.stage.dispose();
    state.stage = null;
  }
  $('#viewer').innerHTML = '<div class="empty">The built duplex appears here.</div>';
  $('#legend').innerHTML = '';
  $('#summary').classList.add('hidden');
  $('#summary').innerHTML = '';
  $('#names').classList.add('hidden');
  $('#nametable').innerHTML = '';
  $('#issues').innerHTML = '';
  $('#status').textContent = '';
  for (const id of ['#dlzip', '#dlpdb', '#link']) $(id).disabled = true;
  state.repr = 'sticks';
  state.color = 'base';
  for (const b of $$('#repr button, #color button')) {
    b.setAttribute('aria-pressed', String(b.dataset.v === state[b.parentElement.id]));
  }
  Object.assign(comp, { backbone: 'DPO', sugar: 'DC2', base: 'DAA', family: 'A', letter: 'A' });
  caret.top = caret.bottom = null;
  $('#c-all').checked = false;
  $('#libfilter').value = '';
  $('#libavail').checked = false;
  $('#libdetails').open = false;
  renderComposer();
  renderLibrary();
  history.replaceState(null, '', location.pathname);
  selectTab('build');
  validateSoon(0);
}

// ------------------------------------------------------------------ wiring
for (const tab of ['build', 'models', 'library', 'docs']) $(`#tab-${tab}`).onclick = () => selectTab(tab);
$('#libfilter').addEventListener('input', renderLibrary);
$('#composer').addEventListener('click', composerClick);
$('#c-all').addEventListener('change', renderComposer);
$('#c-top').onclick = () => insertComposed('top');
$('#c-bottom').onclick = () => insertComposed('bottom');
$('#c-copy').onclick = async () => { try { await navigator.clipboard.writeText(composedCode()); } catch (e) { prompt('Code', composedCode()); } };
$('#gocompose').onclick = () => selectTab('library');
$('#applynames').onclick = applyNames;
$('#resetnames').onclick = () => { state.names = {}; build(); };
for (const which of ['top', 'bottom']) {
  const el = $(`#${which}`);
  const keep = () => { caret[which] = el.selectionStart; };
  for (const ev of ['keyup', 'click', 'select', 'input']) el.addEventListener(ev, keep);
}
$('#libavail').addEventListener('change', renderLibrary);
for (const id of ['#top', '#bottom']) $(id).addEventListener('input', () => validateSoon());
for (const id of ['#model', '#water', '#dnaff', '#box', '#buffer']) $(id).addEventListener('change', () => { showModelInfo(); validateSoon(0); });
$('#build').onclick = build;
$('#reset').onclick = resetAll;
$('#complement').onclick = async () => {
  try {
    const r = await call('complement', { top: $('#top').value, model: +$('#model').value });
    $('#bottom').value = r.bottom;
    validateSoon(0);
  } catch (e) { $('#issues').innerHTML = `<li class="error">${esc(e.message)}</li>`; }
};
$('#clear').onclick = () => { $('#top').value = ''; $('#bottom').value = ''; validateSoon(0); };
function segmented(sel, key) {
  for (const b of $$(`${sel} button`)) b.onclick = () => {
    state[key] = b.dataset.v;
    for (const x of $$(`${sel} button`)) x.setAttribute('aria-pressed', String(x === b));
    applyRepr();
  };
}
segmented('#repr', 'repr');
segmented('#color', 'color');
$('#spin').onclick = () => { if (state.stage) { state.spinning = !state.spinning; state.stage.setSpin(state.spinning); } };
$('#recenter').onclick = () => { if (state.comp) state.comp.autoView(); };
$('#dlzip').onclick = async () => {
  // the Amber/modXNA package is an advanced option: built on demand with the
  // current options of the Advanced panel, for the last built duplex
  if (!state.last || !state.lastReq) return;
  const req = { ...state.lastReq, amber: request().amber, package: true };
  $('#status').innerHTML = '<span class="spinner"></span>Building…';
  try {
    const v = await call('validate', req);
    req.libs = await fetchLibs(v.needed_libs || []);
    if ((v.needed_libs || []).length) req.frcmod = await fetchText('lib/frcmod.modxna');
    const out = await call('build', req);
    if (out.ok && out.zip_b64) download('dsxnabuilder-amber.zip', b64bytes(out.zip_b64), 'application/zip');
    else renderIssues(out.issues || []);
  } catch (e) {
    $('#issues').innerHTML = `<li class="error">${esc(e.message)}</li>`;
  } finally { $('#status').textContent = ''; }
};
$('#dlpdb').onclick = () => state.last && download('dsxnabuilder.pdb', state.last.pdb, 'chemical/x-pdb');
$('#link').onclick = async () => {
  try { await navigator.clipboard.writeText(location.href); $('#status').textContent = 'Link copied'; } catch (e) { prompt('URL', location.href); }
};
// The theme belongs to the surrounding danslab.xyz page. Its toggle writes the shared
// localStorage key, and because the two documents are same-origin the storage event
// reaches us here, so the embedded app follows along without a reload.
addEventListener('storage', (ev) => {
  if (ev.key === 'theme' && (ev.newValue === 'light' || ev.newValue === 'dark')) {
    document.documentElement.setAttribute('data-theme', ev.newValue);
    if (state.stage) { state.stage.setParameters({ backgroundColor: bg() }); applyRepr(); }
  }
});

// ------------------------------------------------------------------ embedding
// On danslab.xyz the app runs in a same-origin iframe (components/AppFrame.tsx). Reporting
// our height lets the host grow the frame, so the page has one scrollbar instead of two.
// The height is that of the <html> box, not scrollHeight: scrollHeight never drops below
// the frame's current height, so the frame could grow (Documentation) but never shrink back.
if (window.parent !== window) {
  document.documentElement.classList.add('embedded');
  let last = 0;
  const report = () => {
    const h = Math.ceil(document.documentElement.getBoundingClientRect().height);
    if (Math.abs(h - last) < 2) return;
    last = h;
    parent.postMessage({ type: 'app:height', height: h }, location.origin);
  };
  if (window.ResizeObserver) new ResizeObserver(report).observe(document.documentElement);
  addEventListener('load', report);
  setInterval(report, 1000);      // catches growth the observer misses (fonts)
  report();
}

renderExamples();
