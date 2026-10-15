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
  // explicit modXNA triplet {BACKBONE/SUGAR/BASE:pairing letter}: phosphorothioate PS2,
  // 2'-O-methylribose and 5-methylcytosine, with an abasic site and overhangs
  { id: 'modxna', label: 'modXNA code', top: 'CG{PS2/OME/M5C:C}GAA*2TTC-CG--', bottom: 'GCGCTT*1AAGCGCAA', model: 4 },
];

const STATUS = {
  loading_python: 'Loading Python runtime…', loading_numpy: 'Loading numpy…', loading_app: 'Loading dsXNAbuilder…',
};
// the loading animation (a turning double helix), as in index.html
const HELIX = '<div class="helix" aria-hidden="true">' +
  Array.from({ length: 10 }, (_, i) => `<span style="--i:${i}"><b></b></span>`).join('') + '</div>';
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
// The deployment stamps PAGE_VERSION and the ?v= of every asset with the version it
// ships (scripts/deploy_danslab.sh); locally they stay "dev".
const PAGE_VERSION = 'dev';
const worker = new Worker(new URL('worker.js?v=dev', import.meta.url), { type: 'module' });
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
    // GitHub Pages lets browsers keep a page for 10 minutes; right after a deployment a
    // cached page can meet the new Python package (the worker reads VERSION afresh). Reload
    // once to get the matching page rather than run the two together.
    const stale = PAGE_VERSION !== 'dev' && d.version !== PAGE_VERSION;
    if (stale) {
      let tried = false;
      try { tried = sessionStorage.getItem('dsxna-reloaded') === d.version; sessionStorage.setItem('dsxna-reloaded', d.version); } catch (e) {}
      if (!tried) { location.reload(); return; }
    }
    state.ready = true;
    $('#ver').textContent = 'v' + d.version;
    $('#status').innerHTML = stale
      ? `<span class="error">A new version (${esc(d.version)}) was published: reload the page (Ctrl+Shift+R).</span>`
      : '';
    call('models').then((info) => { setModels(info); applyURL(); bootDone(); }).catch((e) => {
      $('#bootmsg').textContent = `Could not start the engine: ${e.message}`;
      $('#boot').classList.add('failed');
    });
    if (!$('#view-library').classList.contains('hidden')) loadLibrary();
    if ($('#batchtext').value.trim()) batchValidateSoon(0);
  } else if (d.status === 'error') {
    $('#bootmsg').textContent = `Could not start the engine: ${d.error}`;
    $('#boot').classList.add('failed');
    $('#status').innerHTML = `<span class="error">Could not start the engine: ${esc(d.error)}</span>`;
  } else {
    $('#bootmsg').textContent = STATUS[d.status] || d.status;
  }
}
function bootDone() {
  const boot = $('#boot');
  boot.classList.add('done');
  setTimeout(() => boot.classList.add('hidden'), 400);
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
// The sugars, bases and linkages of the prebuilt library: every combination exists, so the
// tab only offers what can be built, in the short notation of the sequence language.
async function loadLibrary() {
  if (state.pieces || !state.ready) return;
  try {
    state.pieces = await call('pieces');
  } catch (e) {
    $('#piecetable').innerHTML = `<caption class="error">${esc(e.message)}</caption>`;
    return;
  }
  const t = state.pieces;
  const also = (x) => (x.also.length ? ` (also ${x.also.map((a) => `<code>${esc(a)}</code>`).join(', ')})` : '');
  $('#piecetable').innerHTML = '<thead><tr><th></th><th>write</th><th>what it is</th></tr></thead><tbody>' +
    t.sugars.map((x, i) => `<tr><td class="kind">${i ? '' : 'sugar'}</td><td><code>${esc(x.prefix)}</code>${also(x)}</td><td>${esc(x.name)}</td></tr>`).join('') +
    t.bases.map((x, i) => `<tr><td class="kind">${i ? '' : 'base'}</td><td><code>${esc(x.write)}</code></td><td>${esc(x.name)}</td></tr>`).join('') +
    t.linkages.map((x, i) => `<tr><td class="kind">${i ? '' : 'linkage'}</td><td>${x.write ? `<code>${esc(x.write)}</code>` : 'nothing'}</td><td>${esc(x.name)}</td></tr>`).join('') +
    '</tbody>';
  $('#libsource').textContent = `${t.n_residues} residues prebuilt with modXNA ${t.modxna.version} ` +
    `(@${(t.modxna.commit || '').slice(0, 7)}). Other modXNA fragments are not prebuilt yet.`;
  renderComposer();
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
// Click a sugar, a base and a linkage; the nucleotide is written in the short notation
// (e.g. *1+C^5mC) at the cursor of either strand.
const comp = { sugar: 0, base: 0, link: 0 };       // indices into state.pieces
const caret = { top: null, bottom: null };      // last caret position in each strand
function composed(which) {
  const p = state.pieces;
  const nt = p.sugars[comp.sugar].prefix + p.bases[comp.base].write;
  const link = p.linkages[comp.link].write;
  // the linkage is the nucleotide's 5' phosphate: before it in strand 1, after it in
  // strand 2 (typed 3'->5')
  return which === 'top' ? link + nt : nt + link;
}
function renderComposer() {
  const p = state.pieces;
  if (!p) return;
  const chips = (list, key, label) => list.map((x, i) =>
    `<button type="button" class="chip-btn" data-${key}="${i}" aria-pressed="${comp[key] === i}" title="${esc(x.name)}">${esc(label(x))}</button>`).join('');
  $('#c-sugar').innerHTML = chips(p.sugars, 'sugar', (x) => x.prefix);
  $('#c-base').innerHTML = chips(p.bases, 'base', (x) => x.write.replace('C^', ''));
  $('#c-link').innerHTML = chips(p.linkages, 'link', (x) => x.write || 'PO');
  $('#c-top-code').textContent = composed('top');
  $('#c-bottom-code').textContent = composed('bottom');
  $('#c-info').textContent = `${p.bases[comp.base].name} · ${p.sugars[comp.sugar].name} · ${p.linkages[comp.link].name}`;
}
function insertComposed(which) {
  const input = $(which === 'top' ? '#top' : '#bottom');
  const code = composed(which);
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
  for (const key of ['sugar', 'base', 'link']) if (b.dataset[key] !== undefined) comp[key] = +b.dataset[key];
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
  $('#build').textContent = 'Building…';
  $('#status').textContent = '';
  setBusy(true);
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
    setBusy(false);
    $('#build').textContent = 'Build duplex';
    $('#build').disabled = false;
  }
}
// the helix over the viewer while a build runs (it fades in late, so quick builds do not flash)
function setBusy(on) {
  const old = $('#viewer > .busy');
  if (old) old.remove();
  $('#viewer').classList.toggle('is-busy', on);
  if (on) $('#viewer').insertAdjacentHTML('beforeend', `<div class="busy">${HELIX}<div class="msg">Building…</div></div>`);
}
async function fetchText(path) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
  return r.text();
}
const libCache = new Map();         // residue key -> promise of its OFF text
async function fetchLibs(keys) {
  for (const k of keys) {
    if (!libCache.has(k)) {
      const p = fetchText(`lib/residues/${encodeURIComponent(k)}.lib`);
      p.catch(() => libCache.delete(k));
      libCache.set(k, p);
    }
  }
  const texts = await Promise.all(keys.map((k) => libCache.get(k)));
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
const TABS = ['build', 'batch', 'models', 'library', 'docs'];
function selectTab(name) {
  if (!TABS.includes(name)) return;
  if (name === 'library') loadLibrary();
  for (const tab of TABS) {
    $(`#tab-${tab}`).setAttribute('aria-selected', String(tab === name));
    $(`#view-${tab}`).classList.toggle('hidden', tab !== name);
  }
  // the help needs no engine: readable while it loads
  $('#boot').classList.toggle('aside', name === 'docs');
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

// ------------------------------------------------------------------ high-throughput
// Many duplexes from a FASTA-like text (dsxnabuilder/batch.py). Python checks the text and
// every record; the page then builds the records one by one, so that it can show progress
// and stop, and Python packs the PDB files, a report and the input into one zip.
const BATCH_EXAMPLE = `# dsXNAbuilder high-throughput input: one record per duplex
>name=sequence1; model=20
mG*mA*fC*mU*fA*mG*fC*mU*mG*mAdGdG
C*U*GAUCGAC*U--

>name=sequence2; model=4
mG*mA*fC*mU*fA*mG*fC*mU*mG*mA
C*U*GAUCGAC*U

>dickerson
CGCGAATTCGCG
GCGCTTAAGCGC

>name=lna_gapmer; model=22
r(GCAAUGUACUGAAUGC)
+C*+G*+T*dT*dA*dC*dA*dT*dG*dA*dC*dT*dT*+A*+C*+G
`;
const MAX_FILE = 2e6;               // bytes; 50 duplexes take a few kB
const batch = { timer: null, seq: 0, running: false, stop: false, rows: [], results: {}, zip: null };

function batchValidateSoon(delay = 400) {
  clearTimeout(batch.timer);
  batch.timer = setTimeout(batchValidate, delay);
}
async function batchValidate() {
  if (!state.ready || batch.running) return;
  const text = $('#batchtext').value;
  batch.results = {};
  batch.zip = null;
  $('#batchzip').disabled = true;
  $('#batchstatus').textContent = '';
  if (!text.trim()) {
    batch.rows = [];
    renderBatch({ ok: false, issues: [], records: [] });
    return;
  }
  const mine = ++batch.seq;
  try {
    const v = await call('batch_validate', { text });
    if (mine !== batch.seq) return;
    batch.rows = v.records;
    renderBatch(v);
  } catch (e) {
    $('#batchissues').innerHTML = `<li class="error">${esc(e.message)}</li>`;
  }
}
function renderBatch(v) {
  $('#batchbuild').disabled = !v.ok;
  const order = { error: 0, warning: 1, info: 2 };
  $('#batchissues').innerHTML = [...v.issues].sort((a, b) => order[a.level] - order[b.level]).map((i) =>
    `<li class="${i.level}">${esc(i.message)}<span class="code">${esc(i.code)}</span></li>`).join('');
  const n = v.records.length;
  const bad = v.records.filter((r) => !r.ok).length;
  $('#batchcount').textContent = n
    ? `${n} duplex${n === 1 ? '' : 'es'}${bad ? `, ${bad} with errors (reported, not built)` : ''}` : '';
  renderBatchTable();
}
function batchStatus(r) {
  const out = batch.results[r.name];
  if (out) return out.ok ? 'built' : 'failed';
  if (batch.running && batch.current === r.name) return 'building';
  return r.ok ? 'ready' : 'error';
}
function renderBatchTable() {
  if (!batch.rows.length) {
    $('#batchtable').innerHTML = '<caption class="hint">Paste or open the duplexes to see them listed here.</caption>';
    return;
  }
  const rows = batch.rows.map((r, i) => {
    const out = batch.results[r.name];
    const st = batchStatus(r);
    const issues = (out ? out.issues : r.issues) || [];
    const err = issues.find((x) => x.level === 'error');
    const codes = [...new Set(issues.filter((x) => x.level === 'warning').map((x) => x.code))];
    const size = out && out.ok ? `${out.summary.n_atoms} atoms · ` : '';
    const notes = err ? esc(err.message) : size + codes.map((c) => `<span class="code">${esc(c)}</span>`).join('');
    return `<tr><td class="num">${i + 1}</td><td class="name">${esc(r.name)}</td><td class="num">#${r.model}</td>
      <td class="num">${r.n_bp ?? ''}</td><td><span class="st ${st}">${st}</span></td><td class="desc">${notes}</td>
      <td><button class="btn ghost small" type="button" data-view="${i}"${r.ok ? '' : ' disabled'} title="open in the Build tab">view</button></td></tr>`;
  }).join('');
  $('#batchtable').innerHTML = '<thead><tr><th>#</th><th>name</th><th>model</th><th>bp</th><th>status</th><th>notes</th><th></th></tr></thead>' +
    `<tbody>${rows}</tbody>`;
}
function batchProgress(done, total, what) {
  $('#batchbar').style.width = `${total ? Math.round((100 * done) / total) : 0}%`;
  $('#batchmsg').textContent = what;
}
function batchBusy(on) {
  batch.running = on;
  $('#batchbuild').disabled = on;
  $('#batchstop').hidden = !on;
  $('#batchtext').readOnly = on;
  for (const id of ['#batchexample', '#batchopen', '#batchclear']) $(id).disabled = on;
  $('#batchprogress').classList.toggle('hidden', !on);
}
async function batchBuild() {
  if (!state.ready || batch.running) return;
  const text = $('#batchtext').value;
  const mine = ++batch.seq;          // a reset or a new validation meanwhile discards this run
  const t0 = performance.now();
  batch.stop = false;
  batch.results = {};
  batch.zip = null;
  $('#batchzip').disabled = true;
  $('#batchstatus').textContent = '';
  batchBusy(true);
  try {
    batchProgress(0, 1, 'Checking…');
    const v = await call('batch_validate', { text });
    if (mine !== batch.seq) return;
    batch.rows = v.records;
    renderBatch(v);
    if (!v.ok) return;
    const n = v.records.length;
    for (const [i, r] of v.records.entries()) {
      if (batch.stop) break;
      batch.current = r.name;
      batchProgress(i, n, `Building ${i + 1} of ${n} · ${r.name}`);
      renderBatchTable();
      let out;
      if (!r.ok) {
        out = { ok: false, issues: r.issues };
      } else {
        try {
          out = await call('build', { top: r.top, bottom: r.bottom, model: r.model, libs: await fetchLibs(r.needed_libs) });
        } catch (e) {
          out = { ok: false, issues: [{ level: 'error', code: 'ENGINE', message: e.message }] };
        }
      }
      if (mine !== batch.seq) return;
      batch.results[r.name] = out;
    }
    batch.current = null;
    renderBatchTable();
    batchProgress(n, n, 'Packing the zip…');
    const z = await call('batch_zip', { text, results: batch.results });
    if (mine !== batch.seq) return;
    batch.zip = b64bytes(z.zip_b64);
    $('#batchzip').disabled = false;
    const secs = ((performance.now() - t0) / 1000).toFixed(1);
    $('#batchstatus').textContent = `${z.built} of ${n} built${z.failed ? `, ${z.failed} not built` : ''}` +
      `${batch.stop ? ' (stopped)' : ''} in ${secs} s`;
    download('dsxnabuilder_batch.zip', batch.zip, 'application/zip');
  } catch (e) {
    $('#batchissues').innerHTML = `<li class="error">${esc(e.message)}</li>`;
  } finally {
    batch.current = null;
    if (mine === batch.seq) batchBusy(false);
  }
}
function viewRecord(i) {
  const r = batch.rows[i];
  if (!r) return;
  $('#top').value = r.top;
  $('#bottom').value = r.bottom;
  $('#model').value = r.model;
  state.names = {};
  showModelInfo();
  selectTab('build');
  build();
}
async function openBatchFile(file) {
  if (!file) return;
  if (file.size > MAX_FILE) {
    $('#batchissues').innerHTML = `<li class="error">${esc(file.name)} is too large (${(file.size / 1e6).toFixed(1)} MB) for a list of duplexes</li>`;
    return;
  }
  $('#batchtext').value = await file.text();
  batchValidateSoon(0);
}
function resetBatch() {
  clearTimeout(batch.timer);
  batch.seq++;
  batch.stop = true;
  batchBusy(false);
  batch.rows = [];
  batch.results = {};
  batch.zip = null;
  $('#batchtext').value = '';
  $('#batchfile').value = '';
  $('#batchissues').innerHTML = '';
  $('#batchcount').textContent = '';
  $('#batchstatus').textContent = '';
  $('#batchzip').disabled = true;
  $('#batchbuild').disabled = true;
  renderBatchTable();
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
  Object.assign(comp, { sugar: 0, base: 0, link: 0 });
  caret.top = caret.bottom = null;
  renderComposer();
  resetBatch();
  history.replaceState(null, '', location.pathname);
  selectTab('build');
  validateSoon(0);
}

// ------------------------------------------------------------------ wiring
for (const tab of TABS) $(`#tab-${tab}`).onclick = () => selectTab(tab);
$('#composer').addEventListener('click', composerClick);
$('#c-top').onclick = () => insertComposed('top');
$('#c-bottom').onclick = () => insertComposed('bottom');
$('#gocompose').onclick = () => selectTab('library');
$('#applynames').onclick = applyNames;
$('#resetnames').onclick = () => { state.names = {}; build(); };
for (const which of ['top', 'bottom']) {
  const el = $(`#${which}`);
  const keep = () => { caret[which] = el.selectionStart; };
  for (const ev of ['keyup', 'click', 'select', 'input']) el.addEventListener(ev, keep);
}
for (const id of ['#top', '#bottom']) $(id).addEventListener('input', () => validateSoon());
for (const id of ['#model', '#water', '#dnaff', '#box', '#buffer']) $(id).addEventListener('change', () => { showModelInfo(); validateSoon(0); });
$('#build').onclick = build;
$('#reset').onclick = resetAll;
$('#batchtext').addEventListener('input', () => batchValidateSoon());
$('#batchexample').onclick = () => { $('#batchtext').value = BATCH_EXAMPLE; batchValidateSoon(0); };
$('#batchclear').onclick = () => { $('#batchtext').value = ''; batchValidateSoon(0); };
$('#batchopen').onclick = () => $('#batchfile').click();
$('#batchfile').onchange = (ev) => openBatchFile(ev.target.files[0]);
$('#batchbuild').onclick = batchBuild;
$('#batchstop').onclick = () => { batch.stop = true; $('#batchmsg').textContent = 'Stopping after this duplex…'; };
$('#batchzip').onclick = () => batch.zip && download('dsxnabuilder_batch.zip', batch.zip, 'application/zip');
$('#batchtable').addEventListener('click', (ev) => {
  const b = ev.target.closest('[data-view]');
  if (b) viewRecord(+b.dataset.view);
});
$('#batchhelp').onclick = (ev) => {
  ev.preventDefault();
  selectTab('docs');
  $('#batch-format').scrollIntoView({ behavior: 'smooth', block: 'start' });
};
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
