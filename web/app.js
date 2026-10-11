// dsXNAbuilder front end. All computation runs in worker.js (Pyodide); this file is UI only.
import { UI, ISSUES } from './i18n.js?v=0.2.0-dev';

const $ = (s) => document.querySelector(s);
const $$ = (s) => [...document.querySelectorAll(s)];
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

const EXAMPLES = [
  { id: 'dickerson', en: 'Dickerson dodecamer (B-DNA)', es: 'Dodecámero de Dickerson (B-DNA)', top: 'CGCGAATTCGCG', bottom: 'GCGCTTAAGCGC', model: 4 },
  { id: 'arna', en: 'RNA duplex (A-RNA)', es: 'Dúplex de ARN (A-RNA)', top: 'GGCGCUUGCGUC', bottom: 'CCGCGAACGCAG', model: 20 },
  { id: 'hybrid', en: 'RNA:DNA hybrid', es: 'Híbrido ARN:ADN', top: 'r(GACUAGCUGA)', bottom: 'd(CTGATCGACT)', model: 22 },
  { id: 'adna', en: 'A-DNA', es: 'ADN-A', top: 'GGTATACC', bottom: 'CCATATGG', model: 54 },
  { id: 'zdna', en: 'Z-DNA (CG)₃', es: 'ADN-Z (CG)₃', top: 'CGCGCG', bottom: 'GCGCGC', model: 15 },
  { id: 'mismatch', en: 'Mismatch + overhang', es: 'Desapareamiento + extremo', top: 'CGCGAGTTCGCGAA', bottom: 'GCGCTAAAGCGC--', model: 55 },
  { id: 'sirna', en: 'siRNA 2′-OMe/2′-F + PS', es: 'siRNA 2′-OMe/2′-F + PS', top: 'mG*mA*fC*mU*fA*mG*fC*mU*mG*mA', bottom: 'C*U*GAUCGAC*U', model: 20 },
  { id: 'lna', en: 'LNA gapmer : RNA', es: 'Gapmer LNA : ARN', top: 'r(GACUAGCUGA)', bottom: '+C*+T*+G*dA*dT*dC*dG*+A*+C*+T', model: 22 },
  { id: 'moe', en: 'MOE gapmer : RNA', es: 'Gapmer MOE : ARN', top: 'r(GACUAGCUGA)', bottom: 'eC*eT*eG*dA*dT*dC*dG*eA*eC*eT', model: 22 },
  { id: 'mc', en: '5mC + PS on B-DNA', es: '5mC + PS en ADN-B', top: 'CGC^5mCGAA*2TTCGCG', bottom: 'GCGCTT*1AAGCGC', model: 4 },
];

const state = {
  lang: 'en', ready: false, models: [], byId: new Map(), last: null, stage: null, comp: null,
  repr: 'cartoon', spinning: false, validateTimer: null, validSeq: 0,
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
  } else if (d.status === 'error') {
    $('#status').innerHTML = `<span class="error">${esc(t('engine_error'))}${esc(d.error)}</span>`;
  } else {
    $('#status').innerHTML = `<span class="spinner"></span>${esc(t(d.status))}`;
  }
}

// ------------------------------------------------------------------ i18n
function t(key, params = {}) {
  const s = (UI[state.lang] && UI[state.lang][key]) ?? UI.en[key] ?? key;
  return fill(s, params);
}
function fill(s, params) {
  return s.replace(/\{(\w+)\}/g, (m, k) => (params[k] !== undefined && params[k] !== null ? params[k] : m));
}
function issueText(i) {
  const tmpl = ISSUES[state.lang] && ISSUES[state.lang][i.code];
  if (!tmpl) return i.message;
  let s = fill(tmpl, i.params || {});
  const p = i.params || {};
  if (p.count > 1 && p.columns && !tmpl.includes('{columns}') && !tmpl.includes('{count}')) {
    s = (state.lang === 'es' ? `${p.count} residuos (columnas ${p.columns}): ` : `${p.count} residues (columns ${p.columns}): `) + s;
  }
  return s;
}
function setLang(lang) {
  state.lang = lang === 'es' ? 'es' : 'en';
  document.documentElement.lang = state.lang;
  try { localStorage.setItem('dsxna-lang', state.lang); } catch (e) {}
  for (const b of $$('#lang button')) b.setAttribute('aria-pressed', String(b.dataset.v === state.lang));
  for (const el of $$('[data-i18n]')) {
    const v = UI[state.lang][el.dataset.i18n];
    if (v) el.textContent = v;
  }
  for (const a of $$('#view-docs article')) a.classList.toggle('hidden', a.lang !== state.lang);
  renderExamples();
  if (state.models.length) { fillModelSelect(); renderModelTables(); showModelInfo(); }
  if (state.lastValidation) renderValidation(state.lastValidation);
  if (state.last) renderSummary(state.last);
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
  const pol = m.policy === 'any' ? t('m_generic') : m.policy === 'swap' ? t('m_swap') : t('m_pattern');
  return `#${m.id} · ${m.family} · ${m.twist.toFixed(1)}°/${m.rise.toFixed(2)} Å · ${m.native[0]} · ${pol}`;
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
    .map((m) => `<option value="${m.id}">${esc(modelLabel(m))}</option>`).join('')}</optgroup>`).join('');
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
  const pol = t('policy_' + m.policy, { native: m.native[0] });
  $('#modelinfo').textContent = t('info_line', {
    twist: m.twist.toFixed(1), rise: m.rise.toFixed(2), hand: t(m.handedness),
    s1: m.strands[0], s2: m.strands[1], p1: m.pucker[0] || '?', p2: m.pucker[1] || '?',
  }) + ' · ' + pol + ' · ' + m.description;
}
function renderModelTables() {
  const head = `<thead><tr><th>${t('m_id')}</th><th>${t('m_family')}</th><th>${t('m_twist')}</th><th>${t('m_rise')}</th><th>${t('m_hand')}</th><th>${t('m_repeat')}</th><th>${t('m_strands')}</th><th>${t('m_seq')}</th><th>${t('m_desc')}</th><th></th></tr></thead>`;
  const rows = state.models.map((m) => `<tr${m.id === +$('#model').value ? ' class="sel"' : ''}>
    <td class="num">${m.id}</td><td>${esc(m.family)}</td><td class="num">${m.twist.toFixed(2)}</td><td class="num">${m.rise.toFixed(3)}</td>
    <td>${t(m.handedness)}</td><td class="num">${m.k}</td><td>${m.strands.join('/')}</td>
    <td class="seq">${esc(m.native[0])}/${esc(m.native[1])}</td><td class="desc">${esc(m.description)}</td>
    <td><button class="btn ghost small" data-use="${m.id}">${t('use')}</button></td></tr>`).join('');
  $('#modeltable').innerHTML = head + `<tbody>${rows}</tbody>`;
  $('#excludedtable').innerHTML = `<thead><tr><th>${t('m_id')}</th><th>${t('m_family')}</th><th>${t('m_desc')}</th><th>${t('m_reason')}</th></tr></thead><tbody>` +
    (state.excluded || []).map((m) => `<tr><td class="num">${m.id}</td><td>${esc(m.family)}</td><td class="desc">${esc(m.description)}</td><td class="desc">${esc((m.issues || []).join('; '))}</td></tr>`).join('') + '</tbody>';
  for (const b of $$('[data-use]')) b.onclick = () => { $('#model').value = b.dataset.use; showModelInfo(); selectTab('build'); validateSoon(0); };
}

// ------------------------------------------------------------------ request
function request() {
  return {
    top: $('#top').value, bottom: $('#bottom').value, model: +$('#model').value || 4,
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
    $('#grid').innerHTML = `<span class="hint">${esc(t('grid_empty'))}</span>`;
    $('#issues').innerHTML = '';
    $('#build').disabled = true;
    state.lastValidation = null;
    return;
  }
  const mySeq = ++state.validSeq;
  try {
    const v = await call('validate', req);
    if (mySeq !== state.validSeq) return;
    state.lastValidation = v;
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
      <span class="tok${c.top ? (c.top_mod ? ' mod' : '') : ' gap'}">${esc(c.top || '–')}</span>
      <span class="pair ${c.pair || ''}">${c.pair ? sym[c.pair] : ''}</span>
      <span class="tok${c.bottom ? (c.bottom_mod ? ' mod' : '') : ' gap'}">${esc(c.bottom || '–')}</span>
      <span class="num">${(c.i + 1) % 5 === 0 || c.i === 0 ? c.i + 1 : ''}</span></div>`).join('');
    $('#grid').className = 'grid';
    $('#grid').innerHTML = `<div class="end"><div>5′</div><div></div><div>3′</div></div>${cols}<div class="end"><div>3′</div><div></div><div>5′</div></div>`;
  }
  renderIssues(v.issues || []);
}
function renderIssues(list) {
  const order = { error: 0, warning: 1, info: 2 };
  $('#issues').innerHTML = [...list].sort((a, b) => order[a.level] - order[b.level]).map((i) =>
    `<li class="${i.level}">${esc(issueText(i))}<span class="code">${esc(i.code)}</span></li>`).join('');
}

// ------------------------------------------------------------------ build
async function build() {
  if (!state.ready) return;
  const req = request();
  $('#build').disabled = true;
  $('#status').innerHTML = `<span class="spinner"></span>${esc(t('building'))}`;
  try {
    const v = await call('validate', req);
    renderValidation(v);
    if (!v.ok) return;
    req.libs = await fetchLibs(v.needed_libs || []);
    if ((v.needed_libs || []).length) req.frcmod = await fetchText('lib/frcmod.modxna');
    const out = await call('build', req);
    renderIssues(out.issues || []);
    if (!out.ok) return;
    state.last = out;
    renderSummary(out);
    state.lastReq = req;
    try {
      await showStructure(out.pdb, out.summary);
    } catch (e) {
      $('#viewer').innerHTML = `<div class="empty">${esc(t('no_webgl'))}</div>`;
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
  $('#summary').innerHTML = stat(s.n_bp, t('s_bp')) + stat(s.n_atoms, t('s_atoms')) +
    stat(`#${s.model}`, `${t('s_model')} · ${s.family}`) +
    stat(s.o3p ? `${s.o3p[0].toFixed(2)}–${s.o3p[1].toFixed(2)}` : '–', t('s_o3p')) +
    stat(s.clashes, t('s_clashes')) + (s.n_modified ? stat(s.n_modified, t('s_mod')) : '');
}

// ------------------------------------------------------------------ viewer
function bg() { return document.documentElement.getAttribute('data-theme') === 'light' ? '#f2f2f2' : '#0b0b0b'; }
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
    state.stage = new NGL.Stage('viewer', { backgroundColor: bg(), tooltip: true });
    window.addEventListener('resize', () => state.stage.handleResize());
  }
  state.stage.removeAllComponents();
  state.comp = await state.stage.loadFile(new Blob([pdb], { type: 'text/plain' }), { ext: 'pdb', name: 'dsxnabuilder' });
  state.modSele = Object.values(summary.name_map || {}).map((n) => `[${n}]`).join(' or ');
  applyRepr();
  state.comp.autoView();
}
function applyRepr() {
  const c = state.comp;
  if (!c) return;
  c.removeAllRepresentations();
  if (state.repr === 'cartoon') {
    c.addRepresentation('cartoon', { colorScheme: 'chainid', aspectRatio: 4, radiusScale: 1.2 });
    c.addRepresentation('base', { colorScheme: 'resname' });
  } else {
    c.addRepresentation('licorice', { colorScheme: 'element' });
  }
  if (state.modSele) c.addRepresentation('ball+stick', { sele: state.modSele, colorScheme: 'element' });
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
  history.replaceState(null, '', '?' + p.toString());
}
function applyURL() {
  const p = new URLSearchParams(location.search);
  if (p.get('lang')) setLang(p.get('lang'));
  if (p.get('tab')) selectTab(p.get('tab'));
  const set = (id, k) => { if (p.get(k) !== null) $(id).value = p.get(k); };
  set('#top', 'top'); set('#bottom', 'bottom'); set('#model', 'model'); set('#dnaff', 'ff');
  set('#water', 'water'); set('#box', 'box'); set('#buffer', 'buffer');
  showModelInfo();
  if (p.get('top') && p.get('bottom')) build(); else validateSoon(0);
}

// ------------------------------------------------------------------ tabs, examples
function selectTab(name) {
  for (const tab of ['build', 'models', 'docs']) {
    $(`#tab-${tab}`).setAttribute('aria-selected', String(tab === name));
    $(`#view-${tab}`).classList.toggle('hidden', tab !== name);
  }
  if (name === 'build' && state.stage) setTimeout(() => state.stage.handleResize(), 0);
}
function renderExamples() {
  $('#examples').innerHTML = EXAMPLES.map((e) => `<button class="btn ghost small" type="button" data-ex="${e.id}">${esc(e[state.lang] || e.en)}</button>`).join('');
  for (const b of $$('[data-ex]')) b.onclick = () => {
    const e = EXAMPLES.find((x) => x.id === b.dataset.ex);
    $('#top').value = e.top; $('#bottom').value = e.bottom; $('#model').value = e.model;
    showModelInfo(); build();
  };
}

// ------------------------------------------------------------------ wiring
for (const tab of ['build', 'models', 'docs']) $(`#tab-${tab}`).onclick = () => selectTab(tab);
for (const b of $$('#lang button')) b.onclick = () => setLang(b.dataset.v);
$('#theme').onclick = () => {
  const cur = document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
  document.documentElement.setAttribute('data-theme', cur);
  try { localStorage.setItem('theme', cur); } catch (e) {}
  if (state.stage) state.stage.setParameters({ backgroundColor: bg() });
};
for (const id of ['#top', '#bottom']) $(id).addEventListener('input', () => validateSoon());
for (const id of ['#model', '#water', '#dnaff', '#box', '#buffer']) $(id).addEventListener('change', () => { showModelInfo(); validateSoon(0); });
$('#build').onclick = build;
$('#complement').onclick = async () => {
  try {
    const r = await call('complement', { top: $('#top').value, model: +$('#model').value });
    $('#bottom').value = r.bottom;
    validateSoon(0);
  } catch (e) { $('#issues').innerHTML = `<li class="error">${esc(e.message)}</li>`; }
};
$('#clear').onclick = () => { $('#top').value = ''; $('#bottom').value = ''; validateSoon(0); };
for (const b of $$('#repr button')) b.onclick = () => {
  state.repr = b.dataset.v;
  for (const x of $$('#repr button')) x.setAttribute('aria-pressed', String(x === b));
  applyRepr();
};
$('#spin').onclick = () => { if (state.stage) { state.spinning = !state.spinning; state.stage.setSpin(state.spinning); } };
$('#recenter').onclick = () => { if (state.comp) state.comp.autoView(); };
$('#dlzip').onclick = async () => {
  // the Amber/modXNA package is an advanced option: built on demand with the
  // current options of the Advanced panel, for the last built duplex
  if (!state.last || !state.lastReq) return;
  const req = { ...state.lastReq, amber: request().amber, package: true };
  $('#status').innerHTML = `<span class="spinner"></span>${esc(t('building'))}`;
  try {
    const v = await call('validate', req);
    req.libs = await fetchLibs(v.needed_libs || []);
    if ((v.needed_libs || []).length) req.frcmod = await fetchText('lib/frcmod.modxna');
    const out = await call('build', req);
    if (out.ok) download('dsxnabuilder-amber.zip', b64bytes(out.zip_b64), 'application/zip');
    else renderIssues(out.issues || []);
  } catch (e) {
    $('#issues').innerHTML = `<li class="error">${esc(e.message)}</li>`;
  } finally { $('#status').textContent = ''; }
};
$('#dlpdb').onclick = () => state.last && download('dsxnabuilder.pdb', state.last.pdb, 'chemical/x-pdb');
$('#link').onclick = async () => {
  try { await navigator.clipboard.writeText(location.href); $('#status').textContent = t('link_copied'); } catch (e) { prompt('URL', location.href); }
};

let lang = 'en';
try { lang = localStorage.getItem('dsxna-lang') || ((navigator.language || '').startsWith('es') ? 'es' : 'en'); } catch (e) {}
setLang(lang);
