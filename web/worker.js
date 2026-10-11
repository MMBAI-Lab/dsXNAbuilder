// dsXNAbuilder compute worker: self-hosted Pyodide + numpy + the dsxnabuilder package.
// Messages in:  {id, name, payload}  name in models | validate | complement | build
// Messages out: {id, ok, result | error}  and  {status: "..."} while loading.
import { loadPyodide } from './pyodide/pyodide.mjs';

const base = new URL('./', import.meta.url).href;

async function boot() {
  postMessage({ status: 'loading_python' });
  const py = await loadPyodide({ indexURL: base + 'pyodide/' });
  postMessage({ status: 'loading_numpy' });
  await py.loadPackage('numpy');
  postMessage({ status: 'loading_app' });
  const version = (await (await fetch(base + 'VERSION', { cache: 'no-cache' })).text()).trim();
  const zip = await (await fetch(base + 'dsxnabuilder.zip?v=' + encodeURIComponent(version))).arrayBuffer();
  py.unpackArchive(zip, 'zip', { extractDir: '/home/pyodide' });
  py.runPython('import sys; sys.path.insert(0, "/home/pyodide")\nfrom dsxnabuilder.webapi import call, session\nsession()');
  postMessage({ status: 'ready', version });
  return py.globals.get('call');
}

const ready = boot().catch((e) => { postMessage({ status: 'error', error: String(e) }); throw e; });

onmessage = async (ev) => {
  const { id, name, payload } = ev.data;
  try {
    const call = await ready;
    const out = JSON.parse(call(name, JSON.stringify(payload)));
    postMessage({ id, ...out });
  } catch (e) {
    postMessage({ id, ok: false, error: String(e) });
  }
};
