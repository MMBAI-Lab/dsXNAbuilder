"""Local static server for web/ with the MIME types Pyodide needs (.mjs, .wasm).
Usage: python3 scripts/serve.py [port]"""
import functools, http.server, mimetypes, os, sys

mimetypes.add_type('text/javascript', '.mjs')
mimetypes.add_type('application/wasm', '.wasm')
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'web')
port = int(sys.argv[1]) if len(sys.argv) > 1 else 8766
handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=root)
print(f'serving {os.path.abspath(root)} on http://localhost:{port}/')
http.server.ThreadingHTTPServer(('127.0.0.1', port), handler).serve_forever()
