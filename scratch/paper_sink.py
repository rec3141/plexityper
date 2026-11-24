#!/usr/bin/env python3
"""Local sink for papers fetched in Chrome through the library proxy (grid approved 248914999900553216).
POST /<filename> with the file body; saved to discussion/papers/. Binds 127.0.0.1 only."""
import http.server, os, re
OUT = '/data/dev/plexityper/discussion/papers'
class H(http.server.BaseHTTPRequestHandler):
    def cors(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', '*')
        self.send_header('Access-Control-Allow-Private-Network', 'true')
    def do_OPTIONS(self):
        self.send_response(204); self.cors(); self.end_headers()
    def do_POST(self):
        name = re.sub(r'[^A-Za-z0-9._-]', '_', os.path.basename(self.path)) or 'paper.bin'
        n = int(self.headers.get('Content-Length', 0)); data = self.rfile.read(n)
        with open(os.path.join(OUT, name), 'wb') as f: f.write(data)
        self.send_response(200); self.cors(); self.end_headers(); self.wfile.write(f'saved {name} {len(data)}'.encode())
http.server.HTTPServer(('127.0.0.1', 8765), H).serve_forever()
