import http from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.dirname(fileURLToPath(import.meta.url));
const port = Number(process.env.PORT || 4173);
const types = { '.html': 'text/html; charset=utf-8', '.css': 'text/css; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.json': 'application/json; charset=utf-8', '.svg': 'image/svg+xml', '.png': 'image/png', '.pdb': 'chemical/x-pdb', '.mol2': 'chemical/x-mol2', '.sdf': 'chemical/x-mdl-sdfile' };
const server = http.createServer(async (request, response) => {
  if (!['GET', 'HEAD'].includes(request.method)) { response.writeHead(405, { Allow: 'GET, HEAD' }); response.end(); return; }
  try {
    const pathname = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
    const relative = pathname === '/' ? 'index.html' : pathname.replace(/^\/+/, '');
    const allowed = ['index.html', 'styles.css'].includes(relative) || /^(src|assets|data)\/[a-zA-Z0-9_./-]+$/.test(relative);
    const location = path.resolve(root, relative);
    if (!allowed || !location.startsWith(root + path.sep) || relative.split('/').some(part => part.startsWith('.'))) { response.writeHead(404); response.end('Not found'); return; }
    if (!(await stat(location)).isFile()) throw new Error('Not a file');
    const content = await readFile(location);
    response.writeHead(200, { 'Content-Type': types[path.extname(location)] || 'application/octet-stream', 'Cache-Control': 'no-cache', 'X-Content-Type-Options': 'nosniff', 'Content-Length': content.length });
    response.end(request.method === 'HEAD' ? undefined : content);
  } catch { response.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' }); response.end('Not found'); }
});
server.on('error', error => { console.error(`Unable to start: ${error.message}`); process.exitCode = 1; });
server.listen(port, '127.0.0.1', () => console.log(`Polymer Programming\nLocal: http://127.0.0.1:${port}\nPress Ctrl+C to stop.`));
