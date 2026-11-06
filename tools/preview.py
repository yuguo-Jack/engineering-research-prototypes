"""Read-only loopback preview restricted to Git's public candidate file list."""
from pathlib import Path
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import unquote, urlsplit
import argparse
import os
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def public_files(root=ROOT):
    result = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=root, check=True, capture_output=True)
    return frozenset(name.decode("utf-8") for name in result.stdout.split(b"\0") if name)


def preview_path(url, root, allowed):
    relative = unquote(urlsplit(url).path).lstrip("/")
    if not relative or relative.endswith("/"):
        relative += "index.html"
    target = (root / relative).resolve()
    if not target.is_relative_to(root.resolve()) or not target.is_file():
        return None
    if target.relative_to(root.resolve()).as_posix() not in allowed:
        return None
    if target.is_symlink() or getattr(target, "is_junction", lambda: False)():
        return None
    return target


class PreviewHandler(SimpleHTTPRequestHandler):
    allowed_files = frozenset()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def send_head(self):
        self._byte_range = None
        target = preview_path(self.path, ROOT, self.allowed_files)
        if target is None:
            self.send_error(404)
            return None
        requested = self.headers.get("Range")
        if requested:
            size = target.stat().st_size
            match = re.fullmatch(r"bytes=(\d*)-(\d*)", requested.strip())
            start, end = 0, size - 1
            if match and any(match.groups()):
                first, last = match.groups()
                if first:
                    start = int(first)
                    end = min(int(last), end) if last else end
                elif int(last):
                    start = max(0, size - int(last))
                else:
                    start = size
            else:
                start = size
            if start > end or start >= size:
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{size}")
                self.send_header("Content-Length", "0")
                self.end_headers()
                return None
            stream = target.open("rb")
            stream.seek(start)
            self._byte_range = (start, end)
            self.send_response(206)
            self.send_header("Content-Type", self.guess_type(str(target)))
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            self.send_header("Content-Length", str(end - start + 1))
            self.send_header("Last-Modified", self.date_time_string(os.fstat(stream.fileno()).st_mtime))
            self.end_headers()
            return stream
        return super().send_head()

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        super().end_headers()

    def copyfile(self, source, outputfile):
        try:
            if self._byte_range is None:
                return super().copyfile(source, outputfile)
            remaining = self._byte_range[1] - self._byte_range[0] + 1
            while remaining:
                block = source.read(min(65536, remaining))
                if not block:
                    break
                outputfile.write(block)
                remaining -= len(block)
        except (BrokenPipeError, ConnectionResetError):
            # Browsers cancel outstanding media requests when seeking.
            pass

    def translate_path(self, path):
        # Use exactly the same normalized target for validation and actual reads.
        target = preview_path(path, ROOT, self.allowed_files)
        return str(target) if target is not None else str(ROOT / ".preview-denied")

    def list_directory(self, path):
        self.send_error(404)
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    PreviewHandler.allowed_files = public_files()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), PreviewHandler)
    print(f"Overview: http://127.0.0.1:{args.port}/   (Ctrl+C to stop)", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
