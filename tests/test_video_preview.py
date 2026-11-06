"""The local website must support media seeking without exposing private files."""
from http.server import ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen
import sys
import tempfile
import threading
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import preview


class VideoPreview(unittest.TestCase):
    def test_media_ranges_and_private_path_boundary(self):
        class QuietPreview(preview.PreviewHandler):
            allowed_files = frozenset({"sample.mp4"})

            def log_message(self, *args):
                pass

        with tempfile.TemporaryDirectory(dir=ROOT / ".local", prefix="video-http-") as work:
            folder = Path(work)
            (folder / "sample.mp4").write_bytes(b"0123456789")
            (folder / "private.txt").write_text("not public", encoding="utf-8")
            with patch.object(preview, "ROOT", folder):
                server = ThreadingHTTPServer(("127.0.0.1", 0), QuietPreview)
                thread = threading.Thread(target=server.serve_forever, daemon=True)
                thread.start()
                base = f"http://127.0.0.1:{server.server_port}"
                try:
                    for span, expected, content_range in (("bytes=3-7", b"34567", "bytes 3-7/10"),
                                                          ("bytes=-4", b"6789", "bytes 6-9/10"),
                                                          ("bytes=8-999", b"89", "bytes 8-9/10")):
                        with urlopen(Request(base + "/sample.mp4", headers={"Range": span}), timeout=5) as response:
                            self.assertEqual(response.status, 206)
                            self.assertEqual(response.headers["Content-Range"], content_range)
                            self.assertEqual(response.read(), expected)
                    with urlopen(Request(base + "/sample.mp4", headers={"Range": "bytes=3-7"}, method="HEAD"), timeout=5) as response:
                        self.assertEqual(response.headers["Content-Length"], "5")
                        self.assertEqual(response.read(), b"")
                    for path, headers, status in (("/sample.mp4", {"Range": "bytes=10-"}, 416),
                                                  ("/private.txt", {"Range": "bytes=0-"}, 404)):
                        with self.assertRaises(HTTPError) as caught:
                            urlopen(Request(base + path, headers=headers), timeout=5)
                        self.assertEqual(caught.exception.code, status)
                finally:
                    server.shutdown()
                    server.server_close()
                    thread.join(timeout=5)


if __name__ == "__main__":
    unittest.main()
