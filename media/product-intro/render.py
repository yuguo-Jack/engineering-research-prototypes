"""Render the film and package its small, seekable website assets."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = ROOT / ".local/video-production"
OUTPUT = ROOT / "assets/video"
NAME = "engineering-research-prototypes-intro"


def executable(name):
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"{name} is required on PATH; see this folder's README.md")
    return str(Path(found).resolve())


def run(command, **kwargs):
    subprocess.run([str(value) for value in command], check=True, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=HERE / "node_modules/hyperframes/bin/hyperframes.mjs")
    parser.add_argument("--browser", type=Path, help="Optional existing Chromium/Edge executable")
    parser.add_argument("--master", type=Path, default=WORK / "renders/master.mp4")
    parser.add_argument("--package-only", action="store_true", help="Use an already rendered master")
    args = parser.parse_args()
    master = args.master.resolve()
    if not master.is_relative_to(ROOT):
        raise SystemExit("Keep task masters and temporary files inside this repository")
    master.parent.mkdir(parents=True, exist_ok=True)
    scratch = WORK / "tmp"
    scratch.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update(TEMP=str(scratch), TMP=str(scratch), TMPDIR=str(scratch),
               HYPERFRAMES_NO_UPDATE_CHECK="1", HYPERFRAMES_NO_TELEMETRY="1",
               HYPERFRAMES_NO_AUTO_INSTALL="1")
    for name in ("ffmpeg", "ffprobe"):
        env["HYPERFRAMES_" + name.upper() + "_PATH"] = executable(name)
    if args.browser:
        env["HYPERFRAMES_BROWSER_PATH"] = str(args.browser.resolve())
    if not args.package_only:
        if not args.cli.is_file():
            raise SystemExit("Run npm ci in media/product-intro, or supply --cli for an existing shared installation")
        run([sys.executable, "-B", str(HERE / "build_composition.py")], cwd=ROOT)
        run([executable("node"), args.cli.resolve(), "render", "--fps", "30", "--workers", "2", "--crf", "21",
             "--frames-cache-dir", "off", "-o", master], cwd=HERE, env=env)
    if not master.is_file():
        raise SystemExit("The rendered master is missing")
    OUTPUT.mkdir(exist_ok=True)
    video = OUTPUT / f"{NAME}.mp4"
    captions = OUTPUT / f"{NAME}.zh.vtt"
    poster = OUTPUT / f"{NAME}-poster.jpg"
    run([executable("ffmpeg"), "-v", "error", "-y", "-i", master, "-i", captions,
         "-map", "0:v:0", "-map", "0:a:0", "-map", "1:0", "-c", "copy", "-c:s", "mov_text",
         "-metadata:s:s:0", "language=zho", "-metadata:s:s:0", "title=中文字幕",
         "-disposition:s:0", "0", "-metadata", "title=Engineering Research Prototypes",
         "-movflags", "+faststart", video], cwd=ROOT)
    run([executable("ffmpeg"), "-v", "error", "-y", "-ss", "3", "-i", video,
         "-frames:v", "1", "-vf", "scale=1280:720", "-q:v", "2", poster], cwd=ROOT)
    probe = json.loads(subprocess.check_output([executable("ffprobe"), "-v", "error", "-show_format",
                                               "-show_streams", "-of", "json", str(video)], text=True))
    visual = next(stream for stream in probe["streams"] if stream["codec_type"] == "video")
    duration = float(probe["format"]["duration"])
    expected = json.loads((HERE / "timing.json").read_text(encoding="utf-8"))["duration"]
    if (visual["width"], visual["height"], visual["r_frame_rate"]) != (1920, 1080, "30/1") or abs(duration - expected) > .1:
        raise SystemExit("Unexpected video dimensions, frame rate or duration")
    if video.stat().st_size > 10 * 1024 ** 2:
        raise SystemExit("Video exceeds the repository's 10 MiB per-file budget; revise encoding before publication")
    manifest = {"schema_version": 2, "duration_seconds": duration, "width": 1920, "height": 1080,
                "edit": "v2-30s", "bpm": 140,
                "fps": 30, "hyperframes_version": "0.9.19", "video_codec": visual["codec_name"],
                "audio": "AAC stereo 48 kHz; synthetic Chinese narration and procedurally synthesized score",
                "captions": "Dedicated burned-in caption rail; optional WebVTT and mov_text tracks, off by default", "files": []}
    for path in (video, poster, captions):
        data = path.read_bytes()
        manifest["files"].append({"path": path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
