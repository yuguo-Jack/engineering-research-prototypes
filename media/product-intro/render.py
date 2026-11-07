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
OUTPUT = HERE / "video"


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
    video = OUTPUT / "intro.mp4"
    webm = OUTPUT / "intro.webm"
    captions = OUTPUT / "captions.zh.vtt"
    poster = OUTPUT / "poster.jpg"
    run([executable("ffmpeg"), "-v", "error", "-y", "-i", master,
         "-map", "0:v:0", "-map", "0:a:0", "-c:v", "libx264", "-preset", "slow", "-crf", "23",
         "-profile:v", "main", "-level:v", "4.0", "-pix_fmt", "yuv420p", "-g", "60", "-keyint_min", "30",
         "-threads", "4", "-c:a", "aac", "-b:a", "128k", "-ar", "48000", "-ac", "2",
         "-movflags", "+faststart", "-sn", video], cwd=ROOT)
    run([executable("ffmpeg"), "-v", "error", "-y", "-i", master,
         "-map", "0:v:0", "-map", "0:a:0", "-c:v", "libvpx-vp9", "-crf", "34", "-b:v", "0",
         "-deadline", "good", "-cpu-used", "4", "-row-mt", "1", "-tile-columns", "2", "-threads", "4",
         "-g", "60", "-pix_fmt", "yuv420p", "-c:a", "libopus", "-b:a", "96k", "-sn", webm], cwd=ROOT)
    run([executable("ffmpeg"), "-v", "error", "-y", "-ss", "3", "-i", video,
         "-frames:v", "1", "-vf", "scale=1280:720", "-q:v", "2", poster], cwd=ROOT)
    probe = json.loads(subprocess.check_output([executable("ffprobe"), "-v", "error", "-show_format",
                                               "-show_streams", "-of", "json", str(video)], text=True))
    visual = next(stream for stream in probe["streams"] if stream["codec_type"] == "video")
    duration = float(probe["format"]["duration"])
    expected = json.loads((HERE / "timing.json").read_text(encoding="utf-8"))["duration"]
    if (visual["width"], visual["height"], visual["r_frame_rate"]) != (1920, 1080, "30/1") or abs(duration - expected) > .1:
        raise SystemExit("Unexpected video dimensions, frame rate or duration")
    if max(video.stat().st_size, webm.stat().st_size) > 10 * 1024 ** 2:
        raise SystemExit("A video exceeds the repository's 10 MiB per-file budget")
    manifest = {"schema_version": 3, "duration_seconds": duration, "width": 1920, "height": 1080,
                "edit": "v2-30s-web", "bpm": 140,
                "fps": 30, "hyperframes_version": "0.9.19", "video_codec": visual["codec_name"],
                "audio": "AAC stereo 48 kHz; synthetic Chinese narration and procedurally synthesized score",
                "webm": "VP9 Profile 0 + Opus; H.264 Main + AAC fallback in MP4",
                "captions": "Dedicated burned-in caption rail; companion WebVTT transcript", "files": []}
    for path in (video, webm, poster, captions):
        data = path.read_bytes()
        manifest["files"].append({"path": path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
