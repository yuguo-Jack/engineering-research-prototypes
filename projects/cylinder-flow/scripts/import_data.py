"""Import authorized original initial inputs; do not extract arbitrary archive paths."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python scripts/import_data.py AUTHORIZED_SUBMISSION.zip")
    manifest = json.loads((ROOT / "data/manifest.json").read_text(encoding="utf-8"))
    archive = Path(sys.argv[1]).resolve(strict=True)
    if digest(archive) != manifest["archive_sha256"]:
        raise ValueError("Submission SHA256 differs; original input identity cannot be established")
    destination = (ROOT / "data/external").resolve()
    if not destination.is_relative_to(ROOT):
        raise ValueError("External data directory resolves outside this project")
    destination.mkdir(parents=True, exist_ok=True)
    wanted = {row["tar_member"]: row for row in manifest["files"]}
    completed = set()
    with ZipFile(archive) as z, z.open(manifest["zip_member"]) as compressed, tarfile.open(fileobj=compressed, mode="r|gz") as tar:
        for member in tar:
            row = wanted.get(member.name)
            if row is None:
                continue
            if not member.isfile() or member.size != row["bytes"] or Path(row["path"]).name != row["path"]:
                raise ValueError("Invalid expected input member")
            target = destination / row["path"]
            if target.exists():
                if not target.is_file() or target.is_symlink() or digest(target) != row["sha256"]:
                    raise FileExistsError("Existing input differs: " + str(target))
            else:
                partial = target.with_name(target.name + ".partial")
                value = hashlib.sha256()
                created_partial = False
                try:
                    with partial.open("xb") as output, tar.extractfile(member) as stream:
                        created_partial = True
                        for block in iter(lambda: stream.read(4*1024*1024), b""):
                            output.write(block)
                            value.update(block)
                    if value.hexdigest() != row["sha256"]:
                        raise ValueError("Input checksum differs: " + row["path"])
                    partial.replace(target)
                except Exception:
                    if created_partial and partial.is_file():
                        partial.unlink()
                    raise
            completed.add(member.name)
    if completed != set(wanted):
        raise ValueError("Missing original initial inputs")
    print(json.dumps({"status": "verified", "files": len(completed), "bytes": sum(row["bytes"] for row in wanted.values()), "destination": str(destination)}))


if __name__ == "__main__":
    main()
