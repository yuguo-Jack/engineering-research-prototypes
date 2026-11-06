"""Validate the public candidate, source/result identities, metadata, links and generated overview."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import argparse
import hashlib
import json
import re
import subprocess
import sys

from catalog import ROOT, generated_files, load_projects, local_file


def candidate_files(root=ROOT):
    result = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=root, check=True, capture_output=True)
    return sorted({Path(name.decode("utf-8")) for name in result.stdout.split(b"\0") if name})


def linked_file(root, source, target):
    if not target or target.startswith(("#", "mailto:", "data:")):
        return
    url = urlsplit(target)
    if url.scheme or url.netloc:
        return
    target_path = unquote(url.path)
    if not target_path:
        return
    path = (root / target_path.lstrip("/")) if target_path.startswith("/") else source.parent / target_path
    path = path.resolve()
    if not path.is_relative_to(root.resolve()) or not path.exists():
        raise ValueError(f"Broken local link: {source.relative_to(root)} -> {target}")


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.targets = []

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in ("src", "href") and value:
                self.targets.append(value)


def check_schema_fields(project):
    for field in ("id", "title", "summary", "domain", "status", "version", "license", "readme", "reproduce", "results"):
        if not isinstance(project[field], str) or not project[field].strip():
            raise ValueError("Invalid metadata string: " + field)
    for field in ("tags", "language", "commands", "limits"):
        if not isinstance(project.get(field, []), list) or any(not isinstance(x, str) for x in project.get(field, [])):
            raise ValueError("Invalid metadata list: " + field)
    if not isinstance(project.get("release_review_required", True), bool):
        raise ValueError("Invalid release review status")


def check_manifest(project_root, name, records):
    for record in records:
        path = local_file(project_root, record["path"])
        payload = path.read_bytes()
        if len(payload) != record["bytes"] or hashlib.sha256(payload).hexdigest() != record["sha256"]:
            raise ValueError(f"{name}: recorded file changed: {record['path']}")


def check(root=ROOT, release=False):
    _, projects = load_projects(root)
    warnings = []
    for project in projects:
        check_schema_fields(project)
        directory = root / "projects" / project["id"]
        for doc in ("LICENSE.md", "DATA.md"):
            local_file(directory, doc)
        source_manifest = directory / "SOURCE_MANIFEST.json"
        if source_manifest.is_file():
            source = json.loads(source_manifest.read_text(encoding="utf-8"))
            check_manifest(directory, project["id"], [{"path": path, **record} for path, record in source["files"].items()])
        result_manifest = directory / "results/PROVENANCE.json"
        if result_manifest.is_file():
            check_manifest(directory, project["id"], json.loads(result_manifest.read_text(encoding="utf-8"))["files"])
        if project.get("release_review_required", True) or project["license"].startswith("LicenseRef-Author-Review"):
            warnings.append(project["id"] + ": author license/release review pending")
    for name, expected in generated_files(root).items():
        file = root / name
        if not file.is_file() or file.read_text(encoding="utf-8") != expected:
            raise ValueError("Generated overview is stale: " + name)
    files = candidate_files(root)
    total = 0
    text_extensions = {".md", ".json", ".py", ".cpp", ".hpp", ".h", ".txt", ".log", ".yml", ".yaml", ".html", ".css", ".js", ".csv", ".sh"}
    for relative in files:
        path = root / relative
        if not path.is_file():
            raise ValueError("Staged deletion or missing candidate: " + str(relative))
        if path.is_symlink() or getattr(path, "is_junction", lambda: False)():
            raise ValueError("Do not publish filesystem links: " + str(relative))
        if any(part in (".local", "build", "runs", "external", "__pycache__") for part in relative.parts):
            raise ValueError("Local-only file in candidate: " + str(relative))
        if path.suffix.lower() in {".exe", ".dll", ".so", ".o", ".obj", ".a", ".lib", ".zip", ".gz", ".pdb"}:
            raise ValueError("Binary/archive in candidate: " + str(relative))
        if path.stat().st_size > 10*1024*1024:
            raise ValueError("Use external artifact storage for files over 10 MiB: " + str(relative))
        total += path.stat().st_size
        if path.suffix in text_extensions or path.name in ("LICENSE", "CMakeLists.txt"):
            text = path.read_text(encoding="utf-8-sig")
            if path.suffix == ".py":
                compile(text, str(path), "exec")
            if path.suffix == ".json":
                json.loads(text)
            private = (r"/public/home/", r"scnethr\d", r"-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----", r"gh[pousr]_[A-Za-z0-9]{30,}")
            # Do not interpret the checker source's literal detection rules as findings.
            if relative.as_posix() != "tools/check_repo.py" and any(re.search(pattern, text) for pattern in private):
                raise ValueError("Private/local material requires review: " + str(relative))
            if path.suffix == ".md":
                for target in re.findall(r"\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)", text):
                    linked_file(root, path, target)
            if path.suffix == ".html":
                links = Links()
                links.feed(text)
                for target in links.targets:
                    linked_file(root, path, target)
    if release and warnings:
        raise ValueError("Release review incomplete: " + "; ".join(warnings))
    return {"status": "passed", "projects": len(projects), "candidate_files": len(files), "candidate_bytes": total, "release_review": warnings, "scope": "repository_and_recorded_artifacts_not_a_new_physical_solver_run"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", action="store_true")
    args = parser.parse_args()
    try:
        result = check(release=args.release)
    except (ValueError, KeyError, FileNotFoundError, UnicodeError) as error:
        raise SystemExit(str(error)) from error
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
