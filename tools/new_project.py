"""Create an independent draft from the repository template; refuse overwrites."""
from pathlib import Path
from datetime import datetime, timezone, timedelta
import argparse
import json
import re
import shutil

from catalog import markdown_text

ROOT = Path(__file__).resolve().parents[1]


def create_project(root, slug, title, domain):
    if not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", slug):
        raise ValueError("Use a lowercase-kebab project id")
    if not title.strip() or any(character in title for character in "\r\n\0"):
        raise ValueError("A nonempty single-line project title is required")
    config = json.loads((root / "catalog-config.json").read_text(encoding="utf-8"))
    if domain not in config["domains"]:
        raise ValueError("Unknown domain")
    destination = root / "projects" / slug
    if destination.exists():
        raise FileExistsError("Project already exists; nothing was overwritten")
    shutil.copytree(root / "_templates/project", destination)
    date = datetime.now(timezone(timedelta(hours=8))).date().isoformat()
    for path in destination.rglob("*.md"):
        text = path.read_text(encoding="utf-8").replace("__TITLE__", markdown_text(title)).replace("__ID__", slug).replace("__DATE__", date)
        path.write_text(text, encoding="utf-8")
    metadata_path = destination / "project.json"
    project = json.loads(metadata_path.read_text(encoding="utf-8"))
    project.update(id=slug, title=title, domain=domain, updated=date)
    metadata_path.write_text(json.dumps(project, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("id")
    parser.add_argument("--title", required=True)
    parser.add_argument("--domain", default="other")
    args = parser.parse_args()
    try:
        result = create_project(ROOT, args.id, args.title, args.domain)
    except (ValueError, FileExistsError) as error:
        raise SystemExit(str(error)) from error
    print(json.dumps({"created": str(result), "status": "draft", "next": "Complete project files, then run python tools/catalog.py"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
