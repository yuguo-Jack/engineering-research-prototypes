"""Regression checks for generic drafts, literal metadata and readable static docs."""
from pathlib import Path
from html.parser import HTMLParser
import json
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import catalog
import check_repo
import new_project
import preview


class CatalogReview(unittest.TestCase):
    def setUp(self):
        (ROOT / ".local").mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix="catalog-review-", dir=ROOT / ".local")
        self.root = Path(self.temporary.name)
        for folder in ("_templates", "schemas", "docs"):
            shutil.copytree(ROOT / folder, self.root / folder)
        for name in ("CONTRIBUTING.md", "LICENSE"):
            shutil.copyfile(ROOT / name, self.root / name)
        shutil.copyfile(ROOT / "catalog-config.json", self.root / "catalog-config.json")
        (self.root / "README.md").write_text("# Projects\n<!-- PROJECTS:START -->\n<!-- PROJECTS:END -->\n", encoding="utf-8")
        (self.root / "site").mkdir()
        for name in ("styles.css", "app.js"):
            shutil.copyfile(ROOT / "site" / name, self.root / "site" / name)

    def tearDown(self):
        self.temporary.cleanup()

    def project(self, title="Radar prototype"):
        return new_project.create_project(self.root, "radar-demo", title, "signals")

    def write_generated(self):
        files = catalog.generated_files(self.root)
        for name, text in files.items():
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
        return files

    def test_draft_has_no_invented_cfd_evidence_or_empty_copy_button(self):
        self.project()
        files = self.write_generated()
        detail = files["site/projects/radar-demo.html"]
        self.assertNotIn("SOURCE_MANIFEST.json", detail)
        self.assertNotIn("results/VALIDATION.json", detail)
        self.assertNotIn("data-copy=", detail)
        self.assertIn("尚未提供可运行命令", detail)

    def test_all_links_in_a_new_draft_resolve(self):
        self.project()
        files = self.write_generated()
        for name, text in files.items():
            if name.endswith(".html"):
                parser = check_repo.Links()
                parser.feed(text)
                for target in parser.targets:
                    check_repo.linked_file(self.root, self.root / name, target)

    def test_readme_and_overview_use_the_same_new_catalog_in_one_pass(self):
        self.project()
        files = catalog.generated_files(self.root)
        self.assertIn("radar-demo", files["site/docs/README.html"])
        self.assertIn('href="site/docs/projects/radar-demo/REPRODUCE.html"', files["index.html"])

    def test_title_markup_stays_literal_in_the_readme_and_document_page(self):
        self.project('<img src=x onerror="bad()"> [prototype] | one')
        files = catalog.generated_files(self.root)
        self.assertNotIn('<img src="x"', files["site/docs/projects/radar-demo/README.html"])
        self.assertIn("&lt;img", files["site/docs/projects/radar-demo/README.html"])
        self.assertIn(r"\| one", files["README.md"])

    def test_nested_metadata_error_is_actionable(self):
        project = self.project()
        path = project / "project.json"
        metadata = json.loads(path.read_text(encoding="utf-8"))
        metadata["highlights"] = [{"label": "Rate", "value": {"invalid": True}}]
        path.write_text(json.dumps(metadata), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, r"highlights\.0\.value"):
            catalog.generated_files(self.root)

    def test_search_data_contains_display_domain_and_language(self):
        project = self.project()
        path = project / "project.json"
        metadata = json.loads(path.read_text(encoding="utf-8"))
        metadata["language"] = ["MATLAB"]
        path.write_text(json.dumps(metadata), encoding="utf-8")
        overview = catalog.generated_files(self.root)["index.html"]
        self.assertIn("MATLAB signals 通信与雷达", overview)

    def test_project_title_must_be_single_line(self):
        with self.assertRaises(ValueError):
            self.project("Radar\nInjected heading")

    def test_removed_project_exposes_stale_generated_detail(self):
        project = self.project()
        self.write_generated()
        (project / "project.json").unlink()
        obsolete = catalog.obsolete_generated(self.root, catalog.generated_files(self.root))
        self.assertIn(self.root / "site/projects/radar-demo.html", obsolete)

    def test_link_into_an_ignored_local_file_is_rejected(self):
        private = self.root / ".local/private.md"
        private.parent.mkdir()
        private.write_text("private content", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unpublished/ignored"):
            check_repo.linked_file(self.root, self.root / "index.html", ".local/private.md", {"index.html"})

    def test_preview_allows_only_explicit_public_files(self):
        (self.root / "index.html").write_text("public page", encoding="utf-8")
        private = self.root / ".local/private.txt"
        private.parent.mkdir()
        private.write_text("private work", encoding="utf-8")
        allowed = frozenset({"index.html"})
        self.assertEqual(preview.preview_path("/", self.root, allowed), self.root / "index.html")
        self.assertIsNone(preview.preview_path("/.local/private.txt", self.root, allowed))
        self.assertIsNone(preview.preview_path("/%2e%2e/private.txt", self.root, allowed))
        self.assertIsNone(preview.preview_path("/site/", self.root, allowed))


if __name__ == "__main__":
    unittest.main()
