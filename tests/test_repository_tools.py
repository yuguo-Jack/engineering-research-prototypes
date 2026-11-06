"""Behavioral checks for independent project creation, safe imports and catalog paths."""
from pathlib import Path
from unittest.mock import patch
import hashlib
import importlib.util
import io
import json
import shutil
import sys
import tarfile
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import catalog
import new_project


def load_importer():
    spec = importlib.util.spec_from_file_location("project_importer_under_test", ROOT / "projects/cylinder-flow/scripts/import_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RepositoryTools(unittest.TestCase):
    def setUp(self):
        (ROOT / ".local").mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix="unit-", dir=ROOT / ".local")
        self.root = Path(self.temporary.name)
        shutil.copytree(ROOT / "_templates", self.root / "_templates")
        shutil.copyfile(ROOT / "catalog-config.json", self.root / "catalog-config.json")
        (self.root / "README.md").write_text("# Example\n<!-- PROJECTS:START -->\n<!-- PROJECTS:END -->\n", encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def test_new_project_is_a_valid_independent_draft(self):
        project = new_project.create_project(self.root, "radar-demo", 'Radar "demo" <prototype>', "signals")
        _, records = catalog.load_projects(self.root)
        self.assertEqual(records[0]["id"], "radar-demo")
        self.assertEqual(records[0]["status"], "draft")
        self.assertTrue((project / "REPRODUCE.md").is_file())
        generated = catalog.generated_files(self.root)
        self.assertIn("radar-demo", generated["index.html"])
        self.assertNotIn('<prototype>', generated["index.html"])
        self.assertIn("&lt;prototype&gt;", generated["index.html"])

    def test_creation_refuses_overwrite(self):
        project = new_project.create_project(self.root, "my-demo", "First", "other")
        content = (project / "README.md").read_bytes()
        with self.assertRaises(FileExistsError):
            new_project.create_project(self.root, "my-demo", "Second", "other")
        self.assertEqual(content, (project / "README.md").read_bytes())

    def test_creation_rejects_path_traversal(self):
        with self.assertRaises(ValueError):
            new_project.create_project(self.root, "../outside", "Bad", "other")
        self.assertFalse((self.root.parent / "outside").exists())

    def test_catalog_rejects_escaping_document(self):
        project = new_project.create_project(self.root, "my-demo", "Demo", "other")
        (self.root / "private.md").write_text("Not a project file", encoding="utf-8")
        with self.assertRaises(ValueError):
            catalog.local_file(project, "../../private.md")

    def test_generated_output_is_deterministic(self):
        new_project.create_project(self.root, "my-demo", "Demo", "controls")
        self.assertEqual(catalog.generated_files(self.root), catalog.generated_files(self.root))

    def importer_fixture(self):
        project = self.root / "import-project"
        (project / "data").mkdir(parents=True)
        payload = b"toy original initial input\n"
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
            member = tarfile.TarInfo("inputs/toy.bin")
            member.size = len(payload)
            tar.addfile(member, io.BytesIO(payload))
        archive = self.root / "toy.zip"
        with zipfile.ZipFile(archive, "w") as z:
            z.writestr("inner.tar.gz", buffer.getvalue())
        manifest = {"archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(), "zip_member": "inner.tar.gz", "files": [{"path": "toy.bin", "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(), "tar_member": "inputs/toy.bin"}]}
        (project / "data/manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        module = load_importer()
        module.ROOT = project.resolve()
        return module, project, archive, payload

    def test_import_is_verified_and_idempotent(self):
        module, project, archive, payload = self.importer_fixture()
        with patch.object(sys, "argv", ["import_data.py", str(archive)]):
            module.main()
            module.main()
        self.assertEqual(payload, (project / "data/external/toy.bin").read_bytes())

    def test_import_preserves_preexisting_partial(self):
        module, project, archive, _ = self.importer_fixture()
        (project / "data/external").mkdir()
        partial = project / "data/external/toy.bin.partial"
        partial.write_bytes(b"other unfinished work")
        with patch.object(sys, "argv", ["import_data.py", str(archive)]):
            with self.assertRaises(FileExistsError):
                module.main()
        self.assertEqual(b"other unfinished work", partial.read_bytes())

    def test_import_refuses_different_existing_input(self):
        module, project, archive, _ = self.importer_fixture()
        (project / "data/external").mkdir()
        target = project / "data/external/toy.bin"
        target.write_bytes(b"important existing input")
        with patch.object(sys, "argv", ["import_data.py", str(archive)]):
            with self.assertRaises(FileExistsError):
                module.main()
        self.assertEqual(b"important existing input", target.read_bytes())


if __name__ == "__main__":
    unittest.main()
