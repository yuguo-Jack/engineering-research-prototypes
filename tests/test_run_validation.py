"""Validate run-file checking against archived fixtures and deliberate corruptions."""
from pathlib import Path
import importlib.util
import json
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def validator(project_id):
    project = ROOT / "projects" / project_id
    spec = importlib.util.spec_from_file_location(project_id.replace("-", "_") + "_validator", project / "scripts/validate_run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, project


class RunValidation(unittest.TestCase):
    def setUp(self):
        (ROOT / ".local").mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix="run-check-", dir=ROOT / ".local")
        self.output = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def fixture(self, project_id):
        module, project = validator(project_id)
        source = project / "results/raw"
        if project_id == "cylinder-flow":
            source /= "formal-01"
        for name in ("summary.json", "history.csv"):
            shutil.copyfile(source / name, self.output / name)
        return module, project

    def alter_summary(self, key, value):
        path = self.output / "summary.json"
        summary = json.loads(path.read_text(encoding="utf-8"))
        summary[key] = value
        path.write_text(json.dumps(summary), encoding="utf-8")

    def test_archived_cylinder_fixture_passes_as_numerical_data(self):
        module, project = self.fixture("cylinder-flow")
        result = module.validate_directory(project, self.output)
        self.assertEqual(result["status"], "passed")
        self.assertFalse(result["checks_speed"])

    def test_nan_cylinder_residual_is_rejected(self):
        module, project = self.fixture("cylinder-flow")
        self.alter_summary("residual_scaled_max", float("nan"))
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            module.validate_directory(project, self.output)

    def test_changed_drag_is_rejected(self):
        module, project = self.fixture("cylinder-flow")
        self.alter_summary("cd", .5)
        with self.assertRaisesRegex(ValueError, "Numerical regression"):
            module.validate_directory(project, self.output)

    def test_failed_process_record_is_rejected(self):
        module, project = self.fixture("cylinder-flow")
        (self.output / "run-record.json").write_text('{"returncode": 1}', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "process failed"):
            module.validate_directory(project, self.output)

    def test_archived_m6_fixture_recomputes_force_window(self):
        module, project = self.fixture("m6-wing")
        result = module.validate_directory(project, self.output)
        self.assertEqual((result["window_start"], result["window_confirmation"]), (207, 306))
        self.assertFalse(result["official_acceptance_claimed"])

    def test_restarted_m6_is_rejected(self):
        module, project = self.fixture("m6-wing")
        self.alter_summary("restart_iteration_offset", 400)
        with self.assertRaisesRegex(ValueError, "Restarted"):
            module.validate_directory(project, self.output)

    def test_m6_truncated_history_is_rejected(self):
        module, project = self.fixture("m6-wing")
        path = self.output / "history.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            module.validate_directory(project, self.output)

    def test_repeated_probe_preserves_earlier_output_and_log(self):
        script = ROOT / "projects/m6-wing/reproduce.py"
        spec = importlib.util.spec_from_file_location("m6_repeat_test_driver", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ROOT = self.output.resolve()
        module.BUILD = module.ROOT / "build"
        previous = module.BUILD / "skip-postflow-assembly-probe"
        previous.mkdir(parents=True)
        (previous / "evidence.json").write_text('{"earlier": true}', encoding="utf-8")
        (module.BUILD / "ctest.log").write_text("earlier diagnostic", encoding="utf-8")
        module.preserve_probe_output()
        retained = list((module.BUILD / "test-history").iterdir())
        self.assertEqual(len(retained), 1)
        self.assertFalse(previous.exists())
        self.assertEqual((retained[0] / "evidence.json").read_text(encoding="utf-8"), '{"earlier": true}')
        self.assertEqual((retained[0] / "previous-ctest.log").read_text(encoding="utf-8"), "earlier diagnostic")

    def test_probe_output_outside_project_is_preserved_and_rejected(self):
        script = ROOT / "projects/m6-wing/reproduce.py"
        spec = importlib.util.spec_from_file_location("m6_outside_test_driver", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ROOT = self.output / "project"
        module.BUILD = self.output / "outside-build"
        previous = module.BUILD / "skip-postflow-assembly-probe"
        previous.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "escapes"):
            module.preserve_probe_output()
        self.assertTrue(previous.exists())


if __name__ == "__main__":
    unittest.main()
