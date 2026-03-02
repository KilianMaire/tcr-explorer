"""
Tests to validate CI/CD workflow YAML files are well-formed
and contain the expected configuration.
"""

import pathlib
import yaml
import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
WORKFLOWS_DIR = REPO_ROOT / ".github" / "workflows"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_workflow(name: str) -> dict:
    """Load and parse a workflow YAML file."""
    path = WORKFLOWS_DIR / name
    assert path.exists(), f"Workflow file not found: {path}"
    with open(path) as f:
        data = yaml.safe_load(f)
    assert isinstance(data, dict), f"{name} did not parse to a dict"
    return data


def _flatten_steps(workflow: dict) -> list[dict]:
    """Return a flat list of all steps across all jobs."""
    steps = []
    for job in workflow.get("jobs", {}).values():
        steps.extend(job.get("steps", []))
    return steps


def _step_names(steps: list[dict]) -> list[str]:
    """Return list of step 'name' values (empty string for unnamed steps)."""
    return [s.get("name", "") for s in steps]


def _run_blocks(steps: list[dict]) -> str:
    """Concatenate all 'run' blocks into a single string for searching."""
    return "\n".join(s.get("run", "") for s in steps if "run" in s)


# ---------------------------------------------------------------------------
# General: files exist and are valid YAML
# ---------------------------------------------------------------------------

class TestWorkflowFilesExist:
    """Verify that required workflow files are present."""

    def test_test_yml_exists(self):
        assert (WORKFLOWS_DIR / "test.yml").is_file()

    def test_retrain_yml_exists(self):
        assert (WORKFLOWS_DIR / "retrain.yml").is_file()

    def test_workflows_dir_exists(self):
        assert WORKFLOWS_DIR.is_dir()


class TestWorkflowsAreValidYAML:
    """Verify that workflow files parse as valid YAML with required keys."""

    @pytest.mark.parametrize("filename", ["test.yml", "retrain.yml"])
    def test_yaml_parses(self, filename):
        data = _load_workflow(filename)
        assert "name" in data, f"{filename} missing 'name'"
        assert "on" in data or True in data, f"{filename} missing 'on' trigger"
        assert "jobs" in data, f"{filename} missing 'jobs'"


# ---------------------------------------------------------------------------
# test.yml specifics
# ---------------------------------------------------------------------------

class TestTestWorkflow:
    """Validate the test.yml CI workflow configuration."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.wf = _load_workflow("test.yml")
        self.steps = _flatten_steps(self.wf)
        self.names = _step_names(self.steps)
        self.runs = _run_blocks(self.steps)

    def test_name(self):
        assert self.wf["name"] == "Tests"

    def test_triggers_push_and_pr(self):
        triggers = self.wf.get("on") or self.wf.get(True, {})
        assert "push" in triggers
        assert "pull_request" in triggers

    def test_python_matrix(self):
        job = self.wf["jobs"]["test"]
        versions = job["strategy"]["matrix"]["python-version"]
        assert "3.11" in versions

    def test_installs_requirements_dev(self):
        assert "requirements-dev.txt" in self.runs

    def test_does_not_install_pytest_separately(self):
        # Should use requirements-dev.txt, not separate pip install pytest
        for step in self.steps:
            run = step.get("run", "")
            if "pip install" in run and "requirements" not in run:
                assert "pytest" not in run, (
                    "pytest should come from requirements-dev.txt, not a separate pip install"
                )

    def test_cache_key_includes_requirements_dev(self):
        for step in self.steps:
            cache_with = step.get("with", {})
            key = cache_with.get("key", "")
            if "pip" in key and "hashFiles" in key:
                assert "requirements-dev.txt" in key, (
                    "Cache key should hash requirements-dev.txt"
                )

    def test_ignores_broken_test_files(self):
        assert "test_cdr_enricher.py" in self.runs
        assert "test_nl_query.py" in self.runs
        assert "--ignore=" in self.runs

    def test_coverage_flags_present(self):
        assert "--cov" in self.runs

    def test_coverage_upload_step(self):
        assert any("coverage" in n.lower() for n in self.names), (
            "Expected a coverage upload/report step"
        )

    def test_pythonpath_set(self):
        assert "PYTHONPATH" in self.runs

    def test_env_pythondontwritebytecode(self):
        for step in self.steps:
            env = step.get("env", {})
            if "PYTHONDONTWRITEBYTECODE" in env:
                assert env["PYTHONDONTWRITEBYTECODE"] == "1"
                return
        pytest.fail("PYTHONDONTWRITEBYTECODE not set in any step")


# ---------------------------------------------------------------------------
# retrain.yml specifics
# ---------------------------------------------------------------------------

class TestRetrainWorkflow:
    """Validate the retrain.yml model retraining workflow."""

    @pytest.fixture(autouse=True)
    def _load(self):
        self.wf = _load_workflow("retrain.yml")
        self.steps = _flatten_steps(self.wf)
        self.names = _step_names(self.steps)
        self.runs = _run_blocks(self.steps)

    def test_name(self):
        assert self.wf["name"] == "Model Retraining"

    def test_workflow_dispatch_inputs(self):
        triggers = self.wf.get("on") or self.wf.get(True, {})
        dispatch = triggers.get("workflow_dispatch", {})
        inputs = dispatch.get("inputs", {})
        assert "max_epochs" in inputs
        assert "batch_size" in inputs

    def test_push_paths_reference_root_dvc_yaml(self):
        """dvc.yaml is at repo root, not data/dvc.yaml."""
        triggers = self.wf.get("on") or self.wf.get(True, {})
        paths = triggers.get("push", {}).get("paths", [])
        assert "dvc.yaml" in paths, "Should trigger on root dvc.yaml"
        assert "data/dvc.yaml" not in paths, (
            "data/dvc.yaml is outdated; dvc.yaml is at repo root"
        )

    def test_dvc_pull_step(self):
        assert any("dvc pull" in s.get("run", "") for s in self.steps)

    def test_prepare_training_data_correct_args(self):
        assert "--output data/processed/" in self.runs
        assert "--output data/processed/training_data.parquet" not in self.runs

    def test_prepare_step_has_pythonpath(self):
        for step in self.steps:
            run = step.get("run", "")
            if "prepare_training_data.py" in run:
                assert "PYTHONPATH" in run, (
                    "Prepare step must set PYTHONPATH"
                )
                return
        pytest.fail("prepare_training_data.py step not found")

    def test_generate_embeddings_step(self):
        assert "generate_esm2_embeddings.py" in self.runs

    def test_train_uses_script(self):
        assert "scripts/train_model.py" in self.runs
        # Should NOT use inline python -c
        for step in self.steps:
            run = step.get("run", "")
            if "train" in step.get("name", "").lower():
                assert 'python -c' not in run, (
                    "Train step should use scripts/train_model.py, not inline python"
                )

    def test_upload_artifacts_step(self):
        assert any("upload-artifact" in str(s.get("uses", "")) for s in self.steps)

    def test_aws_secrets_in_dvc_pull(self):
        for step in self.steps:
            if "dvc pull" in step.get("run", ""):
                env = step.get("env", {})
                assert "AWS_ACCESS_KEY_ID" in env
                assert "AWS_SECRET_ACCESS_KEY" in env
                return
        pytest.fail("DVC pull step with AWS secrets not found")


# .dvcignore removed (DVC migrated to imgt-ml repo)
