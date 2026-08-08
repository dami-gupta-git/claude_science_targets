import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def pytest_configure(config):
    """Fail with an actionable message when the interpreter lacks gemmi.

    The tests build synthetic structures with gemmi, so an interpreter without
    it errors during collection with a bare ModuleNotFoundError that names the
    test file rather than the cause. gemmi is not pip-installable into every
    environment, so the interpreter matters: run the suite with the conda env
    that has it.
    """
    try:
        import gemmi  # noqa: F401
    except ImportError:
        import pytest
        raise pytest.UsageError(
            "gemmi is not importable by %s.\n"
            "These tests need it. Either run them with an environment that has "
            "gemmi, e.g.\n"
            "    /Users/dgupta/.claude-science/conda/envs/structbio/bin/python "
            "-m pytest tests/ -q\n"
            "or install it first:\n"
            "    conda install -c conda-forge gemmi   # or: pip install gemmi"
            % sys.executable)
