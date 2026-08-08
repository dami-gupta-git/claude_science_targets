import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def pytest_configure(config):
    """Fail with an actionable message when the interpreter lacks statsmodels.

    sensitivity_scan() needs statsmodels for BH correction, which is not on
    every interpreter. Run with the conda env that has it.
    """
    try:
        import statsmodels  # noqa: F401
        import scipy  # noqa: F401
    except ImportError as e:
        import pytest
        raise pytest.UsageError(
            "%s is not importable by %s.\n"
            "These tests need numpy, pandas, scipy and statsmodels. Run them with "
            "an interpreter that has all four, or install them first:\n"
            "    pip install scipy statsmodels"
            % (e.name, sys.executable))
