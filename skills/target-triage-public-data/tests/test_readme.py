"""Tests for the run-README generator.

All inputs are synthesised; no DepMap release or network access is needed.
"""

import numpy as np
import pandas as pd
import pytest

from kernel import (
    FINDING_MAX_WORDS,
    SUMMARY_MAX_WORDS,
    format_p,
    markdown_table,
    triage_readme,
    write_triage_readme,
)

SUMMARY = "A short plain-prose paragraph that a non-specialist can follow."
STEPS = [
    {"name": "Dependency", "finding": "Dispensable in most lines.",
     "table": [{"stratum": "MSI", "n": 75, "mean": -0.854, "p": 4e-12}]},
    {"name": "Sensitivity", "skipped": "No inhibitor is present in the screens on disk."},
]
FILES = [("out.csv", "the numbers")]
SOURCES = ["DepMap 24Q2"]
LIMITS = ["Knockout is not inhibition."]


def render(**overrides):
    kwargs = dict(gene="WRN", summary=SUMMARY, steps=STEPS, files=FILES,
                  data_sources=SOURCES, limits=LIMITS)
    kwargs.update(overrides)
    return triage_readme(**kwargs)


class TestFormatP:
    def test_readable_values_stay_decimal(self):
        assert format_p(0.21) == "0.21"
        assert format_p(0.001) == "0.001"

    def test_small_values_use_superscript_notation(self):
        assert format_p(4.0e-12) == "4.0 × 10⁻¹²"
        assert format_p(2.5e-4) == "2.5 × 10⁻⁴"

    def test_underflowed_zero_reports_a_bound_not_zero(self):
        # A p that underflowed to 0.0 must not print as "0", which claims more
        # precision than the float carries.
        assert format_p(0.0) == "< 1 × 10⁻³⁰⁰"

    def test_missing_and_nan_are_not_significant(self):
        assert format_p(None) == "n/a"
        assert format_p(float("nan")) == "n/a"


class TestMarkdownTable:
    def test_column_order_follows_headers_then_key_order(self):
        rows = [{"b": 2, "a": 1}]
        assert markdown_table(rows).splitlines()[0] == "| b | a |"
        assert markdown_table(rows, headers=["a", "b"]).splitlines()[0] == "| a | b |"

    def test_p_columns_are_formatted_as_p_values(self):
        out = markdown_table([{"p": 4e-12, "p_bh": 1e-3, "other_p": 0.5, "value": 4e-12}])
        assert "4.0 × 10⁻¹²" in out          # the p column
        assert "4e-12" in out                # the non-p float column, plain 3sf

    def test_q_columns_are_formatted_as_p_values(self):
        # sensitivity_scan's BH columns feed straight into these tables; a q of
        # 1e-30 printed as a plain float beside a p in scientific notation is
        # the same number rendered two ways in one row.
        out = markdown_table([{"q_BH": 4e-12, "q_partial_BH": 2.5e-4, "r": 4e-12}])
        assert "4.0 × 10⁻¹²" in out
        assert "2.5 × 10⁻⁴" in out
        assert "4e-12" in out                # the non-q float column, plain 3sf

    def test_columns_merely_containing_p_or_q_are_left_alone(self):
        # The match is on whole underscore-separated components, so a column
        # named for a quantity rather than a probability keeps float rendering.
        out = markdown_table([{"expr": 4e-12, "auc_top": 4e-12, "quantile": 4e-12}])
        assert "×" not in out

    def test_missing_and_nan_render_as_dash(self):
        out = markdown_table([{"a": None, "b": np.nan}])
        assert out.splitlines()[-1] == "| — | — |"

    def test_accepts_a_dataframe(self):
        assert "| n |" in markdown_table(pd.DataFrame([{"n": 3}]))

    def test_empty_rows_raise_rather_than_render_a_headerless_table(self):
        with pytest.raises(ValueError, match="no rows"):
            markdown_table([])


class TestStructure:
    def test_sections_appear_in_the_standard_order(self):
        text = render()
        positions = [text.index(h) for h in ("# WRN target triage", SUMMARY,
                                             "## Result", "## Files",
                                             "## Data sources", "## Limits")]
        assert positions == sorted(positions)

    def test_summary_precedes_any_heading_or_number(self):
        # The opening paragraph is what a non-specialist reads; nothing may
        # come before it but the title.
        body = render().split("\n", 1)[1].strip()
        assert body.startswith(SUMMARY)

    def test_step_tables_render_under_their_step(self):
        text = render()
        assert text.index("### Dependency") < text.index("| stratum |") < text.index("### Sensitivity")

    def test_skipped_step_is_reported_not_omitted(self):
        text = render()
        assert "### Sensitivity" in text
        assert "Not run. No inhibitor" in text

    def test_title_can_be_overridden(self):
        assert render(title="WRN in MSI tumours").startswith("# WRN in MSI tumours")


class TestGuards:
    def test_summary_over_cap_names_the_count_and_the_cap(self):
        with pytest.raises(ValueError, match=f"over the {SUMMARY_MAX_WORDS}-word cap"):
            render(summary="word " * (SUMMARY_MAX_WORDS + 1))

    def test_summary_exactly_at_cap_is_accepted(self):
        # The boundary belongs to the passing side.
        assert render(summary="word " * SUMMARY_MAX_WORDS)

    def test_finding_over_cap_names_its_step(self):
        long_step = [{"name": "Dependency", "finding": "word " * (FINDING_MAX_WORDS + 1)}]
        with pytest.raises(ValueError, match="finding for step 'Dependency'"):
            render(steps=long_step)

    def test_bulleted_summary_is_rejected(self):
        with pytest.raises(ValueError, match="plain prose"):
            render(summary="- a bullet\n- another")

    def test_step_without_finding_or_skipped_is_rejected(self):
        with pytest.raises(ValueError, match="neither 'finding' nor 'skipped'"):
            render(steps=[{"name": "Dependency"}])

    def test_step_without_a_name_is_rejected(self):
        with pytest.raises(ValueError, match="no 'name'"):
            render(steps=[{"finding": "something"}])

    def test_no_steps_is_rejected(self):
        with pytest.raises(ValueError, match="no steps"):
            render(steps=[])


def test_write_triage_readme_writes_what_it_returns(tmp_path):
    path = tmp_path / "README.md"
    text = write_triage_readme(str(path), gene="WRN", summary=SUMMARY, steps=STEPS,
                               files=FILES, data_sources=SOURCES, limits=LIMITS)
    assert path.read_text(encoding="utf-8") == text
    assert text.endswith("\n")


def test_write_triage_readme_is_utf8_regardless_of_locale(tmp_path, monkeypatch):
    # format_p emits '×' and superscript digits for any p below 0.001, which
    # every real run produces. An unpinned open() takes the platform default
    # encoding and raises UnicodeEncodeError on a non-UTF-8 locale, after the
    # analysis has already completed.
    #
    # The default is substituted at builtins.open rather than through the locale
    # module: CPython reads the locale encoding below the Python level, so
    # patching locale.getpreferredencoding leaves open() on UTF-8 and the test
    # passes with the fix reverted.
    import builtins
    real_open = builtins.open

    def ascii_default_open(file, mode="r", *args, encoding=None, **kwargs):
        if encoding is None and "b" not in mode:
            encoding = "ascii"
        return real_open(file, mode, *args, encoding=encoding, **kwargs)

    monkeypatch.setattr(builtins, "open", ascii_default_open)
    path = tmp_path / "README.md"
    text = write_triage_readme(str(path), gene="WRN", summary=SUMMARY, steps=STEPS,
                               files=FILES, data_sources=SOURCES, limits=LIMITS)
    monkeypatch.undo()
    assert "×" in text and "⁻" in text       # the characters that force the issue
    assert path.read_bytes().decode("utf-8") == text
