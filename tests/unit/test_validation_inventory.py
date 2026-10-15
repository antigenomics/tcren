"""Reference spans must be unique at the least-truncated matching level."""

import pytest

from tcren.validation.inventory import compatible_sequence_span


@pytest.mark.parametrize(
    ("sequence", "query", "budget", "expected"),
    [
        ("CDGCDGAC", "ACD", 1, None),  # unique prefix cannot hide repeated suffix
        ("ACGACGCD", "ACD", 1, None),  # repeated prefix cannot hide unique suffix
        ("ACGCD", "ACD", 1, None),  # two distinct equally long matching spans
        ("AAA", "AAB", 1, None),  # overlapping occurrences also make a span ambiguous
        ("XXCDYY", "ACD", 1, {"reference_start": 2, "input_start": 1, "length": 2}),
        ("ACDCD", "ACD", 1, {"reference_start": 0, "input_start": 0, "length": 3}),
        ("XXEFYY", "ACD", 1, None),
    ],
)
def test_reference_span_preserves_all_ambiguity_at_first_matching_level(sequence, query, budget, expected):
    assert compatible_sequence_span(sequence, query, terminal_budget=budget) == expected
