"""Initialize package-owned native pools before test fixtures import numerical libraries.

Already-initialized external Polars pools are exercised in an isolated subprocess.
"""
import tcren  # noqa: F401
