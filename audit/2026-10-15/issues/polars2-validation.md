With current packages (Polars 2.0.0, NumPy 2.5.4, Python 3.12), the existing fast unit audit returns 989 passed and one failure in test_fit_rejects_ambiguous_keys_and_nonbinary_labels_before_fitting[label].

score/fit.py validates manifest labels with manifest['y'].is_in([0, 1]). For a Float64 label column, Polars 2 raises InvalidOperationError because its list literal is Int64. Previously Polars coerced these operand types. The existing test expects the library's explicit ValueError for invalid/nonbinary labels before fitting.

Fix the validation with explicit compatible types while preserving rejection of null, nonfinite, fractional and nonbinary labels, and preserve side-effect-free failure before fitting or writing. Do not silently cast fractional labels to integers. Check sibling label validators for the same operand mismatch. This is a historical fitting API compatibility issue; fitted final models were separately retired by the author.

Reproduction: run the existing test above against current packages. The complete fast suite was bounded and completed in 198.84 seconds; no fix or dependency downgrade was made during the audit.
