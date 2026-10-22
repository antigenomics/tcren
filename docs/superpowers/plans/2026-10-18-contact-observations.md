---
orphan: true
---

# Contact observations implementation plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Export exact symmetry-aware native atom contacts and descriptive summaries through the existing batched CLI.

**Architecture:** Extend the existing neighbour search with opt-in atom indices, then join physical metadata in a dedicated contact-observation module. Route the new quantity through `export_batch`; unchanged quantities retain their current code path. This stage estimates no plausibility model and changes no energy.

**Tech Stack:** Existing Python, NumPy, Polars, SciPy cKDTree, Typer and pytest; no new dependency or native extension.

**Spec:** `CONTACT_PLAUSIBILITY_PLAN.md`, stage 1; author approved proceeding with the design and first implementation stage on 18 October 2026. Implementation is inline in this chat.

## Global constraints

- Inclusive 5 Å heavy-atom contacts, undirected chemical pair counts.
- Retain biological interface/region and physical author atom identities.
- No fabricated glycine Cβ, missing-geometry imputation, fitted score or changed legacy result.
- Whole-set annotation stays batched; default CPU budget one.
- Native/canonical copies are not independent observations.

## Review focus

- Alternate conformers with identical atom names must remain uniquely identifiable.
- Endpoint swaps must swap endpoint geometry and preserve unordered chemistry.
- Missing or degenerate orientation atoms must remain unavailable, never zero.
- Empty contact sets must retain typed schema and successful input accounting.
- MHC helix annotations must not include groove floor/B2M contacts by accident.

### Task 1: Exact observations and CLI integration

**Files:** `src/tcren/contacts/geometry.py`, new `src/tcren/contacts/observations.py`, `src/tcren/exports.py`, `src/tcren/cli.py`, `tests/unit/test_contact_observations.py`, `tests/unit/test_exports.py`.

**Interfaces:** `all_atom_contacts(..., atom_indices=False)` optionally adds flattened heavy-atom endpoint indices without changing defaults. `contact_observations(structure, *, cutoff=5.0, interface="all") -> pl.DataFrame` selects the three interface types, retains physical metadata and unordered chemistry, and adds observed Cβ orientation geometry. `summarize_contact_observations(table) -> pl.DataFrame` exports separate atom and residue frequencies; residue counts deduplicate physical edges.

- [x] Write a synthetic deterministic check for two alternate atom identities, inclusive cutoff, glycine/no Cβ, endpoint symmetry and biological scopes; verify the missing API fails.
- [x] Implement the minimal opt-in indices and observation layer; select one deterministic per-residue alternate conformer for this new quantity and record the selection policy. Keep old contacts untouched.
- [x] Add `export --quantity contact-atoms --interface all` with optional individual interface selection, no potential requirement, exact accounting and summary/hash lineage. Add empty/all-failed and unchanged-default checks.
- [x] Run focused contact/export tests with a 120 s process timeout, then lint and zero-warning docs. Update API skill/docs with actual schema and limitations.

### Task 2: Authorized native summary

**Files:** New dated compact source/count/verification record under `audit/2026-10-18/`; bulk observation deliverables under the established data root, never committed into source.

- [x] Verify the complete eligible raw inventory and frozen exclusion recipe, then run the new batched CLI once into a new output location. This is a new authorized summary, not an incidental rerun of previous analyses.
- [x] Independently check raw-coordinate distances, physical/reversed identities, group counts and all accounting. Report actual eligibility, atom/residue counts and missingness; retain failures explicitly.
- [x] Record source hashes, exact command, time/memory and downstream prior dependencies. Commit scoped source/docs/checks, preserving unrelated edits. Do not claim penalty or classifier validation from this stage.
