# Molecular publication CLI lineage — issue 34

Accepted on 16 October 2026 against library master 457f93f and its byte-matched
immutable 4.0.1 wheel (182 payload files; SHA256
`d24f599861773b9d3ff69c174c6059e14aa245b46e64878e57954611197dceba`).
This is not a release or a new ranking benchmark.

## Acceptance

The owning private benchmark records six sequential, fresh CLI exports for the
complete 13-structure set: position TP/PM, residue TP, contacts TP, placement and
pMHC. Every command accounts for 13 successes and zero errors. Original published
TCRen uses the canonical `karnaukhov2022` potential name; MHC uses `mj`, with a
5-Angstrom cutoff. The existing batch annotation API searches the complete query
set once per reference per command. No outer search or render pool is introduced.

All 15 active coordinate views consume those tables: four CPL views, four TCRvdb
views, four template-stratified VDJdb views and three native pMHC views. Class II
retains both groove chains, separately from class I. No core/register is guessed.
Figure 1 is native TikZ rather than a molecular render. Epitope/template groups
and all matched comparison rows remain unchanged.

The verifier checks physical chain/author number/insertion-code joins, annotated
regions, exact exported contact endpoints, actual atom colours/B values, potential
and searched-index identities, cameras, crop equality within pairs and enabled
scene atom coordinates. Rigid-transform residual is at most
2.11755213404731e-6 Angstrom, consistent with scene float precision. Hashes bind
CLI tables, raw inputs, producer sources, cameras, atom representations and scenes.
Two focused regression tests pass: insertion-aware exported identities/tampering
and high-resolution ray-buffer preservation. Existing exact-API equality checks
remain the library acceptance baseline; no scoring implementation was changed.

All 15 final PNGs were opened individually. The native composite, molecular SI
pages 11/15/16/17, inventory page 69 and final references page 70 were also visually
inspected. The displaced contact-empty receptor pose and unfavourable GLC
contributions remain visible. Exact numeric manifest comparisons are unchanged;
all pixel/crop differences are explicitly recorded. CPL/RAQ pixels match their
predecessors; FLR shares a 1574x1198 crop instead of 1577x1200. Complete original
manifests (including inactive examples and retired fitted fields) and original
image identities remain archived. Three native feature rows match all 14 recorded
columns exactly. No full cohort analysis was rerun.

The clean committed-source paper build has 19 main and 70 SI pages. Structure,
cross-reference and vocabulary gates pass; all 39 active PDF figures pass font
and bounds gates. No undefined reference/citation or overfull box remains.
Existing hyperref bookmark-level warnings are not new scientific failures.

## Production and resources

CLI TSVs feed native R accounting, seven group summaries and eight unchanged ipTM
macros; Python only formats the SI inventory. PyMOL draws authoritative values and
LaTeX composes figures. Normal paper builds consume recorded outputs. Explicit
fresh/render-only targets require new directories. A read-only molecular gate
checks the saved hashes without rewriting evidence.

The complete accepted fresh molecular block took 449.86591312501696 seconds.
Peak sampled simultaneous process-tree RSS was 1,763,311,616 bytes (1.64 GiB),
with 17 live OS threads including bounded auxiliary I/O. Actual MMseqs and PyMOL
processes each reached at most two threads, matching the selected two-CPU budget.
Sampling interval was 20 ms, with 16 transient process-read misses; sampled peaks
are not continuous proofs. The 1200-second/16-GiB bounds were not reached.
Vectorized/table consumption of existing batched APIs sufficed; no new native
code, dependency, fit or replacement score was needed.

## Owning records and remaining gates

Benchmark `MOLECULAR_CLI_LINEAGE.md` and
`results/molecular-panels-accepted-20261016/` own full exports, accounting, raw hashes,
resource/command receipts, comparisons and scenes. Manuscript `FIGURES.md`,
`MOLECULAR_CLI_LINEAGE.md`, `MOLECULAR_RENDER_QA.tsv` and the molecular-render
manifest own the include graph and final asset review. These remain in their
authorized private repositories; no coordinate/assay identifiers are copied here.
Unrelated dirty research outputs were preserved.

Issues 26, 32, 33 and 35 retain independent scope. Public release approval, the
long MD campaign and deterministic topology reductions remain separate gates.
