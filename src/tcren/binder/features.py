"""Native extraction of the 5 binder descriptors from a TCR-pMHC structure.

All compute is native/fast: interface size, chain balance, H-bonds and ΔSASA come from the
``tcren._geom`` C kernel; the CDR1/2-vs-CDR3α potential term is a TCRen sum over the
``ContactMap`` interface. No PyRosetta, no Biopython SASA, no sklearn. The structure is chain-typed
(and TCR V(D)J-annotated for the CDR regions) via :func:`tcren.annotation.classify_chains`.

``binder_features(structure) -> dict`` feeds :func:`tcren.binder.binder_score`.
"""

from __future__ import annotations

import numpy as np

from ..structure.model import PEPTIDE_TYPE, Structure

_TCR = ("TRA", "TRB")
_MHC = ("MHC", "MHCa", "MHCb", "B2M")
POLAR = {"N", "O"}
BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}


def _role_atoms(structure, roles, *, polar_only=False, with_radii=False):
    """Heavy-atom coords (and per-atom residue index / radii) over chains whose type is in ``roles``."""
    xyz, res, rad = [], [], []
    for c in structure.chains:
        if c.chain_type not in roles:
            continue
        for r in c.residues:
            for a in r.atoms:
                if a.element == "H" or (polar_only and a.element not in POLAR):
                    continue
                xyz.append(a.coord)
                res.append(r.seq_index)
                rad.append(BONDI.get(a.element, 1.70))
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    if with_radii:
        return xyz, np.asarray(rad, float)
    return xyz, np.asarray(res, np.int32)


def _dsasa(structure):
    """Interface ΔSASA = (SASA(TCR alone)+SASA(pMHC alone)) − (their SASA in the full complex)."""
    from .. import _geom

    xyz, rad, tag = [], [], []  # 0 = TCR, 1 = pMHC (peptide+MHC)
    for c in structure.chains:
        role = 0 if c.chain_type in _TCR else (1 if c.chain_type in (PEPTIDE_TYPE,) + _MHC else -1)
        if role < 0:
            continue
        for r in c.residues:
            for a in r.atoms:
                if a.element == "H":
                    continue
                xyz.append(a.coord)
                rad.append(BONDI.get(a.element, 1.70))
                tag.append(role)
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    rad = np.asarray(rad, float)
    tag = np.asarray(tag)
    if len(xyz) == 0 or (tag == 0).sum() == 0 or (tag == 1).sum() == 0:
        return float("nan")
    full = _geom.shrake_rupley(xyz, rad, 1.4, 100)
    tcr_alone = _geom.shrake_rupley(xyz[tag == 0], rad[tag == 0], 1.4, 100).sum()
    pmhc_alone = _geom.shrake_rupley(xyz[tag == 1], rad[tag == 1], 1.4, 100).sum()
    return float((tcr_alone + pmhc_alone) - (full[tag == 0].sum() + full[tag == 1].sum()))


def _pp_combo(structure, potential, cutoff):
    """CDR1/2-vs-CDR3α TCRen potential term: z(ΣJ over CDR1/2:pep) − z(ΣJ over CDR3α:pep)."""
    import polars as pl

    from ..contactmap import ContactMap
    from .model import BINDER_MODEL

    mat, idx = potential.as_matrix()

    def sum_e(df):
        if df.height == 0:
            return 0.0
        t = 0.0
        for a, b in zip(df["residue.aa.from"].to_list(), df["residue.aa.to"].to_list()):
            ia, ib = idx.get(a, -1), idx.get(b, -1)
            if ia >= 0 and ib >= 0:
                t += mat[ia, ib]
        return t

    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    df = cm.interface("tcr_peptide", tcr_regions="all")
    reg, ch = df["region.type.from"], df["chain.type.from"]
    e12 = sum_e(df.filter(reg.is_in(["CDR1", "CDR2"])))
    e3a = sum_e(df.filter((reg == "CDR3") & (ch == "TRA")))
    (m12, s12), (m3a, s3a) = BINDER_MODEL["pp_z"]["cdr12"], BINDER_MODEL["pp_z"]["cdr3a"]
    return (e12 - m12) / s12 - (e3a - m3a) / s3a


def binder_features(structure, *, potential=None, organism: str = "human",
                    cutoff: float = 5.0) -> dict[str, float]:
    """Extract the 5 native binder descriptors from a TCR-pMHC structure (path or parsed).

    The structure is chain-typed + TCR-annotated in place. Returns a dict keyed by
    :data:`tcren.binder.FEATURES`, ready for :func:`tcren.binder.binder_score`.
    """
    from .. import _geom
    from ..annotation import classify_chains
    from ..potential import tcren as _tcren_pot
    from ..structure import parse_structure

    if not isinstance(structure, Structure):
        structure = parse_structure(structure)
    if all(c.chain_type is None for c in structure.chains):
        classify_chains(structure, organism=organism)

    tra_xyz, tra_res = _role_atoms(structure, ("TRA",))
    trb_xyz, trb_res = _role_atoms(structure, ("TRB",))
    pep_xyz, _pr = _role_atoms(structure, (PEPTIDE_TYPE,))
    mhc_xyz, _mr = _role_atoms(structure, _MHC)
    tcr_pol, _ = _role_atoms(structure, _TCR, polar_only=True)
    pep_pol, _ = _role_atoms(structure, (PEPTIDE_TYPE,), polar_only=True)

    cd = _geom.contact_descriptors(tra_xyz, tra_res, trb_xyz, trb_res, pep_xyz, mhc_xyz, 5.0, 4.5)
    n_hbond = _geom.interface_hbonds(tcr_pol, pep_pol, 3.5) if len(tcr_pol) and len(pep_pol) else 0
    pot = potential or _tcren_pot()
    return {
        "pm_cov_ntcr": float(cd["pm_cov_ntcr"]),
        "chain_balance": float(cd["chain_balance"]),
        "n_hbond": float(n_hbond),
        "dSASA": _dsasa(structure),
        "pp_combo": _pp_combo(structure, pot, cutoff),
    }
