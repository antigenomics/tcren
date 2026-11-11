"""Fixed experimental CDR3 chemical-neighborhood rules; no fitted parameters.

All heavy-atom witnesses are reduced once per physical residue edge. Peptide
and MHC edges remain separate. These diagnostics do not change screen defaults.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..clashes import BONDI_RADII
from ..contact_types import classify_contacts, stacked_pairs, _charged, _elem
from ..contacts.geometry import BACKBONE_ATOMS
from ..contactmap import ContactMap
from ..contacts.observations import select_conformer

SEVERE_OVERLAP_A = 0.6
CHARGE_DISTANCE_A = 4.0
CHEMISTRY = ('salt_bridge', 'hydrogen_bond', 'cation_pi', 'stacking', 'aromatic', 'hydrophobic')
RULES = {
    'gate_no_cdr3_peptide': 'score 0 when neither CDR3 has a physical peptide edge; otherwise prepared',
    'gate_both_cdr3_pmhc': 'score 0 unless each CDR3 has a physical peptide or MHC edge',
    'gate_chemical_peptide': 'score 0 without any nonconflicting CDR3-sidechain chemical peptide witness',
    'gate_both_chemical_pmhc': 'score 0 unless each CDR3 has a nonconflicting sidechain chemical pMHC witness',
    'pep_support': 'discount topology by fraction of peptide edges having chemical sidechain support',
    'pep_packing': 'discount topology by mean best supported packing per peptide edge',
    'pep_conflict': 'discount topology by one minus conflicted peptide-edge fraction',
    'mhc_conflict': 'discount topology by one minus conflicted MHC-edge fraction',
    'mhc_unanchored_conflict': 'discount topology by conflicted MHC edges at residues without chemical peptide support, divided by all MHC edges',
    'combined_conflict': 'discount topology by product of peptide and MHC conflict-free fractions',
}
EDGE = ['chain.id.from', 'residue.index.from', 'chain.id.to', 'residue.index.to']
SOURCE = ['chain.id.from', 'residue.index.from']


def neighborhood_counts(atoms: pl.DataFrame, *, stacking=None) -> dict:
    """Reduce directed CDR3 atom pairs; interface is peptide or mhc, d in (0,5].

    A chemical witness uses existing v2 types, a CDR3 sidechain atom, and no
    severe clash or same-sign charge at that atom pair. Edge conflict means ANY
    witness clashes or repels, even if another witness supports that edge.
    Packing min(d/r,r/d) is bounded by 1; r sums the two Bondi radii. An edge
    takes its best supported witness, so atom multiplicity supplies no weight.
    """
    required = set(EDGE+['chain.type.from', 'region.type.from', 'residue.aa.from',
                         'residue.aa.to', 'atom.from', 'atom.to', 'dist', 'interface'])
    if required-set(atoms.columns):
        raise ValueError('complete directed atom contact columns required')
    atoms = atoms.filter(pl.col('region.type.from') == 'CDR3')
    if set(atoms['interface'])-{'peptide','mhc'} or set(atoms['chain.type.from'])-{'TRA','TRB'}:
        raise ValueError('CDR3 atom contacts require TRA/TRB and peptide/mhc interfaces')
    distance = atoms['dist'].to_numpy()
    if not np.isfinite(distance).all() or ((distance <= 0) | (distance > 5)).any():
        raise ValueError('atom distance must be finite in (0,5] Angstrom')
    if atoms.select(EDGE+['atom.from','atom.to']).is_duplicated().any():
        raise ValueError('duplicate physical atom witness')
    typed = classify_contacts(atoms, scheme='v2', stacking=stacking)
    radius = np.zeros(atoms.height)
    charges = []
    for side in ('from','to'):
        # A small dictionary of unique atom types, rather than per-edge Python typing.
        keys = list(zip(atoms['residue.aa.'+side], atoms['atom.'+side]))
        lookup = {key: (_charged(*key), BONDI_RADII.get(_elem(key[1]), 1.70)) for key in set(keys)}
        values = np.array([lookup[key][0] for key in keys], dtype=bool).reshape(-1,2)
        charges.append(values)
        radius += np.array([lookup[key][1] for key in keys])
    same_charge = np.any(charges[0] & charges[1], axis=1) & (distance <= CHARGE_DISTANCE_A)
    bad = ((radius-distance) > SEVERE_OVERLAP_A) | same_charge
    chemical = typed.select(pl.any_horizontal(pl.col('is_'+name) for name in CHEMISTRY)).to_series().to_numpy()
    supported = chemical & ~bad & ~atoms['atom.from'].is_in(BACKBONE_ATOMS).to_numpy()
    packing = np.zeros(atoms.height)
    packing[supported] = np.minimum(distance[supported]/radius[supported], radius[supported]/distance[supported])
    edges = atoms.with_columns(pl.Series('supported',supported), pl.Series('conflicted',bad),
                              pl.Series('packing',packing)).group_by(EDGE+['interface','chain.type.from']).agg(
                                  pl.col('supported').any(), pl.col('conflicted').any(), pl.col('packing').max())
    anchored = edges.filter((pl.col('interface')=='peptide') & pl.col('supported')).select(SOURCE).unique()
    edges = edges.join(anchored.with_columns(anchored=pl.lit(True)), on=SOURCE, how='left', validate='m:1').with_columns(pl.col('anchored').fill_null(False)).sort(EDGE)
    result = dict(pmhc_chains=edges['chain.type.from'].n_unique(),
                  chemical_pmhc_chains=edges.filter(pl.col('supported'))['chain.type.from'].n_unique())
    for interface in ('peptide','mhc'):
        part = edges.filter(pl.col('interface')==interface)
        result.update({interface+'_edges':part.height,
                       interface+'_supported':int(part['supported'].sum()),
                       interface+'_conflicted':int(part['conflicted'].sum()),
                       interface+'_packing_sum':float(part['packing'].sum())})
    result['mhc_unanchored_conflicted'] = edges.filter((pl.col('interface')=='mhc') & pl.col('conflicted') & ~pl.col('anchored')).height
    return result


def _neighborhood_context(structure, *, cutoff=5.0):
    """Select one conformer for both physical contacts and stacking chemistry."""
    if cutoff != 5.0:
        raise ValueError('experimental neighborhood rules use the fixed 5 Angstrom cutoff')
    structure = select_conformer(structure)
    cm = ContactMap.from_structure(structure, cutoff=cutoff, atom_pairs=True)
    frames = [cm.interface(interface).filter(pl.col('region.type.from')=='CDR3')
              .with_columns(pl.lit(name).alias('interface'))
              for name,interface in [('peptide','tcr_peptide'),('mhc','tcr_mhc')]]
    return structure, pl.concat(frames)


def neighborhood_atoms(structure, *, cutoff=5.0) -> pl.DataFrame:
    """Fresh conformer-selected atoms on both complete peptide/MHC interfaces."""
    return _neighborhood_context(structure, cutoff=cutoff)[1]


def neighborhood_features(structure, *, cutoff=5.0) -> dict:
    structure, atoms = _neighborhood_context(structure, cutoff=cutoff)
    return {'pdb.id':structure.pdb_id, **neighborhood_counts(atoms, stacking=stacked_pairs(structure))}


def neighborhood_contacts(structure, *, cutoff=5.0) -> pl.DataFrame:
    """All CDR3 atom witnesses, complete peptide/MHC interfaces, v2 chemistry."""
    structure, atoms = _neighborhood_context(structure, cutoff=cutoff)
    return classify_contacts(atoms, scheme='v2', stacking=stacked_pairs(structure)).with_columns(
        (~pl.col('atom.from').is_in(BACKBONE_ATOMS)).alias('sidechain.from'))


def neighborhood_scores(models: pl.DataFrame) -> pl.DataFrame:
    """Fixed counterfactuals from prepared components; outcomes are never read.

    Soft rules keep w*ipTM unchanged and multiply only (1-w)*topology by Q.
    Empty conflict denominators imply Q=1; empty support implies Q=0. Hard
    gates set the complete score to 0. Models and finite score coverage persist.
    """
    required = ['pdb.id','iptm','topology','confidence_weight','pmhc_chains','chemical_pmhc_chains',
                'peptide_edges','peptide_supported','peptide_conflicted','peptide_packing_sum',
                'mhc_edges','mhc_conflicted','mhc_unanchored_conflicted']
    if set(required)-set(models.columns) or models['pdb.id'].is_duplicated().any():
        raise ValueError('unique models with complete components and neighborhood counts required')
    numeric = models.select(required[1:]).to_numpy()
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ValueError('finite nonnegative model components and counts required')
    if models.select(pl.any_horizontal(pl.col('iptm','topology','confidence_weight') > 1).any()).item():
        raise ValueError('prepared components must be bounded by 1')
    for interface in ('peptide','mhc'):
        for name in (('supported','conflicted','packing_sum') if interface=='peptide' else ('conflicted',)):
            if (models[interface+'_'+name] > models[interface+'_edges']+1e-12).any():
                raise ValueError('edge summaries exceed interface edge counts')
    if (models['mhc_unanchored_conflicted'] > models['mhc_conflicted']).any():
        raise ValueError('unanchored conflicts exceed MHC conflicts')
    def fraction(numerator, denominator, empty):
        return pl.when(pl.col(denominator)>0).then(pl.col(numerator)/pl.col(denominator)).otherwise(empty)
    pfree = 1-fraction('peptide_conflicted','peptide_edges',0.)
    mfree = 1-fraction('mhc_conflicted','mhc_edges',0.)
    quality = {'pep_support':fraction('peptide_supported','peptide_edges',0.),
               'pep_packing':fraction('peptide_packing_sum','peptide_edges',0.),
               'pep_conflict':pfree, 'mhc_conflict':mfree,
               'mhc_unanchored_conflict':1-fraction('mhc_unanchored_conflicted','mhc_edges',0.),
               'combined_conflict':pfree*mfree}
    confidence = pl.col('confidence_weight')*pl.col('iptm')
    topology = (1-pl.col('confidence_weight'))*pl.col('topology')
    result = models.with_columns((confidence+topology).alias('prepared'))
    result = result.with_columns([(confidence+topology*q).alias(name) for name,q in quality.items()])
    gates = {'gate_no_cdr3_peptide':pl.col('peptide_edges')>0,
             'gate_both_cdr3_pmhc':pl.col('pmhc_chains')==2,
             'gate_chemical_peptide':pl.col('peptide_supported')>0,
             'gate_both_chemical_pmhc':pl.col('chemical_pmhc_chains')==2}
    return result.with_columns([pl.when(condition).then(pl.col('prepared')).otherwise(0.).alias(name)
                                for name,condition in gates.items()])
