import polars as pl
import pytest

from tcren.contacts.language import contact_language


def test_language_uses_contact_and_noncontact_shell_observations():
    opportunities = pl.DataFrame({
        "pdb.id": ["m1"] * 3, "interface": ["tcr_peptide"] * 3,
        "chain.id.from": ["T"] * 3, "chain.id.to": ["P"] * 3,
        "chain.type.from": ["TRA"] * 3, "chain.type.to": ["PEPTIDE"] * 3,
        "residue.index.from": [0, 0, 0], "residue.index.to": [0, 1, 2],
        "aa.from": ["A"] * 3, "aa.to": ["G", "D", "A"],
        "atom.from": ["CA", None, "CB"], "atom.to": ["N", None, "O"],
        "region.type.from": ["CDR3"] * 3,
        "contacted": [True, False, True], "within_shell": [True, True, False],
        "n_atom_pairs": [2, 0, 1], "n_heavy_from": [4] * 3, "n_heavy_to": [4] * 3,
    })
    effects = pl.DataFrame({
        "interface": ["tcr_peptide", "tcr_peptide"],
        "block": ["A_engagement", "A_engagement"],
        "cell": ["A|G", "A|D"],
        "effect": [1., 6.], "complement": [0., -2.],
    })
    pooled = {"intercept": .2, "tcr_peptide|A_engagement|mean": 2., "sweep_max": .5}
    intercept_only = {"intercept": .8}
    result = contact_language(opportunities, [(effects, pooled), (effects, intercept_only)])
    assert result["id"].to_list() == ["m1"]
    assert result["contact_language"][0] == pytest.approx(.25)
    assert result["n_language_models"][0] == 2


def test_contact_language_reference_is_content_identified(tmp_path):
    from hashlib import sha256
    import json

    from tcren.contacts.language import load_contact_language_reference

    path = tmp_path / "reference.json"
    payload = {
        "format": "tcren-contact-language-1",
        "version": "native387-contact-language-v1",
        "models": [{
            "name": "synthetic",
            "coefficients": {"intercept": 0.5},
            "effects": [["tcr_peptide", "A_engagement", "A|C", 0.3, -0.2]],
        }],
        "atom_effects": {"tcr_peptide/CA|O": 0.75},
        "residue_effects": {"tcr_peptide/A|C": 0.25},
    }
    path.write_text(json.dumps(payload))
    reference = load_contact_language_reference(path)
    assert reference.sha256 == sha256(path.read_bytes()).hexdigest()
    assert reference.models[0][0]["cell"].to_list() == ["A|C"]
    assert reference.models[0][1] == {"intercept": 0.5}
    assert reference.atom_effects == payload["atom_effects"]
    with pytest.raises(ValueError, match="SHA-256"):
        load_contact_language_reference(path, expected_sha256="0" * 64)
