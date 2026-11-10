import hashlib
import io
import tarfile

from tcren.complementarity.command import _structures
from tcren.structure.model import Structure


def test_archive_member_identities_do_not_collapse_or_become_output_paths(tmp_path, monkeypatch):
    archive = tmp_path / "models.tar.gz"
    names = ["a/b.pdb", "a__b.pdb", "../../outside.pdb"]
    with tarfile.open(archive, "w:gz") as output:
        for name in names:
            info = tarfile.TarInfo(name)
            info.size = 6
            output.addfile(info, io.BytesIO(b"record"))
    observed = []

    def parse(path, pdb_id):
        observed.append(path)
        return Structure(pdb_id, [])

    monkeypatch.setattr("tcren.structure.parse_structure", parse)
    structures = list(_structures(archive, tmp_path))
    assert [s.pdb_id for s in structures] == ["a/b", "a__b", "../../outside"]
    assert len(set(s.pdb_id for s in structures)) == 3
    assert all(str(tmp_path) in path for path in observed)
    assert all(s.annotation_provenance["visualization_input"]["sha256"] == hashlib.sha256(b"record").hexdigest()
               for s in structures)
    assert not (tmp_path.parent / "outside.pdb").exists()
