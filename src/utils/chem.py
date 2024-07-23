from rdkit import Chem


def substructure_search(mols, mol):
    chem_mol = Chem.MolFromSmiles(mol)
    return [
        m for m in mols if Chem.MolFromSmiles(m["smile"]).HasSubstructMatch(chem_mol)
    ]


def valid_smile(str: str):
    mol = Chem.MolFromSmiles(str)
    if mol is None:
        return False
    return True


if __name__ == "__main__":
    pass
