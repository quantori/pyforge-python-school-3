from src.utils.chem import substructure_search
from src.models.molecule import Molecule
from random import randint

molecules_table = []


def find_by_id(id: int):
    for m in molecules_table:
        if m["id"] == id:
            return m
    return None


def find_index_by_id(id: int):
    for i, m in enumerate(molecules_table):
        if m["id"] == id:
            return i
    return -1


def delete_by_id(id: int):
    molecule = find_by_id(id)
    if molecule:
        molecules_table.remove(molecule)
        return True
    else:
        return False


def get_all():
    return molecules_table


def get_filtered(substructre: str):
    return substructure_search(molecules_table, substructre)


def create(new_molecule: Molecule):
    molecule = new_molecule.model_dump()
    molecule["id"] = randint(1, 10000)
    molecules_table.append(molecule)
    return molecule


def update_by_id(id: int, new_molecule: Molecule):
    index = find_index_by_id(id)

    if index == -1:
        return None

    molecule = new_molecule.model_dump()
    molecule["id"] = id
    molecules_table[index] = molecule
    return molecule
