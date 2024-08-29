from sqlalchemy.orm import Session
from rdkit import Chem
import crud


def substructure(db: Session):
    # Fetch all molecules from the database
    molecules = crud.get_molecules(db)

    # Return an empty list if no molecules are found
    if not molecules:
        return []

    # Convert the list of molecule objects to a dictionary with identifiers as keys
    molecules_dict = {mol.identifier: mol for mol in molecules}

    # Create RDKit molecule objects for all molecules
    molecule_objects = {}
    for identifier, mol in molecules_dict.items():
        try:
            molecule_objects[identifier] = Chem.MolFromSmiles(mol.smile)
        except Exception as e:
            # Log or handle the invalid SMILES string here
            print(f"Error processing SMILES for {identifier}: {e}")
            molecule_objects[identifier] = None

    # Initialize results and substructure matches
    results = []
    substructure_matches = {identifier: [] for identifier in molecule_objects}

    # Check for substructure matches
    for identifier, mol in molecule_objects.items():
        if mol is None:
            continue

        for sub_id, sub_mol in molecule_objects.items():
            if identifier == sub_id or sub_mol is None:
                continue

            if mol.HasSubstructMatch(sub_mol):
                substructure_matches[identifier].append({
                    "identifier": sub_id,
                    "smile": molecules_dict[sub_id].smile
                })

    # Prepare the final result list
    for identifier, matches in substructure_matches.items():
        if matches:
            results.append({
                "identifier": identifier,
                "substructures": matches
            })

    return results