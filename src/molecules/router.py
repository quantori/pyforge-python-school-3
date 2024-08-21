from fastapi import APIRouter, HTTPException
from rdkit import Chem
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeResponse, MoleculeAdd, MoleculeUpdate


router = APIRouter()

#Add molecule (smiles) and its identifier
@router.post("/molecules", status_code=201, tags=["Molecules"], summary="Add new molecules to the DB", response_description="Molecule added successfully")
async def add_molecule(molecule: MoleculeAdd) -> dict:
    '''
    Create a molecule with the following details:

    **smiles**: Molecule in SMILES format  
    '''
    substructure = Chem.MolFromSmiles(molecule.smiles)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")
    molecule_data = molecule.model_dump()
    try:
        new_molecule_id = await MoleculeDAO.add_molecule(**molecule_data)
        return {"message": "The molecule is added!", "molecule": molecule}
    except HTTPException as e:
        raise e #(checks if smiles already exists in db and raises 409 error)
    except Exception:
        raise HTTPException(status_code=500, detail="Error adding the molecule")
    
#Get molecule by identifier
@router.get("/molecules/{molecule_id}",  tags=["Molecules"], summary="Retrieve molecule by ID", response_description="Molecule retrieved successfully")
async def get_molecule_by_id(molecule_id: int):
    '''
    Get molecule by identifier:

    **id**: Required
    '''
    rez = await MoleculeDAO.find_full_data(molecule_id=molecule_id)
    if rez is None:
        raise HTTPException(status_code=404, detail=f"Molecule with id {molecule_id} does not exist!")
    return rez

#Updating a molecule by identifier
@router.put("/molecules/{molecule_id}", tags=["Molecules"], summary="Update molecule by ID", response_model=MoleculeResponse)
async def update_molecule(molecule_id: int, updated_molecule: MoleculeUpdate):
    '''
    Updating a molecule by identifier:

    **id**: Required
    **smiles**: Molecule in SMILES format 
    '''
    substructure = Chem.MolFromSmiles(updated_molecule.smiles)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")
    try:
        updated = await MoleculeDAO.update(molecule_id, updated_molecule)
        
        if updated is None:
            raise HTTPException(status_code=404, detail="Molecule not found")
        
        return updated
    
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail="Unexpected error occurred")

#Delete a molecule by identifier
@router.delete("/molecules/{molecule_id}",  tags=["Molecules"], summary="Delete molecule by ID", response_description="Molecule Deleted")
async def delete_molecule(molecule_id: int) -> dict:
    '''
    Delete a molecule by identifier

    **id**: Required
    '''
    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        return {"message": f"The molecule with id {molecule_id} is deleted!"}
    else:
        raise HTTPException(status_code=404, detail=f"Molecule with id {molecule_id} does not exist!")


#List all molecules
@router.get("/molecules", tags=["Search"], summary="Retrieve all molecules", response_description="Molecules retrieved successfully")
async def get_all() -> list[MoleculeResponse]:
    return await MoleculeDAO.find_all_molecules()


#Substructure search for all added molecules
@router.get("/substructures", tags=["Search"], summary="Substructure search molecules", response_description="Substructure match molecules retrieved successfully")
async def substructure_search(smiles: str):
    '''
    Substructure search

    **substructure**: SMILES molecule required
    '''
    try:
        matches = await MoleculeDAO.substructure_search(smiles)
        if not matches:
            raise HTTPException(status_code=404, detail="No matching molecules found")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    return matches