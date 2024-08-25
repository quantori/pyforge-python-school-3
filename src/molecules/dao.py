import json
from sqlalchemy import delete
from sqlalchemy import update as sqlalchemy_update
from sqlalchemy.future import select
from molecules.models import Molecule
from database import async_session_maker
from dao_dir.base import BaseDAO
from rdkit import Chem
from sqlalchemy.exc import SQLAlchemyError
from typing import List, Dict


class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(cls):
        async with async_session_maker() as session:
            query = select(cls.model)
            molecules = await session.execute(query)
            return molecules.scalars().all()

    @classmethod
    async def find_full_data(cls, mol_id):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(id=mol_id)
            result = await session.execute(query)
            mol_info = result.scalar_one_or_none()

            if not mol_info:
                return None

            mol_data = mol_info.to_dict()
            return mol_data

    @classmethod
    async def add_mol(cls, **mol_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                new_mol = cls.model(**mol_data)
                session.add(new_mol)
                await session.flush()
                new_mol_id = new_mol.id
                await session.commit()
                return new_mol_id
            
    @classmethod
    async def update_mol(cls, mol_id: int, name: str):
        async with async_session_maker() as session:
            async with session.begin():
                query = sqlalchemy_update(cls.model).where(cls.model.id == mol_id).values(name=name)
                result = await session.execute(query)
                await session.commit()
                if result.rowcount == 0:
                    raise ValueError("Molecule not found")
                
    @classmethod
    async def delete_mol_by_id(cls, mol_id: int):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=mol_id)
                result = await session.execute(query)
                mol_to_delete = result.scalar_one_or_none()

                if not mol_to_delete:
                    return None

                await session.execute(delete(cls.model).filter_by(id=mol_id))

                await session.commit()
                return mol_id
            
    @classmethod
    async def find_by_substructure(cls, substructure_smiles: str) -> List[Dict]:
        async with async_session_maker() as session:
            try:
                query = select(cls.model)
                result = await session.execute(query)
                molecules = result.scalars().all()

                if not substructure_smiles:
                    raise ValueError("Substructure SMILES string cannot be empty")

                substructure_mol = Chem.MolFromSmiles(substructure_smiles)
                if substructure_mol is None:
                    raise ValueError("Invalid substructure SMILES string")

                matches = []
                for molecule in molecules:
                    mol = Chem.MolFromSmiles(molecule.name)
                    if mol and mol.HasSubstructMatch(substructure_mol):
                        mol_data = {
                            "id": molecule.id,
                            "name": molecule.name,
                        }
                        matches.append(mol_data)

                return matches
            except SQLAlchemyError as e:
                raise Exception("Database error occurred") from e
            
    @classmethod
    async def upload_file(cls, file_content: str) -> int:
        """
        Parse the JSON content of the file and add molecules to the database.
        """
        try:
            molecules = json.loads(file_content)
            added_count = 0

            async with async_session_maker() as session:
                async with session.begin():
                    for molecule in molecules:
                        mol_id = molecule.get("mol_id")
                        name = molecule.get("name")
                        if await cls.find_molecule_by_id(mol_id):
                            continue
                        if not Chem.MolFromSmiles(name):
                            continue
                        new_mol = cls.model(id=mol_id, name=name)
                        session.add(new_mol)
                        added_count += 1
                    await session.flush()
                    await session.commit()

            return added_count

        except json.JSONDecodeError as e:
            raise ValueError("Invalid file format") from e
        except SQLAlchemyError as e:
            raise Exception("Database error occurred") from e