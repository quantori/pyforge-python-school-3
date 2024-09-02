import json
from sqlalchemy import delete, update as sqlalchemy_update
from sqlalchemy.future import select
from sqlalchemy.exc import SQLAlchemyError
from molecules.models import Molecule
from database import async_session_maker
from dao_dir.base import BaseDAO
from rdkit import Chem
from typing import List, Dict, AsyncIterator


class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(
        cls,
        limit: int = 100,
        offset: int = 0
    ) -> List[Dict]:
        async with async_session_maker() as session:
            query = select(cls.model).limit(limit).offset(offset)
            result = await session.execute(query)
            return [
                cls._model_to_dict(mol) for mol in result.scalars().all()
            ]

    @classmethod
    async def find_all_molecules_iterator(
        cls,
        limit: int
    ) -> AsyncIterator[Dict]:
        offset = 0
        while True:
            molecules_group = await cls.find_all_molecules(limit, offset)
            if not molecules_group:
                break
            for molecule in molecules_group:
                yield molecule
            offset += limit

    @classmethod
    async def find_full_data(cls, mol_id: int) -> Dict:
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(id=mol_id)
            result = await session.execute(query)
            mol_info = result.scalar_one_or_none()
            if not mol_info:
                return None
            return cls._model_to_dict(mol_info)

    @staticmethod
    def _model_to_dict(model_instance) -> Dict:
        if model_instance is None:
            return None
        return {
            column.name: getattr(model_instance, column.name)
            for column in model_instance.__table__.columns
        }

    @classmethod
    async def add_mol(cls, **mol_data: dict) -> int:
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
                query = (
                    sqlalchemy_update(cls.model)
                    .where(cls.model.id == mol_id)
                    .values(name=name)
                )
                result = await session.execute(query)
                await session.commit()
                if result.rowcount == 0:
                    raise ValueError("Molecule not found")

    @classmethod
    async def delete_mol_by_id(cls, mol_id: int) -> int:
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=mol_id)
                result = await session.execute(query)
                mol_to_delete = result.scalar_one_or_none()
                if not mol_to_delete:
                    raise ValueError("Molecule not found")
                await session.execute(delete(cls.model).filter_by(id=mol_id))
                await session.commit()
                return mol_id

    @classmethod
    async def find_by_substructure(
        cls,
        substructure_smiles: str,
        limit: int = 100,
        offset: int = 0
    ) -> List[Dict]:
        if not substructure_smiles:
            raise ValueError("Substructure SMILES string cannot be empty")
        substructure_mol = Chem.MolFromSmiles(substructure_smiles)
        if substructure_mol is None:
            raise ValueError("Invalid substructure SMILES string")

        async with async_session_maker() as session:
            try:
                query = select(cls.model).offset(offset).limit(limit)
                result = await session.execute(query)
                molecules = result.scalars().all()

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
    async def find_by_substructure_iterator(
        cls,
        substructure_smiles: str,
        limit: int
    ) -> AsyncIterator[Dict]:
        if not substructure_smiles:
            raise ValueError("Substructure SMILES string cannot be empty")
        substructure_mol = Chem.MolFromSmiles(substructure_smiles)
        if substructure_mol is None:
            raise ValueError("Invalid substructure SMILES string")

        async with async_session_maker() as session:
            offset = 0
            while True:
                try:
                    query = select(cls.model).offset(offset).limit(limit)
                    result = await session.execute(query)
                    molecules_batch = result.scalars().all()

                    if not molecules_batch:
                        break

                    matches = []
                    for molecule in molecules_batch:
                        mol = Chem.MolFromSmiles(molecule.name)
                        if mol and mol.HasSubstructMatch(substructure_mol):
                            mol_data = {
                                "id": molecule.id,
                                "name": molecule.name,
                            }
                            matches.append(mol_data)
                    if matches:
                        for match in matches:
                            yield match
                    offset += limit

                except Exception:
                    raise

    @classmethod
    async def upload_file(cls, file_content: str) -> int:
        try:
            molecules = json.loads(file_content)
            if not isinstance(molecules, list):
                raise ValueError(
                    "Invalid file format: expected a list of molecules"
                )

            added_count = 0

            async with async_session_maker() as session:
                async with session.begin():
                    for molecule in molecules:
                        mol_id = molecule.get("mol_id")
                        name = molecule.get("name")

                        if not mol_id or not name:
                            continue

                        if await cls.find_full_data(mol_id):
                            continue

                        mol = Chem.MolFromSmiles(name)
                        if not mol:
                            continue

                        new_mol = cls.model(id=mol_id, name=name)
                        session.add(new_mol)
                        added_count += 1

                    await session.flush()
                    await session.commit()

            return added_count

        except json.JSONDecodeError as e:
            raise ValueError(
                "Invalid file format: unable to decode JSON"
            ) from e
        except SQLAlchemyError as e:
            raise Exception("Database error occurred") from e
        except Exception as e:
            raise Exception(
                f"An error occurred while processing the file: {str(e)}"
            ) from e
