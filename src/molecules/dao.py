from sqlalchemy import delete
from sqlalchemy.future import select
from molecules.models import Molecule
from database import async_session_maker
from dao.base import BaseDAO


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
