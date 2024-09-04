from functools import lru_cache

from sqlalchemy import text, and_, select
from sqlmodel import Session

from src.molecules.model import Molecule
from src.molecules.schema import SearchParams
from src.repository import SQLAlchemyRepository


class MoleculeRepository(SQLAlchemyRepository):
    def __init__(self):
        super().__init__(Molecule)

    def find_all(self, session: Session, page=0, page_size=1000, search_params: SearchParams = None):
        """
        If name is provided, then fuzzy search with trigrams is performed, results are ordered by similarity,
        notice that in that case, order_by and order are ignored. Filtering by mass is still possible.

        If name is not provided, then normal search is performed, results are ordered by order_by and order.

        Return type for now is a collection of ScalarResult, might include extra field "sml" for similarity,
        so be careful when using this method.

        //TODO: This append is not a good practice, but I will keep it for now, I saw there are nice patterns
        //TODO: in sql alchemy querying guide, I will try to implement them later.
        """
        query_builder = []

        if search_params.name:
            query_builder.append(
                f"""
                select similarity('{search_params.name}', m.name) as sml,
                       m.molecule_id,
                       m.name,
                       m.smiles,
                       m.mass,
                       m.created_at,
                       m.updated_at
                from molecules m
                where m.name % '{search_params.name}'
                """
            )
            if search_params.min_mass:
                query_builder.append(f" and {search_params.min_mass} <= m.mass")
            if search_params.max_mass:
                query_builder.append(f" and {search_params.max_mass} >= m.mass")

            query_builder.append(" order by sml desc")

        else:
            query_builder.append("select * from molecules where 1 = 1")

            if search_params.min_mass:
                query_builder.append(f" and {search_params.min_mass} <= mass")
            if search_params.max_mass:
                query_builder.append(f" and {search_params.max_mass} >= mass")

            if search_params.order_by:
                query_builder.append(f" order by {search_params.order_by} {search_params.order}")

        query_builder.append(f" limit {page_size} offset {page * page_size}")

        query = ''.join(query_builder)
        return session.execute(text(query))


@lru_cache
def get_molecule_repository():
    return MoleculeRepository()
