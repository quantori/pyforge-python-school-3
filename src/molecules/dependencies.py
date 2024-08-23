from functools import lru_cache
from src.molecules.schemas import PaginationQueryParams


@lru_cache
def get_pagination_query_params(page: int = 0, page_size: int = 1000):
    return PaginationQueryParams(page=page, page_size=page_size)
