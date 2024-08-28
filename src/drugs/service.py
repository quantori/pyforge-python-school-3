from src.drugs import mapper
from src.drugs.repository import DrugRepository
from src.drugs.schema import DrugRequest


class DrugService:

    def __init__(self, drug_repository: DrugRepository):
        self._drug_repository = drug_repository

    def save(self, drug: DrugRequest):
        drug, session = self._drug_repository.save(drug.model_dump())
        ans = mapper.drug_to_response(drug)
        session.close()
        return ans
