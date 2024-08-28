from sqlalchemy.exc import IntegrityError

from src.drugs import mapper
from src.drugs.repository import DrugRepository
from src.drugs.schema import DrugRequest
from src.exceptions import BadRequestException, UnknownIdentifierException


class DrugService:

    def __init__(self, drug_repository: DrugRepository):
        self._drug_repository = drug_repository

    def save(self, drug: DrugRequest):
        try:
            drug, session = self._drug_repository.save(drug.model_dump())
            ans = mapper.drug_to_response(drug)
            session.close()
        except IntegrityError as e:
            """
            Here most probably if this exception is raised, it means that the molecule_id is not found in the database.
            But in case of any other exception, I will still provide other exception details to the user.
            """
            raise BadRequestException("Check if the molecules exist in the database. \n" + str(e))
        return ans

    def find_by_id(self, drug_id: int):
        drug = self._drug_repository.find_by_id(drug_id)
        if drug is None:
            raise UnknownIdentifierException(drug_id)
        return mapper.drug_to_response(drug)

