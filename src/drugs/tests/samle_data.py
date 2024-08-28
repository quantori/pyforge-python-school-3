from src.drugs.model import QuantityUnit
from src.drugs.schema import DrugMoleculeRequest, DrugRequest
from src.molecules.schemas import MoleculeRequest

caffeine_request = MoleculeRequest.model_validate(
    {"name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}
)

water_request = MoleculeRequest.model_validate({"name": "Water", "smiles": "O"})

sugar_request = MoleculeRequest.model_validate(
    {"name": "Sugar", "smiles": "C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O"}
)

coffe_request = DrugRequest(
    name="Coffe",
    description="Best drug ever",
    molecules=[
        DrugMoleculeRequest(molecule_id=1,
                            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C", quantity=0.2, quantity_unit=QuantityUnit.MASS
                            )
    ],
)
