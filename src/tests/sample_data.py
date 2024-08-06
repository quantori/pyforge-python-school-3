from src.models.molecule_models import Molecule


def molecule_model_aspirin():
    return Molecule("CC(=O)Oc1ccccc1C(=O)O", molecule_name="Aspirin",
                    description="Aspirin is in a group of medications called salicylates... so on and so forth.")


def molecule_model_carbon_no_name_no_description():
    return Molecule("c")


def molecule_model_methane_no_description_custom_id():
    return Molecule("C", molecule_id=999, molecule_name="Methane")