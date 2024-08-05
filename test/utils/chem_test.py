from rdkit import Chem
from src.utils.chem import substructure_search, valid_smile

def test_valid_smile():
    assert valid_smile("CCO") == True, "CCO is a valid SMILES"
    assert valid_smile("C1CCCCC1") == True, "Cyclohexane is a valid SMILES"
    assert valid_smile("InvalidSMILES") == False, "InvalidSMILES is not a valid SMILES"

def test_substructure_search():
    mols = [
        {"smile": "CCO"},
        {"smile": "C1CCCCC1"},
        {"smile": "CCN"},
        {"smile": "CCOCC"},
    ]
    
    result = substructure_search(mols, "CCO")
    expected = [
        {"smile": "CCO"},
        {"smile": "CCOCC"},
    ]
    assert result == expected, f"Expected {expected}, but got {result}"
    
    result = substructure_search(mols, "C1CCCCC1")
    expected = [
        {"smile": "C1CCCCC1"},
    ]
    assert result == expected, f"Expected {expected}, but got {result}"
    
    result = substructure_search(mols, "CCN")
    expected = [
        {"smile": "CCN"},
    ]
    assert result == expected, f"Expected {expected}, but got {result}"