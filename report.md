# Report for the current state of the project

## Introduction

This report is a summary of the current state of the project. I'll try to explain everything simply.

I tried to implement pratices from fastapis official documentation. There was a nice tutorial there. Also, I had to 
look up to the documentations of pydantic and sqlalchemy latest versions. 

Everything in api is synchronous now. Despite the lesson presented in the lecture 8 by Yaroslav Korkodinov, I decided not to use asyncronous database. I think
at this point, I would rather focus on the basics. I'll definitely try to implement asyncronous database soon.

## Roles

    - admin

    - physician

    - patient

    - none


## endpoints



###  `api/molecules/`
Expose this endpoint for curious people. Maybe it can be useful for physicians, but mostly for 
curious people, who want to know what substances are there in the medices they consume, or maybe to
learn more about molecules


##### **for drug_admin**:

    - POST: Create a new molecule

    - DELETE: /{molecule_id} Delete a molecule by id

    - PUT: /{molecule_id} Update a molecule by id

    - GET: /{molecule_id} Get a molecule by id

    - POST: /upload Upload a csv file containing molecules info

##### **for everyone**:

    - GET: /?page={page}&page_size={page_size}?name={name} Get all molecules, you can filter by name

    - GET: /search/smiles_containing_substructure?substructure={substructure_smiles} Get all molecules containing a given substructure

    - GET: /search/substructures_of?smiles={smiles} Get all substructures of a molecule

### `api/drugs/`

##### **for drug_admin**:
  
    - GET: /{drug_id} Get a drug by id

    - POST: Create a new drug

    - DELETE: /{drug_id} Delete a drug by id

    - PUT: /{drug_id} Update a drug by id

    
###### **for everyone**

    - GET: /?page={page}&page_size={page_size}?name={drug_name} Get all drugs, you can filter by name

### `api/patients/`
    
    - GET: /{patient_id} Get a patient by id
    
    - POST: Create a new patient

    - DELETE: /{patient_id} Delete a patient by id

    - PUT: /{patient_id} Update a patient by id




