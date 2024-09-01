
# Current Progress Overview

Currently main routes include:

- `/molecules` 
- `/drugs` 

These routes support basic CRUD operations(drugs is missing update).

Additionally, molecules support substructure search and file upload.

Everything is synchronous.

# Database Schema

Molecule groups are not implemented yet.

![img.png](img.png)


# API routes

## Molecules

    - POST: Create a new molecule

    - DELETE: /{molecule_id} Delete a molecule by id

    - PUT: /{molecule_id} Update a molecule by id

    - POST: /upload Upload a csv file containing molecules info

    - GET: /search?name={name} Search for a molecule by name

    - GET: /search/substructures?smiles={smiles}?limit={limit} Search for a molecule by substructure
    
    - GET: /search/superstructures?smiles={smiles}?limit={limit} Search for a molecule by superstructure






