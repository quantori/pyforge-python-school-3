from src.molecules.tests.testing_utils import alkanes

# import csv writer
import csv


def generate_csv_file_alkanes():
    with open("alkanes.csv", mode="w") as file:
        writer = csv.writer(file)
        writer.writerow(["name", "smiles", "extra_column_are_ok"])
        for molecule in alkanes.values():
            writer.writerow([molecule["name"], molecule["smiles"]])


def generate_csv_file_invalid_header():
    with open("invalid_header.csv", mode="w") as file:
        writer = csv.writer(file)
        writer.writerow(["name", "invalid_column"])
        for molecule in alkanes.values():
            writer.writerow([molecule["name"], molecule["smiles"]])


def generate_csv_file_alkanes_decane_and_nonane_have_invalid_smiles():
    with open("decane_nonane_invalid_smiles.csv", mode="w") as file:
        writer = csv.writer(file)
        writer.writerow(["name", "smiles"])
        for molecule in alkanes.values():
            if molecule["name"] in ["Decane", "Nonane"]:
                writer.writerow([molecule["name"], "incontnentia"])
            else:
                writer.writerow([molecule["name"], molecule["smiles"]])


def generate_testing_files():
    generate_csv_file_alkanes()
    generate_csv_file_invalid_header()
    generate_csv_file_alkanes_decane_and_nonane_have_invalid_smiles()
