import logging
import pytest
from HTTP_database_client.src.database_client import HTTPDatabaseClient
import json

"""
This integration tests are supposed to be run against a mock server that simulates the behavior of the actual server.
Data is cleared after each test run.
"""


@pytest.fixture
def database_client() -> HTTPDatabaseClient:
    return HTTPDatabaseClient("http://localhost:6900")  # Assuming the mock server runs on port 6900


def test_home_page(database_client: HTTPDatabaseClient):
    response = database_client.request("GET", "/")
    assert response.status == 200


def test_create_collection(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    response = database_client.create_collection("asdasdasd")
    assert response.status == 201
    #     clean up
    response = database_client.clean_up()
    assert response.status == 200


def test_get_collection(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    response = database_client.create_collection("collection1")
    assert response.status == 201
    response = database_client.request("GET", "/collections/collection1")
    assert response.status == 200
    response = database_client.clean_up()
    assert response.status == 200


def test_add_document(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    database_client.create_collection("collection1")
    document = {"data": {"field1": "value1", "field2": "value2"}}
    response = database_client.add_document("collection1", document)
    print(response.body)
    assert "_document_id" in response.body
    added_document = database_client.get_document("collection1", response.body["_document_id"])
    assert response.status == 201
    assert document["data"] == added_document.body["data"]


def test_add_document_extra_fields(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    #     adding documents with extra fields other than the field data: dict
    # this should raise an error
    database_client.create_collection("collection1")
    document = {"field1": "value1", "field2": "value2", "extra_field": "extra_value"}
    response = database_client.add_document("collection1", document)
    assert response.status == 422


def test_clean_up(database_client: HTTPDatabaseClient):
    database_client.create_collection("collection1")
    database_client.add_document("collection1", {"data": {"field1": "value1", "field2": "value2"}})
    response = database_client.clean_up()
    assert response.status == 200
    #     check if the collection is deleted
    response = database_client.request("GET", "/collections/collection1")
    assert response.status == 404


def test_get_documents(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    database_client.create_collection("collection1")
    document1 = {"data": {"field1": "value1", "field2": "value2"}}
    document2 = {"data": {"field1": "value3", "field2": "value4"}}
    database_client.add_document("collection1", document1)
    database_client.add_document("collection1", document2)
    response = database_client.get_documents("collection1")
    documents = response.body["documents"]
    assert len(documents) == 2
    assert document1["data"] in [document["data"] for document in documents]
    assert document2["data"] in [document["data"] for document in documents]


#
#
def test_delete_collection(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    database_client.create_collection("collection1")
    document1 = {"data": {"field1": "value1", "field2": "value2"}}
    document2 = {"data": {"field1": "value3", "field2": "value4"}}
    database_client.add_document("collection1", document1)
    database_client.add_document("collection1", document2)
    response = database_client.delete_collection("collection1")
    assert response.status == 200
    #     check if the collection is deleted
    response = database_client.request("GET", "/collections/collection1")
    assert response.status == 404


def test_find_documents_by_field(database_client: HTTPDatabaseClient):
    database_client.clean_up()
    database_client.create_collection("molecules")
    database_client.add_document("molecules", {"data": {"name": "Methane", "smiles": "C"}})
    database_client.add_document("molecules", {"data": {"smiles": "CCO"}})
    database_client.add_document("molecules", {"data": {"name": "Ethanol", "smiles": "CCO"}})
    database_client.add_document("molecules", {"data": {"name": "Methane"}})
    database_client.add_document("molecules", {"data": {"description": "stupid document"}})
    response = database_client.get_documents("molecules", "name", "Methane")
    print(response.body)
    assert len(response.body["documents"]) == 2
    assert response.status == 200
    assert {"name": "Methane", "smiles": "C"} in [document["data"] for document in response.body["documents"]]
    assert {"name": "Methane"} in [document["data"] for document in response.body["documents"]]




