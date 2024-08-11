import http.client
import json
from .schemas import CustomResponse


def _parse_base_url(url):
    # Assuming base_url is of the form "http://localhost:6969"
    if "http://" in url:
        url = url.replace("http://", "")
    host, port = url.split(":")
    return host, int(port)


class HTTPDatabaseClient:
    """The main method of this class is request, which can be used to make any HTTP request to the server.
    The other methods make this class a useful interface for interacting with the server."""

    def __init__(self, base_url: str):
        self.base_url = base_url.rstrip('/')
        self.host, self.port = _parse_base_url(self.base_url)

    def request(self, method, path, body=None) -> CustomResponse:
        """This method can be used to make any HTTP request to the server."""
        conn = http.client.HTTPConnection(self.host, self.port)
        headers = {'Content-Type': 'application/json'}
        json_body = json.dumps(body) if body else None
        conn.request(method, path, body=json_body, headers=headers)
        response = conn.getresponse()
        # when the connection is closed, the response attributs are no longer available, so we need to store them
        # in a custom response object
        custom_response = CustomResponse.from_http_response(response)
        conn.close()
        return custom_response

    def create_collection(self, collection_name: str) -> CustomResponse:
        return self.request("POST", "/collections", {"name": collection_name})

    def delete_collection(self, collection_name: str) -> CustomResponse:
        return self.request("DELETE", f"/collections/{collection_name}")

    def add_document(self, collection_name: str, document: dict) -> CustomResponse:
        return self.request("POST", f"/collections/{collection_name}/documents", document)

    def update_document(self, collection_name: str, document_id: str, document: dict) -> CustomResponse:
        return self.request("PUT", f"/collections/{collection_name}/documents/{document_id}", document)

    def delete_document(self, collection_name, document_id) -> CustomResponse:
        return self.request("DELETE", f"/collections/{collection_name}/documents/{document_id}")

    def get_documents(self, collection_name, field=None, value=None) -> CustomResponse:
        if field and value:
            return self.request("GET", f"/collections/{collection_name}/documents?field={field}&value={value}")
        return self.request("GET", f"/collections/{collection_name}/documents")

    def get_document(self, collection_name, document_id) -> CustomResponse:
        return self.request("GET", f"/collections/{collection_name}/documents/{document_id}")

    def clean_up(self) -> CustomResponse:
        return self.request("DELETE", "/collections")
