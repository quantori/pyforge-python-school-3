import json

from pydantic import BaseModel


class CustomResponse(BaseModel):
    status: int
    reason: str = ""
    headers: dict = {}
    body: dict = {}

    # create from a HTTPResponse object, pay attention to the types!!!
    @classmethod
    def from_http_response(cls, response):
        body = json.loads(response.read().decode())
        return cls(status=response.status, reason=response.reason, headers=dict(response.getheaders()),
                   body=body)

    class Config:
        description = "A custom response object. status is always present. reason, headers, and body are optional."