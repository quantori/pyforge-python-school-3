import json

from pydantic import BaseModel


class CustomResponse(BaseModel):
    status: int
    reason: str = None
    headers: dict = None
    body: dict = None

    # create from a HTTPResponse object, pay attention to the types!!!
    @classmethod
    def from_http_response(cls, response):
        return CustomResponse(
            status=response.status,
            reason=response.reason,
            headers=dict(response.getheaders()),
            body=json.loads(response.read().decode('utf-8')) if response.getheader('Content-Type') == 'application/json'
            else None

        )

    class Config:
        description = "A custom response object. status is always present. reason, headers, and body are optional."
