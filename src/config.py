from pydantic_settings import BaseSettings


class Config(BaseSettings):
    """Configuration settings for the application.
    Currently, the only setting is the database URL, that is supposed to be provided as an environment variable,
    to connect to the HTTP_database that I developed separately to this project.
    """

    DATABASE_URL: str
