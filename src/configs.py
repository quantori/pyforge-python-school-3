from functools import lru_cache

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):

    SECRET_KEY: str
    ALGORITHM: str
    ACCESS_TOKEN_EXPIRE_MINUTES: int

    DB_USER: str
    DB_PASSWORD: str
    DB_HOST: str
    DB_PORT: int
    DB_NAME: str
    TEST_MODE: bool
    TEST_DB_URL: str
    CONTEXT_PATH: str

    SUPERADMIN_EMAIL: str
    SUPERADMIN_PASSWORD: str

    model_config = SettingsConfigDict(
        env_file="/home/gaioz/quantori/pyforge-python-school-3/.env"
    )

    @property
    def database_url(self) -> str:
        if self.TEST_MODE:
            return self.TEST_DB_URL
        return f"postgresql://{self.DB_USER}:{self.DB_PASSWORD}@{self.DB_HOST}:{self.DB_PORT}/{self.DB_NAME}"


@lru_cache
def get_settings() -> Settings:
    return Settings()
