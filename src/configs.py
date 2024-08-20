from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    DB_USER: str
    DB_PASSWORD: str
    DB_HOST: str
    DB_PORT: int
    DB_NAME: str
    TEST_MODE: bool
    TEST_DB_URL: str

    model_config = SettingsConfigDict(env_file=".env")

    @property
    def database_url(self) -> str:
        if self.TEST_MODE:
            return self.TEST_DB_URL
        return f"postgresql://{self.DB_USER}:{self.DB_PASSWORD}@{self.DB_HOST}:{self.DB_PORT}/{self.DB_NAME}"
