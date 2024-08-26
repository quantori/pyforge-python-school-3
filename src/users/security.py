import enum
from datetime import datetime, timezone, timedelta
import jwt
from fastapi.security import OAuth2PasswordBearer
from passlib.context import CryptContext
from src.configs import get_settings


def get_superuser_email_and_password():
    return get_settings().SUPERADMIN_EMAIL, get_settings().SUPERADMIN_PASSWORD


class Scope:
    MOLECULES_WRITE = "molecules:write"
    DRUGS_WRITE = "drugs:write"
    USERS_READ = "users:read"
    USERS_WRITE = "users:write"


pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
oauth2_scheme = OAuth2PasswordBearer(
    tokenUrl="/users/token",
    scopes={
        Scope.MOLECULES_WRITE: "Access to every molecule endpoint",
        Scope.DRUGS_WRITE: "Access to every drug endpoint",
        Scope.USERS_READ: "Access to read user details",
        Scope.USERS_WRITE: "Access to write user details",
    },
)


class Role(enum.Enum):
    SUPER_ADMIN = "SUPER_ADMIN"
    LAB_ADMIN = "LAB_ADMIN"
    HOSPITAL_ADMIN = "HOSPITAL_ADMIN"


role_scopes = {
    Role.SUPER_ADMIN: [
        Scope.MOLECULES_WRITE,
        Scope.DRUGS_WRITE,
        Scope.USERS_READ,
        Scope.USERS_WRITE,
    ],
    Role.LAB_ADMIN: [Scope.MOLECULES_WRITE, Scope.DRUGS_WRITE],
    Role.HOSPITAL_ADMIN: [],
}


def verify_password(plain_password, hashed_password):
    return pwd_context.verify(plain_password, hashed_password)


def get_password_hash(password):
    return pwd_context.hash(password)


def create_access_token(data: dict):
    to_encode = data.copy()
    expire = datetime.now(timezone.utc) + timedelta(
        minutes=get_settings().ACCESS_TOKEN_EXPIRE_MINUTES
    )
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(
        to_encode, get_settings().SECRET_KEY, algorithm=get_settings().ALGORITHM
    )
    return encoded_jwt


def decode_access_token(token):
    return jwt.decode(
        token, get_settings().SECRET_KEY, algorithms=[get_settings().ALGORITHM]
    )
