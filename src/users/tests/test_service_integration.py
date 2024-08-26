import pytest

from src.database import Base
from src.users.exceptions import DuplicateEmailException
from src.users.schemas import RegisterRequest
from src.database import get_session_factory
from src.users.service import get_user_service
from src.users.repository import get_user_repository
from src.users.tests.sample_data import lab_admin_register_request_data, hospital_admin_register_request_data


@pytest.fixture
def init_db():
    Base.metadata.drop_all(get_session_factory()().get_bind())
    Base.metadata.create_all(get_session_factory()().get_bind())
    yield
    Base.metadata.drop_all(get_session_factory()().get_bind())


@pytest.fixture
def user_service():
    return get_user_service(get_user_repository(get_session_factory()))


def test_test(init_db, user_service):
    assert True


def test_register_lab_admin(init_db, user_service):
    register_request = RegisterRequest.model_validate(lab_admin_register_request_data)
    response = user_service.register(register_request)
    assert response.email == lab_admin_register_request_data["email"]

    find_by_email_response = user_service.find_by_email(lab_admin_register_request_data["email"])
    print(find_by_email_response)

    assert find_by_email_response.email == lab_admin_register_request_data["email"]
    assert find_by_email_response.full_name == lab_admin_register_request_data["full_name"]
    assert find_by_email_response.role.value == lab_admin_register_request_data["role"]
    assert find_by_email_response.is_active


def test_register_twice(init_db, user_service):
    register_request = RegisterRequest.model_validate(lab_admin_register_request_data)
    response = user_service.register(register_request)
    assert response.email == lab_admin_register_request_data["email"]

    with pytest.raises(DuplicateEmailException):
        user_service.register(register_request)


def test_get_by_id(init_db, user_service):
    register_request = RegisterRequest.model_validate(lab_admin_register_request_data)
    response = user_service.register(register_request)
    assert response.email == lab_admin_register_request_data["email"]

    find_by_email_response = user_service.find_by_email(lab_admin_register_request_data["email"])
    assert find_by_email_response.email == lab_admin_register_request_data["email"]

    user_id = find_by_email_response.user_id

    find_by_id_response = user_service.find_by_id(user_id)
    assert find_by_id_response.email == lab_admin_register_request_data["email"]
    assert find_by_id_response.full_name == lab_admin_register_request_data["full_name"]
    assert find_by_id_response.role.value == lab_admin_register_request_data["role"]
    assert find_by_id_response.is_active


def test_get_all_users(init_db, user_service):
    register_request = RegisterRequest.model_validate(lab_admin_register_request_data)
    response = user_service.register(register_request)
    assert response.email == lab_admin_register_request_data["email"]

    register_request1 = RegisterRequest.model_validate(hospital_admin_register_request_data)
    response1 = user_service.register(register_request1)
    assert response1.email == hospital_admin_register_request_data["email"]

    all_users = user_service.find_all()
    assert len(all_users) == 2
    assert all_users[0].email == lab_admin_register_request_data["email"]
    assert all_users[1].email == hospital_admin_register_request_data["email"]



