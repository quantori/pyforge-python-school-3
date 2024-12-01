import pytest
from fastapi.testclient import TestClient
from unittest.mock import MagicMock, patch

from src.main import app, get_cached_result, set_cache
from src.database import get_db
from src.models import Molecule


@pytest.fixture
def mock_postgres(mocker):
    """Mock для SQLAlchemy сессии."""
    mock_session = mocker.MagicMock()
    return mock_session


@pytest.fixture
def mock_redis(mocker):
    """Mock для Redis клиента."""
    mock_redis_client = mocker.MagicMock()
    mocker.patch("src.main.redis_client", mock_redis_client)
    return mock_redis_client


@pytest.fixture
def fastapi_client(mock_postgres, mock_redis):
    """Создаёт клиент FastAPI с подменёнными зависимостями."""
    def override_get_db():
        yield mock_postgres

    app.dependency_overrides[get_db] = override_get_db
    return TestClient(app)


@pytest.fixture
def sample_molecule():
    sample_molecule = Molecule(
        id=23,
        name="Benzene",
        smiles="c1ccccc1",
        weight=78.11,
        formula="C6H6"
    )
    return sample_molecule


@pytest.fixture
def sample_molecules():
    """Фикстура для тестовых данных."""
    return [
        Molecule(id=1, name="Water", smiles="O", weight=18.015, formula="H2O"),
        Molecule(id=2, name="Ethanol", smiles="CCO", weight=46.07, formula="C2H6O"),
    ]


def test_get_molecule_by_id_with_cache(fastapi_client, mock_postgres, mock_redis, sample_molecule):

    # Мокаем запрос в базу данных
    mock_postgres.query.return_value.filter.return_value.first.return_value = sample_molecule

    # Мокаем Redis (нет кеша для первого вызова)
    mock_redis.get.return_value = None

    # Первый запрос (обращение к базе данных)
    response = fastapi_client.get(f"/molecule/{sample_molecule.id}")
    assert response.status_code == 200

    data = response.json()
    assert data["id"] == sample_molecule.id
    assert data["name"] == sample_molecule.name
    assert data["smiles"] == sample_molecule.smiles
    assert data["weight"] == sample_molecule.weight
    assert data["formula"] == sample_molecule.formula

    # Проверка, что данные были закешированы
    cache_key = f"molecule:{sample_molecule.id}"
    mock_redis.setex.assert_called_once_with(
        cache_key,
        60,  # Время кеширования
        '{"id": 23, "name": "Benzene", "smiles": "c1ccccc1", "weight": 78.11, "formula": "C6H6"}'
    )

    # Мокаем Redis (данные есть в кеше для второго вызова)
    mock_redis.get.return_value = '{"id": 23, "name": "Benzene", "smiles": "c1ccccc1", "weight": 78.11, "formula": "C6H6"}'

    # Второй запрос (данные берутся из кеша)
    response_cached = fastapi_client.get(f"/molecule/{sample_molecule.id}")
    assert response_cached.status_code == 200

    cached_data_from_response = response_cached.json()
    assert cached_data_from_response == data

    # Проверяем, что второй запрос не вызывал обращение к базе
    mock_postgres.query.return_value.filter.return_value.first.assert_called_once()


def test_get_all_molecules_with_cache(fastapi_client, mock_postgres, mock_redis, sample_molecules):
    # Мокаем запрос в базу данных
    mock_postgres.query.return_value.offset.return_value.limit.return_value.all.return_value = sample_molecules

    # Мокаем Redis (нет кеша для первого вызова)
    mock_redis.get.return_value = None

    # Первый запрос (обращение к базе данных)
    response = fastapi_client.get("/molecules_list?skip=0&limit=10")
    assert response.status_code == 200

    data = response.json()
    assert len(data) == 2
    for i, molecule in enumerate(sample_molecules):
        assert data[i]["id"] == molecule.id
        assert data[i]["name"] == molecule.name
        assert data[i]["smiles"] == molecule.smiles
        assert data[i]["weight"] == molecule.weight
        assert data[i]["formula"] == molecule.formula

    # Проверка, что данные были закешированы
    cache_key = "molecules_list:skip=0:limit=10"
    mock_redis.setex.assert_called_once_with(
        cache_key,
        60,  # Время кеширования
        '[{"id": 1, "name": "Water", "smiles": "O", "weight": 18.015, "formula": "H2O"}, '
        '{"id": 2, "name": "Ethanol", "smiles": "CCO", "weight": 46.07, "formula": "C2H6O"}]'
    )

    # Мокаем Redis (данные есть в кеше для второго вызова)
    mock_redis.get.return_value = '[{"id": 1, "name": "Water", "smiles": "O", "weight": 18.015, "formula": "H2O"}, ' \
                                  '{"id": 2, "name": "Ethanol", "smiles": "CCO", "weight": 46.07, "formula": "C2H6O"}]'

    # Второй запрос (данные берутся из кеша)
    response_cached = fastapi_client.get("/molecules_list?skip=0&limit=10")
    assert response_cached.status_code == 200

    cached_data_from_response = response_cached.json()
    assert cached_data_from_response == data

    # Проверяем, что второй запрос не вызывал обращение к базе
    mock_postgres.query.return_value.offset.return_value.limit.return_value.all.assert_called_once()
