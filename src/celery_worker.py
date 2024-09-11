from celery import Celery
from dotenv import load_dotenv
import os

load_dotenv(".env")

REDIS_URL = os.getenv('REDIS_URL')
if not REDIS_URL:
    raise ValueError("REDIS_URL is not set in the environment variables")


celery = Celery(
    "tasks",
    broker=REDIS_URL,
    backend=REDIS_URL
)

celery.conf.update(
    task_routes={
        'tasks.substructure_search_task': {'queue': 'search_queue'},
    }
)

try:
    from redis import Redis
    redis_client = Redis.from_url(REDIS_URL)
    redis_client.ping()
    print("Redis connection successful.")
except Exception as e:
    print(f"Failed to connect to Redis: {e}")
    raise