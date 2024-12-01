# This is a Celery configuration file

from celery import Celery

celery = Celery(
    'tasks',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0',
)

celery.autodiscover_tasks(['src']) # couldn't find the tasks without this parameter

celery.conf.update(task_track_started=True)