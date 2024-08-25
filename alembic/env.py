import asyncio
from logging.config import fileConfig
from alembic import context
from sqlalchemy.ext.asyncio import create_async_engine
import os
from os.path import dirname, abspath
import sys
from src.database import Base

sys.path.insert(0, dirname(dirname(abspath(__file__))))
# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config


db_url = os.getenv("DB_URL")

if db_url:
    config.set_main_option("sqlalchemy.url", db_url)
# Interpret the config file for Python logging.
# This line sets up loggers basically.
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# add your model's MetaData object here
# for 'autogenerate' support
# from myapp import mymodel
# target_metadata = mymodel.Base.metadata
target_metadata = Base.metadata

# other values from the config, defined by the needs of env.py,
# can be acquired:
# my_important_option = config.get_main_option("my_important_option")
# ... etc.


async def run_async_migrations() -> None:
    connectable = create_async_engine(
        config.get_main_option("sqlalchemy.url"), 
        future=True
    )

    async with connectable.connect() as connection:
        await connection.run_sync(lambda conn: context.configure(
            connection=conn,
            target_metadata=target_metadata,
            process_bind_parameters=True,
            include_schemas=True,
        ))

        async with connection.begin():
            await connection.run_sync(lambda conn: context.run_migrations())


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.
    """
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    asyncio.run(run_async_migrations())


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
    