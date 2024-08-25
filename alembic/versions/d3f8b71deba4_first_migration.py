"""First migration

Revision ID: d3f8b71deba4
Revises: 
Create Date: 2024-08-23 18:04:19.926772

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision: str = 'd3f8b71deba4'
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        'molecules',
        sa.Column('id', sa.Integer, primary_key=True),
        sa.Column('name', sa.String(50), nullable=False),
    )


def downgrade() -> None:
    op.drop_table('molecules')
