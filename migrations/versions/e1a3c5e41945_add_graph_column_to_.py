"""Add graph column to AdvancedRetrosynthesisResult

Revision ID: e1a3c5e41945
Revises: 
Create Date: 2024-08-04 20:07:33.694242

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = 'e1a3c5e41945'
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass
