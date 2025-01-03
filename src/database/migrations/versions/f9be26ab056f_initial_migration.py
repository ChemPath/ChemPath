"""Initial migration

Revision ID: f9be26ab056f
Revises: 
Create Date: 2024-09-10 15:29:38.769483

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'f9be26ab056f'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('compound_types')
    op.drop_table('traditional_knowledge')
    op.drop_table('literature_references')
    with op.batch_alter_table('plant_compounds', schema=None) as batch_op:
        batch_op.add_column(sa.Column('systematic_name', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('synonyms', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('cas_number', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('average_mass', sa.Float(), nullable=True))
        batch_op.add_column(sa.Column('monoisotopic_mass', sa.Float(), nullable=True))
        batch_op.add_column(sa.Column('chemical_formula', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('iupac_name', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('inchi_key', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('inchi_identifier', sa.String(), nullable=True))
        batch_op.add_column(sa.Column('solubility', sa.Float(), nullable=True))
        batch_op.add_column(sa.Column('log_s', sa.Float(), nullable=True))
        batch_op.add_column(sa.Column('log_p', sa.Float(), nullable=True))
        batch_op.alter_column('name',
               existing_type=sa.VARCHAR(),
               nullable=True)
        batch_op.drop_column('inchi')
        batch_op.drop_column('formula')

    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('plant_compounds', schema=None) as batch_op:
        batch_op.add_column(sa.Column('formula', sa.VARCHAR(), autoincrement=False, nullable=True))
        batch_op.add_column(sa.Column('inchi', sa.VARCHAR(), autoincrement=False, nullable=True))
        batch_op.alter_column('name',
               existing_type=sa.VARCHAR(),
               nullable=False)
        batch_op.drop_column('log_p')
        batch_op.drop_column('log_s')
        batch_op.drop_column('solubility')
        batch_op.drop_column('inchi_identifier')
        batch_op.drop_column('inchi_key')
        batch_op.drop_column('iupac_name')
        batch_op.drop_column('chemical_formula')
        batch_op.drop_column('monoisotopic_mass')
        batch_op.drop_column('average_mass')
        batch_op.drop_column('cas_number')
        batch_op.drop_column('synonyms')
        batch_op.drop_column('systematic_name')

    op.create_table('literature_references',
    sa.Column('id', sa.INTEGER(), autoincrement=True, nullable=False),
    sa.Column('pubmed_id', sa.VARCHAR(), autoincrement=False, nullable=False),
    sa.Column('title', sa.VARCHAR(), autoincrement=False, nullable=False),
    sa.Column('abstract', sa.TEXT(), autoincrement=False, nullable=True),
    sa.Column('authors', sa.VARCHAR(), autoincrement=False, nullable=True),
    sa.Column('publication_date', sa.DATE(), autoincrement=False, nullable=True),
    sa.Column('journal', sa.VARCHAR(), autoincrement=False, nullable=True),
    sa.Column('compound_id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.PrimaryKeyConstraint('id', name='literature_references_pkey'),
    sa.UniqueConstraint('pubmed_id', name='literature_references_pubmed_id_key')
    )
    op.create_table('traditional_knowledge',
    sa.Column('id', sa.INTEGER(), autoincrement=True, nullable=False),
    sa.Column('plant_compound_id', sa.INTEGER(), autoincrement=False, nullable=False),
    sa.Column('culture', sa.VARCHAR(length=100), autoincrement=False, nullable=False),
    sa.Column('use_description', sa.TEXT(), autoincrement=False, nullable=False),
    sa.Column('source', sa.VARCHAR(length=200), autoincrement=False, nullable=True),
    sa.PrimaryKeyConstraint('id', name='traditional_knowledge_pkey')
    )
    op.create_table('compound_types',
    sa.Column('id', sa.INTEGER(), autoincrement=True, nullable=False),
    sa.Column('name', sa.VARCHAR(length=50), autoincrement=False, nullable=False),
    sa.PrimaryKeyConstraint('id', name='compound_types_pkey'),
    sa.UniqueConstraint('name', name='compound_types_name_key')
    )
    # ### end Alembic commands ###
