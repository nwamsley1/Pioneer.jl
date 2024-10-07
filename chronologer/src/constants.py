import os


epsilon = 1e-7
train_fdr = 0.01

min_peptide_len = 6
max_peptide_len = 50


src_dir = os.path.dirname(os.path.abspath(__file__))
chronologer_db_loc = os.path.join( src_dir, '..', 'data', 'Chronologer_DB_220308.gz' )

default_chronologer_model_file = os.path.join( src_dir, '..', 'models', 'Chronologer_20220601193755.pt' )

seed = 2447

validation_fraction = 0.2

