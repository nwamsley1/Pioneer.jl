import sys, time, argparse, time, concurrent
import numpy as np, pandas as pd
import torch
from tqdm import tqdm
import src.constants as constants
from src.local_io import read_rt_database
from src.tensorize import codedseq_to_array
from src.chronologer.model import initialize_chronologer_model
from src.kde_alignment import KDE_align
from src.utils import timer


def parse_args(args):
    parser = argparse.ArgumentParser( usage = 'Predict_RT.py [optional arguments] input.tsv source_name output.tsv',
                                      description = 'Predict RTs using Chronologer. Input file should be tab-'
                                                    'separated values with at least one column PeptideModSeq ('
                                                    'additional columns will be preserved) and a new Pred_HI column '
                                                    'added with predictions. Incompatible peptides will be dropped.',
                                     formatter_class = argparse.RawTextHelpFormatter, )
    parser.add_argument( 'input_file',
                         type = str,
                         help = 'Path to input file (tab-separated)')
    parser.add_argument( 'output_file',
                         type = str,
                         help = 'Path to output file (tab-separated)' )
    parser.add_argument( '--gpu',
                         action = 'store_true',
                         help = 'Use GPU (CUDA) for Cartographer/Chronologer computations', )
    parser.add_argument( '--chronologer_model',
                         type = str,
                         default = constants.default_chronologer_model_file,
                         help = 'Path to Chronologer model file'
                                '\nDefault: ' + constants.default_chronologer_model_file, )
    parser.add_argument( '--batch_size',
                         type = int,
                         default = 2048,
                         #choices = range( 1, 10241, ),
                         help = 'Number of peptides per computation batch'
                                '\nDefault: 2048', )    
    return parser.parse_args(args)


def generate_df_with_pred( input_tsv, output_tsv, chronologer_model, device, batch_size, ):
    chronologer = initialize_chronologer_model( chronologer_model, ).to( device )
    df = read_rt_database( input_tsv, chronologer_compatible_only=True, )
    seq_array = np.asarray( [ codedseq_to_array( p ) for p in df.CodedPeptideSeq ],
                              'int64', )
    pred_batches = [ ] 
    for i in tqdm(range( 0, len(seq_array), batch_size )):
        peptide_tensor = torch.Tensor( seq_array[ i : i+batch_size, ] ).to( torch.int64 ).to( device )
        predictions = chronologer( peptide_tensor ).cpu().T[0].detach().numpy()
        pred_batches.append( predictions )
    df['iRT'] = np.concatenate( pred_batches )
    df.drop('CodedPeptideSeq', axis = 1, inplace = True)
    df.drop('PeptideModSeq', axis = 1, inplace = True)
    df.drop('chronologer_sequence', axis = 1, inplace = True)
    df.drop('PeptideLength', axis = 1, inplace = True)
    df.to_csv(output_tsv, sep='\t', header=True, index=False, )
    return df
        

def main():
    start_time = time.time()
    args = parse_args(sys.argv[1:])
    
    device = 'cuda' if args.gpu else 'cpu'

    generate_df_with_pred( args.input_file,
                           args.output_file,
                           args.chronologer_model, 
                           device, 
                           args.batch_size, )
    
    print( 'Prediction of retention times completed (' + timer(start_time) + ')' )
    return 0

if __name__ == "__main__":
    main()
    sys.exit()
