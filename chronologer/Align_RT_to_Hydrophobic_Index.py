import sys, time, argparse, time, concurrent
import numpy as np, pandas as pd
import torch

import chronologer.src.constants as constants
from chronologer.src.local_io import read_rt_database
from chronologer.src.tensorize import codedseq_to_array
from chronologer.src.chronologer.model import initialize_chronologer_model
from chronologer.src.kde_alignment import KDE_align
from chronologer.src.utils import timer


def parse_args(args):
    parser = argparse.ArgumentParser( usage = 'Align_RT_to_Hydrophobic_Index.py [optional arguments] input.tsv source_name output.tsv',
                                      description = 'Align retention times into hydrophobic index space for '
                                                    'Chronologer (re-)training. Input files should be tab-'
                                                    'separated values with two columns: PeptideModSeq and RT. '
                                                    'A source name for the dataset must be provided (used as '
                                                    'part of the Chronologer loss function)',
                                     formatter_class = argparse.RawTextHelpFormatter, )
    parser.add_argument( 'input_file',
                         type = str,
                         help = 'Path to input file (tab-separated)')
    parser.add_argument( 'source',
                         type = str,
                         help = 'Database name for Chronologer')
    parser.add_argument( 'output_file',
                         type = str,
                         help = 'Path to output file (tab-separated)' )
    parser.add_argument( '--align_with_user_modifications',
                         action = 'store_true',
                         help = 'Consider custom user modifications when performing KDE alignment', )
    parser.add_argument( '--kde_n',
                         type = int,
                         default = 3000,
                         help = 'Number of points per axis in array for KDE alignment'
                                '\nDefault: 3000', ) 
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


def generate_RT_to_HI_spline( input_tsv, chronologer_model, device, include_user_mods, batch_size, kde_n, ):
    chronologer = initialize_chronologer_model( chronologer_model, ).to( device )
    df = read_rt_database( input_tsv, chronologer_compatible_only=True, )
    if not include_user_mods: # Exclude with coded seqs that contain numbers
        df = df[ [ c[1:-1].isalpha() for c in df.CodedPeptideSeq ] ]
    df = df.groupby('CodedPeptideSeq').mean().reset_index() # Ensure 1 entry per peptide
    #print( df.head(6) )
    seq_array = np.asarray( [ codedseq_to_array( p ) for p in df.CodedPeptideSeq ],
                              'int64', )
    pred_batches = [ ] 
    for i in range( 0, len(seq_array), batch_size ):
        peptide_tensor = torch.Tensor( seq_array[ i : i+batch_size, ] ).to( torch.int64 ).to( device )
        predictions = chronologer( peptide_tensor ).cpu().T[0].detach().numpy()
        pred_batches.append( predictions )
    df[ 'Pred_HI' ] = np.concatenate( pred_batches )
    spline, bounds, kde_map, fit_RT, fit_HI = KDE_align( np.array( df.RT ),
                                                         np.array( df.Pred_HI ),
                                                         kde_n, )
    return spline, bounds
        
def generate_output( input_tsv, output_tsv, source, spline, bounds, ):
    df = read_rt_database( input_tsv, chronologer_compatible_only=False, )
    df = df.drop( 'CodedPeptideSeq', axis=1, )
    clipped_rts = [ np.clip( rt, *bounds, ) for rt in df.RT ]
    df[ 'HI' ] = [ spline(rt) for rt in clipped_rts ]
    df[ 'Source' ] = source
    df.to_csv( output_tsv, sep='\t', header=True, index=False, )
    return 0



def main():
    start_time = time.time()
    args = parse_args(sys.argv[1:])
    
    device = 'cuda' if args.gpu else 'cpu'
    spline, bounds = generate_RT_to_HI_spline( args.input_file, 
                                               args.chronologer_model, 
                                               device, 
                                               args.align_with_user_modifications, 
                                               args.batch_size, 
                                               args.kde_n, )
    generate_output( args.input_file,
                     args.output_file,
                     args.source,
                     spline,
                     bounds, )
    
    print( 'KDE mapping to HI space complete (' + timer(start_time) + ')' )
    return 0

if __name__ == "__main__":
    main()
    sys.exit()
