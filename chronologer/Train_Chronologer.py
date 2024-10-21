
## Libraries
import os, sys, argparse, time
from datetime import datetime
import numpy as np
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import chronologer.src.constants as constants 
from chronologer.src.chronologer.training_functions import train_chronologer 

def parse_args(args):
    src_dir = os.path.dirname(os.path.abspath(__file__))
    timestamp = datetime.now().strftime( '%Y%m%d%H%M%S' )
    default_out_filename = 'Chronologer_'+timestamp+'.pt'
    parser = argparse.ArgumentParser( usage = 'Train_Chronologer.py [optional arguments]',
                                      ###
                                      ### DESCRIPTION NEEDS TO BE UPDATED 
                                      ###
                                      description = 'Train Chronologer using Chronologer DB -/+ user data in a'
                                                    'properly formatted user TSV file. By default, training will '
                                                    'begin from the published Chronologer model and is primarily '
                                                    'intended for augmented training with new data. For training '
                                                    'from scratch, the --start_model parameter can be set to None.',
                                     formatter_class = argparse.RawTextHelpFormatter, )
    parser.add_argument( '--output_file',
                         type = str,
                         help = 'Model filename',
                         default = default_out_filename, )
    parser.add_argument( '--output_dir',
                         type = str,
                         help = 'Directory to save model'
                                '\n\tDefault: chronologer/models',
                         default = os.path.join(src_dir,'models'), )
    parser.add_argument( '--exclude_chronologer_db',
                         action = 'store_true',
                         help = 'Do not include 2M+ peptide Chronologer database during training (not recommended)' )
    parser.add_argument( '--training_data',
                         type = str,
                         nargs = '+',
                         default = [ ],
                         help = 'Path(s) to TSV file(s) with additional training data'
                                '\n\tDefault: None', )
    parser.add_argument( '--start_model',
                         type = str,
                         help = 'Pre-trained chronologer model for training initialization' +
                                '\n\tDefault: '+constants.default_chronologer_model_file,
                         default = constants.default_chronologer_model_file, )
    parser.add_argument( '--validation_fraction',
                         type = float,
                         help = 'Fraction of data to withhold for model evaluation' +
                                '\n\tDefault: ' + format( constants.validation_fraction,'.2f' ),
                         default = constants.validation_fraction, )
    parser.add_argument( '--lr_weight',
                         type = float,
                         help = 'Scaler weight applied to the learning rate during training '
                                '(useful to set to 0.5 for re-training)'
                                '\n\tDefault: 1.0',
                         default = 1.0, )
    parser.add_argument( '--seed',
                         type = int,
                         help = 'Random seed for splitting training vs evaluation data' +
                                '\n\tDefault: ' + str(constants.seed),
                         default = constants.seed, )
    parser.add_argument( '--gpu',
                         action = 'store_true',
                         help = 'Use GPU (CUDA) for Chronologer training', )
    return parser.parse_args(args)


def main():
    start_time = time.time()
    args = parse_args(sys.argv[1:])
    model_out_file = os.path.join( args.output_dir, args.output_file, )
    if not args.exclude_chronologer_db: 
        data = [ constants.chronologer_db_loc ] + args.training_data
    else:
        data = args.training_data
    train_chronologer( output_file_name = model_out_file,
                       training_data = data,
                       test_frac = args.validation_fraction,
                       random_seed = args.seed,
                       start_time = start_time,
                       device = 'cuda' if args.gpu else 'cpu',
                       start_model = args.start_model if args.start_model != 'None' else None,
                       initial_batch_scaler = 1.0 / args.lr_weight, )


if __name__ == "__main__":
    main()
    sys.exit()
