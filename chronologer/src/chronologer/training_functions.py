
import os, sys
from datetime import datetime
import numpy as np
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#from torch.utils.data import ConcatDataset
import pandas as pd


import chronologer.src.constants as constants 
from chronologer.src.local_io import read_rt_database
from chronologer.src.chronologer.settings import training_parameters
from chronologer.src.tensorize import hi_db_to_tensors
from chronologer.src.chronologer.model import initialize_chronologer_model
from chronologer.src.loss_functions import RT_masked_negLogL, LogL_Loss
from chronologer.src.training_loop import train_model
from chronologer.src.utils import timer


def rt_to_tensor_database( data_files, test_frac, random_seed, ):
    dbs = { 'train' : [ ], 'test' : [ ], }
    #n_sources = 0
    sources = set()
    for data_file in data_files:
        db = read_rt_database( data_file, )
        for source, source_db in db.groupby('Source'):
            split_idx = int( np.round( test_frac * len(source_db) ) )
            if split_idx >= 2:
                source_db = source_db.sample( frac=1, random_state=random_seed, )
                source_db = source_db.reset_index() # Reset now that its shuffled
                dbs['train'].append( source_db.iloc[ split_idx:, ] )
                dbs['test'].append( source_db.iloc[ :split_idx, ] )
                sources.add( source )
            else:
                print( 'Insufficient observations of ' + source + ' for training' )
    n_sources = len( sources )
    data = dict( [ ( t, hi_db_to_tensors( pd.concat( dbs[t] ) ) ) for t in dbs ] )
    return data, n_sources
    

def train_chronologer( output_file_name, training_data, test_frac, random_seed, start_time, 
                       device = None, start_model = None, initial_batch_scaler = 1.0 ):
    # Prepare data
    print( 'Chronologer training initiated (' + timer( start_time ) + ')' )
    datasets, n_sources = rt_to_tensor_database( training_data, test_frac, random_seed ) #merge_datasets( datasets )
    print( 'Training data successfully loaded (' + timer( start_time ) + ')' )
    model = initialize_chronologer_model( model_file = start_model, )
    if start_model:
        print( 'Previous Chronologer model loaded: '+ start_model )
    
    
    loss_fx = LogL_Loss( n_sources=n_sources ) #RT_masked_negLogL( n_sources, )
    
    parameters = list(model.parameters()) + list(loss_fx.parameters())
    optimizer = training_parameters[ 'optimizer' ]( parameters,
                                                    lr = training_parameters[ 'learning_rate' ], )
    
    if device:
        train_device, eval_device = [device]*2
    else:
        train_device = training_parameters[ 'train_device' ]
        eval_device = training_parameters[ 'eval_device' ]
    
    print( 'Ready to begin Chronologer training (' + timer( start_time ) + ')' )
    final_loss, _ = train_model( model, 
                                 datasets, 
                                 training_parameters[ 'initial_batch_size' ] * initial_batch_scaler, 
                                 training_parameters[ 'max_batch_size' ], 
                                 training_parameters[ 'epochs_to_2x_batch' ],
                                 loss_fx, 
                                 optimizer,
                                 training_parameters[ 'n_epochs'], 
                                 train_device, 
                                 eval_device,
                                 output_file_name,
                                 start_time, )
    return final_loss
    