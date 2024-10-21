import os
import torch

from src.chronologer.model import initialize_chronologer_model

def serialize_chronologer_model( model_file=None, ):
    # Will operate on latest timestamped model by default
    if not model_file:
        model_dir = os.path.join( '..', 'models' )
        
        models = [ m for m in os.listdir( model_dir ) 
                   if m[:11] == 'Chronologer' and m[-3:] == '.pt' ]
        models.sort( reverse=True, )
        model_file = os.path.join( model_dir, models[0] )
    # Output will be same as file name but with spt extension
    output_file = model_file[ : model_file.rfind('.') ] + '.spt'
    model = initialize_chronologer_model( model_file, )
    serialized_model = torch.jit.script( model )
    serialized_model.save( output_file )
    return None
    
