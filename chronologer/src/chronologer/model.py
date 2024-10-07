
import torch
import torch.nn as nn

import src.constants as constants
from src.tensorize import residues
from src.core_layers import resnet_block
from src.chronologer.settings import hyperparameters, training_parameters


class chronologer_model( nn.Module ):
    def __init__( self, vec_length, n_states, embed_dim, n_blocks, kernel, drop_rate, act_fx, ):
        super().__init__()
        self.seq_embed = nn.Embedding( n_states, embed_dim, padding_idx=0, )
        torch.nn.init.kaiming_normal_( self.seq_embed.weight, nonlinearity='linear', )
        self.resnet_blocks = nn.Sequential( *[ resnet_block( embed_dim, 
                                                             embed_dim,
                                                             kernel, 
                                                             d+1, 
                                                             act_fx, ) for d in range(n_blocks) ] )
        self.dropout = nn.Dropout( drop_rate, )
        self.flatten = nn.Flatten()
        self.output = nn.Linear( vec_length * embed_dim, 1, )
        nn.init.xavier_normal_( self.output.weight, )
        nn.init.constant_( self.output.bias.data, 0.0, )
        
    def forward( self, x, ):
        x = self.seq_embed( x ).transpose( 1, -1, )
        x = self.resnet_blocks( x )
        x = self.dropout( x )
        x = self.flatten( x )
        return self.output( x )
        
        
def initialize_chronologer_model( model_file = None, gpu = False, ):
    
    device = 'cuda' if gpu else 'cpu'
    
    model = chronologer_model( constants.max_peptide_len + 2,
                               len( residues ) + 1,
                               hyperparameters[ 'embed_dimension' ],
                               hyperparameters[ 'n_resnet_blocks' ],
                               hyperparameters[ 'kernel_size' ],
                               training_parameters[ 'dropout_rate'],
                               hyperparameters[ 'activation_function' ], )
    if model_file:
        model.load_state_dict( torch.load( model_file, 
                                           map_location=torch.device(device), ), 
                               strict=True, )
    model.eval()
    return model