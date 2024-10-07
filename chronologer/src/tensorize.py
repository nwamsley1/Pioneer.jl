
import re
import numpy as np

import torch
from torch.utils.data import TensorDataset

import src.constants as constants
from src.custom import return_custom_modifications

residues = [ # Primary 20 amino acids            
             'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
             'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
             
             # Modified amino acids, see masses.py for details
             'c', 'm', 'd', 'e', 's', 't', 'y', 'a', 'b', 'u', 
             'n', 'o', 'p', 'q', 'r', 'x', 'z',
             
             # N/C-terminal states
             # '-' = free N-term, '_' = free C-term
             # '^' = N-acetyl, 
             # '(' = S-carbamidomethylcysteine cyclization
             # ')' = pyroglutamate, 
             # '&' = N-terminal TMT0, '*' = N-terminal TMT10
             '-', '^', '(', ')', '&', '*', '_', 
             
             # Include 10 user assignable values
             '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
          ]

n_states = len(residues)+1

aa_to_int = dict(zip(residues,range(1,len(residues)+1)))


mod_regex_keys = { r'M\[\+15\.99.{,6}\]':'m', r'C\[\+57\.02.{,6}\]':'c', 
                   r'C\[\+39\.99.{,6}\]':'d', r'E\[\-18\.01.{,6}\]':'e', r'Q\[\-17\.02.{,6}\]':'e',
                   r'S\[\+79\.96.{,6}\]':'s', r'T\[\+79\.96.{,6}\]':'t', r'Y\[\+79\.96.{,6}\]':'y',
                   r'K\[\+42\.01.{,6}\]':'a', r'K\[\+100\.0.{,6}\]':'b', r'K\[\+114\.0.{,6}\]':'u', 
                   r'K\[\+14\.01.{,6}\]':'n', r'K\[\+28\.03.{,6}\]':'o', r'K\[\+42\.04.{,6}\]':'p',
                   r'R\[\+14\.01.{,6}\]':'q', r'R\[\+28\.03.{,6}\]':'r', 
                   r'K\[\+224\.1.{,6}\]':'z', r'K\[\+229\.1.{,6}\]':'x', }

# Apply custom modifications
for code, residue, mod_mass in return_custom_modifications():
    sign = '\+' if mod_mass >= 0 else '-'
    mod_str = residue + '\[' + sign + str(abs(mod_mass)).replace('.','\.') + '\]'
    mod_regex_keys[ mod_str ] = code
    

nterm_keys = { '+42.01' : '^', '+224.1' : '&', '+229.1' : '*' }
## pyroglu = (e, cyclocys = )d

def modseq_to_codedseq( seq ):
    mod_seq = seq
    for mod_id in mod_regex_keys:
        seq = re.sub( mod_id, mod_regex_keys[mod_id], seq )
        
    # N/C mods
    if seq[0] == 'd' : seq = ')' + seq
    elif seq[0] == 'e' : seq = '(' + seq
    elif seq[0] == '[' : seq = nterm_keys[ seq[1:7] ] + seq[ seq.find(']')+1: ]
    else: seq = '-' + seq
    seq = seq + '_'  
    
    # Ensure that there are no additional modifications
    return seq if seq.count('[') == 0 else None

    #if seq.count('[') > 0: #or\
    #    print( mod_seq, seq )
    #   #len(seq)-2 < constants.min_peptide_len or\
    #   #len(seq)-2 > constants.max_peptide_len:
    #    return None
    #else:
    #    return seq


def codedseq_to_array(seq, max_size=constants.max_peptide_len+2):
    seq_by_int = [ aa_to_int[ seq[i] ] for i in range( len( seq ) ) ]
    seq_by_int += [0]*(max_size - len(seq_by_int))
    return np.asarray( seq_by_int, 'int64' )

def modseqs_to_tensor( mod_seqs, ):
    coded_peptides = [ modseq_to_codedseq( p ) for p in mod_seqs ]
    peptide_array = np.array( [ codedseq_to_array( p ) for p in coded_peptides ] )
    peptide_tensor = torch.Tensor( peptide_array ).to( torch.int64 )
    return peptide_tensor


def hi_db_to_tensors( hi_db ):
    seq_array = np.asarray( [ codedseq_to_array( p ) for p in hi_db.CodedPeptideSeq ],
                            'int64', )
    hi_array = np.asarray( [ [ hi ] for hi in hi_db.HI ], 'float32', )
    sources = sorted( set( hi_db.Source ) )
    source_array = np.asarray( [ [ float( s == sx ) for sx in sources ] for s in hi_db.Source ],
                               'float32', )
    tensors = [ torch.Tensor( x ) for x in [ seq_array, hi_array, source_array] ]
    tensors[0] = tensors[0].to(torch.int64) # Need to ensure seq tensor are ints for embedding layer
    
    return TensorDataset( *tensors )


def precursor_charge_to_array( precursor_charge, copies=0 ):
    array = np.zeros( constants.n_charges, )
    array[ precursor_charge - constants.min_precursor_charge ] = 1
    if copies == 0: # Return a 1d array
        return array
    else: # Expand to 2D array
        return np.tile( array, (copies,1), )
        

def pad_array( array, n_ion_types=4, ): #peplen, precursor_z, ):
    array = np.array( array ) # Ensure its an array coming in
    pad_size = constants.max_peptide_len - 1 - array.shape[1]
    padded_array = np.pad( array, ( (0,0), (0,pad_size) ), constant_values=-1, )
    
    #padded_array = np.full( ( n_ion_types, constants.max_peptide_len-1 ), -1.0, 'float32' )
    #for i, row in enumerate( array ):
    #    padded_array[ i, :len(row)] = row
    return padded_array 

def spectra_to_tensors( spectra_dicts ): # Input is list of spectral dictionaries
    peptide_tensor = torch.Tensor( np.array( [ codedseq_to_array( d['coded_seq'] ) 
                                               for d in spectra_dicts ] ) ).to( torch.int64 )
    charge_tensor = torch.Tensor( np.array( [ precursor_charge_to_array( d['precursor_charge'] )
                                             for d in spectra_dicts ] ) ).to( torch.float32 )
    intensity_tensor = torch.tensor( np.array( [ pad_array( d['ladder_int'] )
                                                 for d in spectra_dicts ] ) ).to( torch.float32 )
    return peptide_tensor, charge_tensor, intensity_tensor


#def precursor_charges_to_array( precursor_charge, batch_size, ):
#    charge_ohe = [ precursor_charge == z for z in range( min_precursor_charge, 
#                                                         max_precursor_charge+1, ) ]
#    return np.array( [ charge_ohe ] * batch_size )
    