'''
Code for loading custom amino acid modifications 
that are supplied in user_modifications.txt
'''

import os

def return_custom_modifications():
    src_dir = os.path.dirname(os.path.abspath(__file__))
    custom_modifications = [ ]
    for line in open( os.path.join( src_dir, '..', 'user_modifications.txt', ) ):
        if line[0] != '#':
            parts = line.rstrip().split('\t')
            assert len(parts) == 3, 'Invalid modification: '+line.rstrip()
            code, residue, mod_str = parts
            custom_modifications.append( ( code, residue, float(mod_str) ) )
    return custom_modifications
            
