'''
Peptide manipulation functions
'''

import re, itertools, random
import numpy as np, pandas as pd
from pyteomics.parser import cleave, expasy_rules

import src.constants as constants
from src.masses import masses
from src.tensorize import residues, modseq_to_codedseq



def mod_str_formatter(seq, mod_dict, start_site, end_site, mod_offset=0):
    output_str = ''
    for site in mod_dict:
        if start_site <= site <= end_site: # ensure that mods are not outside substring
            mod_str = str(int(site-mod_offset)) + seq[site-1-mod_offset] +\
                      format(int(mod_dict[site]), '+d') + ' '
            output_str += mod_str
    return output_str[:-1]

def mod_seq_generator(seq, mod_dict):
    mod_seq = ''
    start = 0
    for i in mod_dict:
        mod = '['+format(mod_dict[i], '+f')+']'
        mod_seq = mod_seq + seq[start:i] + mod
        start = i
    mod_seq += seq[start:]
    return mod_seq

def modseq_toModDict(mod_seq):
    mods = {}
    temp_seq = str(mod_seq)
    for i in range(mod_seq.count('[')):
        mod_start_index = temp_seq.find('[')
        mod_end_index = temp_seq.find(']')
        mod_mass = np.float64( temp_seq[mod_start_index+1:mod_end_index] )
        mods[mod_start_index] = mod_mass
        temp_seq = temp_seq[:mod_start_index] + temp_seq[mod_end_index+1:]
    return mods

def modseq_to_seq( mod_seq ):
    return re.sub(r'\[.+?\]','',mod_seq)


def mass_calc( modseq, initial_mass = masses['H2O'], ):
    seq = modseq_to_seq( modseq )
    mods = modseq_toModDict( modseq )
    residue_masses = [ masses[r] for r in seq ]
    for i in mods: 
        residue_masses[ i - 1 ] += mods[ i ]
    return initial_mass + np.sum( residue_masses ) 

def mz_calc( modseq, precursor_charge, ):
    mass = mass_calc( modseq )
    return ( mass + masses['proton']*precursor_charge ) / precursor_charge

def ladder_mz_generator(mod_seq, charge=2, max_fragment_z=3, phos_loss=False, ):
    seq = modseq_to_seq( mod_seq )
    ## Identify modifications
    mods = modseq_toModDict(mod_seq)

    residue_masses = np.array( [ masses[a] for a in seq ] )
    mod_mass_shifts = np.array( [ mods[i] if i in mods else 0.0 
                             for i in range(1,len(seq)+1) ] )
    residue_wMod_shifts = residue_masses + mod_mass_shifts
    
    peptide_mass = mass_calc( mod_seq )
    frag_masses = {}
    frag_masses['b'] = np.cumsum( residue_wMod_shifts[:-1] )
    frag_masses['y'] = peptide_mass - frag_masses['b']
    
    if phos_loss:
        phos_loc = np.array( [ m == masses['Phospho'] for m in mod_mass_shifts ] )
        phos_mask = np.cumsum( phos_loc ) >= 1
        b_phos_mask = phos_mask[:-1]
        frag_masses['b-PO4'] = b_phos_mask *\
                            (frag_masses['b'] - b_phos_mask*(masses['Phospho']+masses['H2O']))
        y_phos_mask = phos_mask[::-1][:-1]
        frag_masses['y-PO4'] = y_phos_mask *\
                            (frag_masses['y'] - y_phos_mask*(masses['Phospho']+masses['H2O']))
                           
    mz_array = []
    for frag_z in range( 1, max_fragment_z+1, ):
        for frag_type in frag_masses:
            mzs = (frag_masses[frag_type]+masses['proton']*frag_z) / frag_z
            mzs = mzs * ( frag_z <= charge ) * ( frag_masses[frag_type] > 0 ) # Mask out impossible ions
            mz_array.append( mzs )
    return np.array( mz_array ) # Current structure is (12, peplen-1)), consider flattening



def peptide_generator( seq, protease, miscleavages, min_peptide_len, max_peptide_len, ):
    protease_to_expasy = { 'trypsin' : 'trypsin',
                           'chymotrypsin' : 'chymotrypsin high specificity',
                           'aspn' : 'asp-n',
                           'gluc' : 'glutamyl endopeptidase', }
    peptides = cleave( seq.upper(), 
                       expasy_rules[ protease_to_expasy[ protease ] ], 
                       missed_cleavages = miscleavages, )
    # Limit peptides based on min/max length and only canonical amino acids
    peptides = set( [ p for p in peptides 
                      if len(p) >= min_peptide_len and 
                         len(p) <= max_peptide_len and
                         set(residues[:20]) | set(p) == set(residues[:20]) ] )
    return sorted( peptides )

def protein_to_peptides( accession, sequence, protease, miscleavages, min_peptide_len, max_peptide_len, ):
    peptides = peptide_generator( sequence, protease, miscleavages, min_peptide_len, max_peptide_len, )
    df = pd.DataFrame( { 'PeptideSeq' : peptides,
                         'isDecoy' : False,
                         'ProteinAccession' : accession, } )
    return df

def peptide_df_to_decoys( df, reverse = True, ):
    df[ 'isDecoy' ] = False # CURRENTLY FIX TO ENTRAPMENT MODE
    df[ 'ProteinAccession' ] = [ 'ENTRAP_' + p for p in df.ProteinAccession ]
    df[ 'PeptideSeq' ] = [ peptide_to_decoy(p, reverse=reverse,) for p in df.PeptideSeq ]
    return df.dropna()
                


def filter_modified_miscleavages( modified_peptides, protease, n_miscleavages, ):
    # Need to remove anything with a terminal modification
    # or peptides that have multiple real miscleavages w/o mods
    protease_rules = { 'trypsin' :      { 'terminus' : 'C', 
                                          'cleavage_residues' : [ 'R', 'K' ], 
                                          'exclude_cleavage' : [ 'RP', 'KP' ] },
                       'chymotrypsin' : { 'terminus' : 'C', 
                                          'cleavage_residues' : [ 'F', 'Y', 'W', ],
                                          'exclude_cleavage' : [ 'FP', 'YP', 'WP', 'WM', ] },
                       'gluc' :         { 'terminus' : 'C', 
                                          'cleavage_residues' : [ 'E' ],
                                          'exclude_cleavage' : [ ], },
                       'aspn' :         { 'terminus' : 'N', 
                                          'cleavage_residues' : [ 'D' ],
                                          'exclude_cleavage' : [ ], } }
    rules = protease_rules[protease]
    
    filtered_peptides = []
    for mod_peptide in modified_peptides:
        coded_peptide = modseq_to_codedseq( mod_peptide )
        if coded_peptide: # Compatible with Cartographer/Chronologer
            peptide = modseq_to_seq( mod_peptide )
            terminal_residue = peptide[-1] if rules['terminus'] == 'C' else peptide[0]
            enzymatic_cut = terminal_residue in rules['cleavage_residues'] # Assume protein N/C-term if false
            coded_terminus = coded_peptide[-2] if rules['terminus'] == 'C' else coded_peptide[1]
            cleavage_at_mod = coded_terminus != terminal_residue and enzymatic_cut
            if not cleavage_at_mod:
                n_miscleaves = np.sum( [ coded_peptide.count(r) for r in rules['cleavage_residues'] ] ) -\
                               np.sum( [ coded_peptide.count(r) for r in rules['exclude_cleavage'] ] ) -\
                               enzymatic_cut
                if n_miscleaves <= n_miscleavages:
                    filtered_peptides.append( mod_peptide )
    return filtered_peptides

def reverse_peptide( peptide ):
    # Reverse the peptide sequence but hold the termini constant
    return peptide[0]+peptide[1:-1][::-1]+peptide[-1]


def scramble_peptide( peptide ):
    # Shuffle all but the terminal residues
    # Ensure that it doesn't match either the regular or reverse peptides 
    rev_peptide = reverse_peptide( peptide )
    if rev_peptide == peptide:
        return np.nan
    else:
        scrambled_peptide = peptide
        while scrambled_peptide in [ peptide, rev_peptide ]:
            #print(scrambled_peptide)
            scrambled_peptide = peptide[0] +\
                                ''.join( random.sample( peptide[1:-1], 
                                                        len(peptide)-2 ) ) +\
                                peptide[-1]
            #print(scrambled_peptide)
            #print()
        return scrambled_peptide
        


def peptide_to_decoy( peptide, reverse = True, ):
    if reverse:
        return reverse_peptide( peptide )
    else:
        return scramble_peptide( peptide )
        
        



def apply_mods( peptides, fixed_mods, var_mods, max_n_var_mods, ):
    modified_peptides = []
    for peptide in peptides:
        fixed_mod_list = [ ( i+1, fixed_mods[r] ) 
                           for i, r in enumerate(peptide) 
                           if r in fixed_mods ]
        var_mod_list = [ ( i+1, var_mods[r] ) 
                         for i, r in enumerate(peptide) 
                         if r in var_mods ]
        var_mod_lists = [ list(s) for n in range(1, max_n_var_mods+1) 
                          for s in itertools.combinations( var_mod_list, n ) ]
        for var_mod_list in [[]] + var_mod_lists:
            mod_dict = dict( sorted( fixed_mod_list + var_mod_list, key=lambda x:x[0], ) )
            modified_peptides.append( mod_seq_generator( peptide, mod_dict, ) )
    return modified_peptides

