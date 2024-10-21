'''
Molecular masses
'''

from src.custom import return_custom_modifications

masses = { 'proton' : 1.00727646688, 'hydrogen' : 1.007825035, 'carbon' : 12.000000,
           'nitrogen' : 14.003074, 'oxygen' : 15.99491463, 'phosphorus' : 30.973762,
           'sulfur' : 31.9720707, 
           
           'A' : 71.037113805,  'C' : 103.009184505, 'D' : 115.026943065, 'E' : 129.042593135,
           'F' : 147.068413945, 'G' : 57.021463735,  'H' : 137.058911875, 'I' : 113.084064015,
           'K' : 128.094963050, 'L' : 113.084064015, 'M' : 131.040484645, 'N' : 114.042927470,
           'P' : 97.052763875,  'Q' : 128.058577540, 'R' : 156.101111050, 'S' : 87.032028435,
           'T' : 101.047678505, 'V' : 99.068413945,  'W' : 186.079312980, 'Y' : 163.063328575, }

# Modifications
masses['H2O'] = masses['hydrogen']*2 + masses['oxygen']
masses['NH3'] = masses['nitrogen'] + masses['hydrogen']*3
masses['Ox'] = masses['oxygen']
masses['Cam'] = masses['hydrogen']*3 + masses['carbon']*2 + masses['nitrogen'] + masses['oxygen']
masses['Phospho'] = masses['hydrogen'] + masses['oxygen']*3 + masses['phosphorus']
masses['Ac'] = masses['hydrogen']*2 + masses['carbon']*2 + masses['oxygen']
masses['Me'] = masses['carbon'] + masses['hydrogen']*2
masses['Succ'] = masses['carbon']*4 + masses['oxygen']*3 + masses['hydrogen']*5
masses['Ub'] = masses['G']*2
masses['TMT0'] = 224.152478
masses['TMT10'] = 229.162932

# Modified residues
masses['c'] = masses['C'] + masses['Cam']       # Carbamidomethyl cysteine
masses['m'] = masses['M'] + masses['Ox']        # Oxidized methionine
masses['d'] = masses['c'] - masses['NH3']       # S-carbamidomethylcysteine cyclization
masses['e'] = masses['E'] - masses['H2O']       # Pyroglutamate
masses['s'] = masses['S'] + masses['Phospho']   # Phosphoserine
masses['t'] = masses['T'] + masses['Phospho']   # Phosphothreonine
masses['y'] = masses['Y'] + masses['Phospho']   # Phosphotyrosine
masses['a'] = masses['K'] + masses['Ac']        # Acetylated lysine
masses['b'] = masses['K'] + masses['Succ']      # Succinylated lysine
masses['u'] = masses['K'] + masses['Ub']        # Ubiquitinated lysine
masses['n'] = masses['K'] + masses['Me']        # Monomethyl lysine
masses['o'] = masses['K'] + masses['Me']*2      # Dimethyl lysine
masses['p'] = masses['K'] + masses['Me']*3      # Trimethyl lysine
masses['q'] = masses['R'] + masses['Me']        # Monomethyl arginine
masses['r'] = masses['R'] + masses['Me']*2      # Dimethyl arginine
masses['z'] = masses['K'] + masses['TMT0']      # TMT0-modified lysine
masses['x'] = masses['K'] + masses['TMT10']     # TMT10-modified lysine


# Apply custom modifications
for code, residue, mod_mass in return_custom_modifications():
    masses[ code ] = masses[ residue ] + mod_mass
