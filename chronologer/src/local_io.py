import os, sqlite3
import numpy as np, pandas as pd
from pyteomics.parser import cleave, expasy_rules

import src.constants as constants
from src.masses import masses
from src.tensorize import modseq_to_codedseq



def patch_modseqs_with_cyclo_N_term( entries ):
    # The current version of EncyclopeDIA (1.12) currently does not provide monoisotopic masses for
    # cyclization of glutamate or carbamidocysteine, so need to include a patch to properly annotate
    # these ions
    modseqs = list( entries.chronologer_sequence )
    modseqs = [ p.replace( 'E[-18.0]', 'E[-' + str(masses['H2O']) + ']' ) for p in modseqs ]
    modseqs = [ p.replace( 'C[-17.0]', 'C[+' + str(masses['Cam'] - masses['NH3']) + ']' ) for p in modseqs ]
    entries['PeptideModSeq'] = modseqs
    return 0

def patch_modseqs_with_N_acetyl( entries ):
    # EncyclopeDIA reports N-acetylation as a monoisotopic mass to the right of the first amino acid,
    # which Cartographer/Chronologer will read as side chain modifications rather than N-terminal,
    # and so it needs to be shifted to the first position to be properly parsed.
    modseqs = list( entries.chronologer_sequence )
    for i, p in enumerate( modseqs ):
        if p[1:13] == '[+42.010565]':
            modseqs[i] = '[+42.010565]' + p[0] + p[13:]
    entries['PeptideModSeq'] = modseqs
    return 0


def read_rt_database( db_loc, chronologer_compatible_only=True, ):
    delimiter = '\t'
    if db_loc.split('.')[-1] == "csv":
        delimiter = ','
    df = pd.read_csv( db_loc, delimiter=delimiter, )
    patch_modseqs_with_cyclo_N_term( df )
    patch_modseqs_with_N_acetyl( df )
    coded_seqs = [ modseq_to_codedseq( p ) for p in df.chronologer_sequence]
    df['CodedPeptideSeq'] = [ c_seq if c_seq else np.nan for c_seq in coded_seqs ]
    if chronologer_compatible_only:
        #df = df.dropna() # Drop any entries that do not produce a compatible coded seq
        for (i, p) in enumerate(df.CodedPeptideSeq):
            if p!=p:
                print(p, df.sequence[i], df.PeptideModSeq[i], df.mods[i])
        df['PeptideLength'] = [ len(p)-2 for p in df.CodedPeptideSeq ]
        df = df[ ( df.PeptideLength.between( constants.min_peptide_len, 
                                             constants.max_peptide_len, ) ) &
                ( pd.Series([ p.count('[') == 0 for p in df.CodedPeptideSeq ]) ) ]
    return df


def read_table( sqlite_file, table_name ):
    con = sqlite3.connect( sqlite_file )
    c = con.cursor()
    df = pd.read_sql_query("SELECT * from " + table_name, con)
    con.close()
    if table_name == 'entries':
        patch_modseqs_with_cyclo_N_term( df )
        patch_modseqs_with_N_acetyl( df )
    return df


def read_fasta( file_name ): # Return dict with keys = names, values = seqs
    seq_dict = {}
    trigger = 0
    for line in open( file_name, 'r' ):
        if line[0] == '>': # New seq
            if trigger != 0:
                seq_dict[ name ] = seq
            name = line[1:].rstrip() #.split(' ')[0] # PATCH NAMES FOR UNIPROT
            seq = ''
            trigger = 1
        else:
            seq += line.rstrip()
    seq_dict[name] = seq # Add the final sequence since 
    return seq_dict
    

    


dlib_dtypes= { 'metadata' :           { 'Key' :                                'string',
                                        'Value' :                              'string',
                                      },
                      
               'entries' :            { 'PrecursorMz' :                        'double', 
                                        'PrecursorCharge' :                    'int', 
                                        'PeptideModSeq' :                      'string',
                                        'PeptideSeq' :                         'string',
                                        'Copies' :                             'int',
                                        'RTInSeconds' :                        'double',
                                        'Score':                               'double',
                                        'MassEncodedLength' :                  'int',
                                        'MassArray' :                          'blob',
                                        'IntensityEncodedLength' :             'int',
                                        'IntensityArray' :                     'blob',
                                        'CorrelationEncodedLength' :           'int',
                                        'CorrelationArray' :                   'blob',
                                        'RTInSecondsStart' :                   'double',
                                        'RTInSecondsStop' :                    'double',
                                        'MedianChromatogramEncodedLength' :    'int',
                                        'MedianChromatogramArray' :            'blob',
                                        'SourceFile' :                         'string', 
                                       },
                       
              'peptidetoprotein' :     { 'PeptideSeq' :                        'string',
                                         'isDecoy' :                           'boolean',
                                         'ProteinAccession' :                  'string', 
                                       },

              'peptidequants' :        { 'PrecursorCharge' :                   'int',
                                         'PeptideModSeq' :                     'string',
                                         'PeptideSeq' :                        'string',
                                         'SourceFile' :                        'string',
                                         'RTInSecondsCenter' :                 'double',
                                         'RTInSecondsStart' :                  'double',
                                         'RTInSecondsStop' :                   'double',
                                         'TotalIntensity' :                    'double',
                                         'NumberOfQuantIons' :                 'int',
                                         'QuantIonMassLength' :                'int',
                                         'QuantIonMassArray' :                 'blob',
                                         'QuantIonIntensityLength' :           'int',
                                         'QuantIonIntensityArray' :            'blob',
                                         'BestFragmentCorrelation' :           'double',
                                         'BestFragmentDeltaMassPPM' :          'double',
                                         'MedianChromatogramEncodedLength' :   'int',
                                         'MedianChromatogramArray' :           'blob',
                                         'MedianChromatogramRTEncodedLength' : 'int',
                                         'MedianChromatogramRTArray' :         'blob',
                                         'IdentifiedTICRatio' :                'double',
                                       },

              'peptidelocalizations' : { 'PrecursorCharge' :                   'int',
                                         'PeptideModSeq' :                     'string',
                                         'PeptideSeq' :                        'string',
                                         'SourceFile' :                        'string',
                                         'LocalizationPeptideModSeq' :         'string',
                                         'LocalizationScore' :                 'double',
                                         'LocalizationFDR' :                   'double',
                                         'LocalizationIons' :                  'string',
                                         'NumberOfMods' :                      'int',
                                         'NumberOfModifiableResidues' :        'int',
                                         'IsSiteSpecific' :                    'boolean',
                                         'IsLocalized' :                       'boolean',
                                         'RTInSecondsCenter' :                 'double',
                                         'LocalizedIntensity' :                'double',
                                         'TotalIntensity' :                    'double',
                                       },

              'peptidescores' :        { 'PrecursorCharge' :                   'int',
                                         'PeptideModSeq' :                     'string',
                                         'PeptideSeq' :                        'string',
                                         'SourceFile' :                        'string',
                                         'QValue' :                            'double',
                                         'PosteriorErrorProbability' :         'double',
                                         'IsDecoy' :                           'boolean',
                                       },
                                    
              'proteinscores' :        { 'ProteinGroup' :                      'int',
                                         'ProteinAccession' :                  'string',
                                         'SourceFile' :                        'string',
                                         'QValue' :                            'double',
                                         'MinimumPeptidePEP' :                 'double',
                                         'IsDecoy' :                           'boolean',
                                       },

              'retentiontimes' :       { 'SourceFile' :                        'string',
                                         'Library' :                           'float',
                                         'Actual' :                            'float',
                                         'Predicted' :                         'float',
                                         'Delta' :                             'float',
                                         'Probability' :                       'float',
                                         'Decoy' :                             'boolean',
                                         'PeptideModSeq' :                     'string',
                                       },
             }

def create_dlib( file_name, overwrite=True, ):
    if os.path.isfile( file_name ):
        if overwrite: 
            os.remove( file_name )
        else:
            assert False, 'DLIB already exists and overwrite set to False'
    con = sqlite3.connect( file_name )
    cursor = con.cursor()
    for table in dlib_dtypes:
        creation_string = 'CREATE TABLE ' + table + '( '
        for column in dlib_dtypes[ table ]:
            dtype_str = column + ' ' + dlib_dtypes[table][column] + ', '
            creation_string += dtype_str
        creation_string = creation_string[:-2] + ' )'
        cursor.execute( creation_string )
    cursor.execute( "INSERT INTO metadata (Key, Value) VALUES ('version', '0.1.14')" )
    con.commit()
    con.close()
    return 0

    
    
def append_table_to_dlib( table, table_name, dlib_file ):
    # Add any missing columns
    #for column in dlib_dtypes[ table_name ]:
    #    if column not in table.columns:
    #        table[ column ] = None
    
    # Append to existing dlib
    con = sqlite3.connect(dlib_file)
    cursor = con.cursor()
    table.to_sql(table_name, con, if_exists='append', index=False, )
    con.commit()
    con.close()
    return 0

def overwrite_table_to_dlib( table, table_name, dlib_file ):
    con = sqlite3.connect(dlib_file)
    cursor = con.cursor()
    table.to_sql(table_name, con, if_exists='replace', index=False, )
    con.commit()
    con.close()
    return 0

def table_exists_check( cursor, table_name ):
    cursor.execute( "SELECT count(name) FROM sqlite_master WHERE type='table' AND name='" + 
                    table_name +"'" )
    if cursor.fetchone()[0]==1:
        return True
    else:
        return False

def build_dlib_index( dlib_file, ):
    con = sqlite3.connect(dlib_file)
    cursor = con.cursor()
    if table_exists_check( cursor, 'metadata', ):
        cursor.execute( "create index if not exists 'Key_Metadata_index' on 'metadata' ('Key' ASC)" )
    if table_exists_check( cursor, 'entries', ):
        cursor.execute( "create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Entries_index' on 'entries' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)" )
        cursor.execute( "create index if not exists 'PeptideSeq_Entries_index' on 'entries' ('PeptideSeq' ASC)" )
        cursor.execute( "create index if not exists 'PrecursorMz_Entries_index' on 'entries' ('PrecursorMz' ASC)" )
    if table_exists_check( cursor, 'peptidequants', ):
        cursor.execute( "create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Peptides_index' on 'peptidequants' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)" )
        cursor.execute( "create index if not exists 'PeptideSeq_Peptides_index' on 'peptidequants' ('PeptideSeq' ASC)" )
    if table_exists_check( cursor, 'peptidelocalizations', ):
        cursor.execute( "create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Localizations_index' on 'peptidelocalizations' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)" )
        cursor.execute( "create index if not exists 'PeptideSeq_Localizations_index' on 'peptidelocalizations' ('PeptideSeq' ASC)" )
    if table_exists_check( cursor, 'peptidescores', ):
        cursor.execute( "create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Scores_index' on 'peptidescores' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)" )
        cursor.execute( "create index if not exists 'PeptideSeq_Scores_index' on 'peptidescores' ('PeptideSeq' ASC)" )
    if table_exists_check( cursor, 'proteinscores', ):
        cursor.execute( "create index if not exists 'ProteinGroup_ProteinScores_index' on 'proteinscores' ('ProteinGroup' ASC)" )
        cursor.execute( "create index if not exists 'ProteinAccession_ProteinScores_index' on 'proteinscores' ('ProteinAccession' ASC)" )
    if table_exists_check( cursor, 'fragmentquants', ):
        cursor.execute( "create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Fragments_index' on 'fragmentquants' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)" )
        cursor.execute( "create index if not exists 'PeptideSeq_Fragments_index' on 'fragmentquants' ('PeptideSeq' ASC)" )
    if table_exists_check( cursor, 'peptidetoprotein', ):
        cursor.execute( "create index if not exists 'ProteinAccession_PeptideToProtein_index' on 'peptidetoprotein' ('ProteinAccession' ASC)" )
        cursor.execute( "create index if not exists 'PeptideSeq_PeptideToProtein_index' on 'peptidetoprotein' ('PeptideSeq' ASC)" )
    con.commit()
    con.close()
    return 0