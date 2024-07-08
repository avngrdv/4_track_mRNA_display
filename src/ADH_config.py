experiment = '4-track_data_analysis_demo'

class constants:
    '''
    Star symbol (*) is internally reserved for stop codons that terminate
    translation.
    
    Plus and underscore symbols (+ and _) are internally reserved tokens.
    Numerals (1234567890) are internally reserved for library design 
    specification. These symbols (123456790+_) should not be used to
    encode amino acids.
    
    Other symbols are OK.
    
    If reading through stop codons is desired (returning sequences containing
    stop codons in the middle), edit the translation table below to encode
    the relevant codons as something other than (*)! The library designs will
    also need to be edited to include new stop codon encodings as possible
    monomers.
    '''    
    translation_table = {
                    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'d',
                    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
                    'TGC':'a', 'TGT':'a', 'TGA':'*', 'TGG':'W',
                   }
    
    #SMILES strings corresponding to all encoded amino acids;
    #sorted alphabetically by amino acid one-letter codes; used
    #to create ECFP feature matrices
    aa_SMILES = [  
                    'N[C@@H](C)C(=O)',           
                    'N[C@@H](CC(=O)O)C(=O)',
                    'N[C@@H](CCC(=O)O)C(=O)',
                    'N[C@@H](Cc1ccccc1)C(=O)',
                    'NCC(=O)',
                    'N[C@@H](Cc1c[nH]cn1)C(=O)',
                    'N[C@@H]([C@H](CC)C)C(=O)',
                    'N[C@@H](CCCCN)C(=O)',        
                    'N[C@@H](CC(C)C)C(=O)',
                    'N[C@@H](CC(=O)N)C(=O)',
                    'O=C[C@@H]1CCCN1',  
                    'N[C@@H](CCC(=O)N)C(=O)',
                    'N[C@@H](CCCNC(=N)N)C(=O)',
                    'N[C@@H](CO)C(=O)',
                    'N[C@@H]([C@H](O)C)C(=O)',
                    'N[C@@H](C(C)C)C(=O)',
                    'N[C@@H](Cc1c[nH]c2c1cccc2)C(=O)',
                    'N[C@@H](Cc1ccc(O)cc1)C(=O)',
                    'CN[C@@H](C)C(=O)', #N-Me-Ala
                    'NC(C=O)=C' #Dha
                ]

    
    #nucleotide complement table;
    #bases are represented by their ascii symbol numbers
    complement_table = {65: 84,
                        67: 71,
                        84: 65, 
                        71: 67, 
                        78: 78
                       }
          
class LibaryDesigns:
      
    #DNA library design parameters     
    dna_templates = [
                        'atgagtgatattacggctgagaacctctacttccagagc112112112112112112112112345112112112112112112112112ggcagctacccatatgacgtgcccgactatgcaggccgatagtgacggggggcggaaa'.upper(),
                        'atgagtgatattacggctgagaacctctacttccagagc112112112112112345112112112112112112112112112112ggcagctacccatatgacgtgcccgactatgcaggccgatagtgacggggggcggaaa'.upper(),
                        'atgagtgatattacggctgagaacctctacttccagagc112112112112112112345112112112112112112112112ggcagctacccatatgacgtgcccgactatgcaggccgatagtgacggggggcggaaa'.upper(),
                        'atgagtgatattacggctgagaacctctacttccagagc112112112112112112112112345112112112112112112ggcagctacccatatgacgtgcccgactatgcaggccgatagtgacggggggcggaaa'.upper(),
                        'atgagtgatattacggctgagaacctctacttccagagc112112112112112112112112112112345112112112112112ggcagctacccatatgacgtgcccgactatgcaggccgatagtgacggggggcggaaa'.upper()
                    ]

    dna_monomers = {
                        1: ('A', 'G', 'T', 'C'),
                        2: ('G', 'T'),
                        3: ('A'),
                        4: ('T'),
                        5: ('G'),
                        
                    }
                    
    #peptide library design parameters
    pep_templates = [
                     'dSDITAENLYFQS11111111211111111GSYPYDVPDYAGR',
                     'dSDITAENLYFQS1111121111111111GSYPYDVPDYAGR',
                     'dSDITAENLYFQS111111211111111GSYPYDVPDYAGR',
                     'dSDITAENLYFQS111111112111111GSYPYDVPDYAGR',
                     'dSDITAENLYFQS1111111111211111GSYPYDVPDYAGR',
                    ]
            
    pep_monomers = {
                        1: ('A', 'a', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'd', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'),
                        2: ('d')
                    }


class TrackerConfig:
    
    #directory holding sequencing data files (fastq or fastq.gz)
    seq_data = '../sequencing_data'
        
    #directory for writing logs to
    logs = '../logs'
    
    #directory that stores fastqparser outputs
    parser_out = '../parser_outputs'
    
    #directory that stores outputs of data analysis ops
    analysis_out = '../parser_outputs'

    #ml_data
    ml_data = '../ml_data'
    
    #ml_data
    model = '../tf_trained_models'    


class LoggerConfig:
       
    #verbose loggers print to the console
    verbose = True
    
    #write logs to file
    log_to_file = True

    #logger level; accepted any of ('DEBUG', INFO', 'WARNING', 'ERROR')
    level = 'INFO'
    
class FastqParserConfig:
    
    #a regex pattern that has to match in order to initiate the ORF
    #used when translation is performed as force_translation=False
    # utr5_seq = 'AGGAGAT......ATG'
    utr5_seq = 'ATG'
 
class PreproConfig:
    
    #TODO: this is not a good idea, have the method take X_val as an arg instead
    import numpy as np
    X_val = np.load('../npy files/h.105_138_pepos_NUM.npy')
    
   
class ClassifierConfig:

    #-----------------------------------
    #model architecture
    #-----------------------------------
    from tf.cnn_model import cnn_vm_J2
    model = cnn_vm_J2
    
    #NN architecture input_dimension
    inp_dim = (18, 204)
    
    #network dropout rate
    drop = 0.15

    #experiment name
    experiment_name = experiment

    #-----------------------------------
    #training metaparameters
    #-----------------------------------
    
    #learning rate
    from tf.schedules import NoamSchedule
    lr = NoamSchedule(inp_dim[-1], warmup_steps=4000)

    #model training optimizer object
    import tensorflow.keras as K
    optimizer = K.optimizers.Adam(
                                  learning_rate=lr,
                                  beta_1=0.9,
                                  beta_2=0.98, 
                                  epsilon=1e-9
                                  )
        
    #training loss function
    # loss = 'binary_crossentropy'
    loss = K.losses.BinaryCrossentropy()
    
    #training metrics
    metrics=[K.losses.MeanSquaredError()]
    
    #max number of epochs; usually smaller if early_stopping callback is set
    epochs = 1000
    
    #fitting process verbosity
    verbosity = 1       
    
    #training batch size
    batch_size = 2048
    
    #tf.Dataset param 
    shuffle_buffer = 10
    
    #callbacks
    from tf.callbacks import EarlyStop, Checkpoint, AdditionalValidationSets
    
    import os
    save_dir = os.path.join(TrackerConfig.model, 'checkpoint_models')
    model_name = experiment + '_model_checkpoint_{epoch:02d}.h5'
    filepath = os.path.join(save_dir, model_name)
        
    callbacks = [
                 EarlyStop(patience=8), 
                 Checkpoint(filepath=filepath, save_best_only=True),
                ]    

    
    
    
    
    
    
    