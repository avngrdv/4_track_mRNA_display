# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 22:42:26 2022
@author: Alex Vinogradov
"""

import os, inspect
import numpy as np
import clibas.plotters as Plotter
from clibas.baseclasses import Handler


class DataAnalysisTools(Handler):
    '''
    Crafted wrappers around various data analysis tools to simplify their
    calling during the pipeline creation process.
    '''
    
    def __init__(self, *args):
        super(DataAnalysisTools, self).__init__(*args)
    
        self._validate_designs()
        self._validate_constants()
        return

    def __repr__(self):
        return '<DataAnalysisTools object>'
    
    def length_analysis(self, where=None, save_txt=False):
        
        self._where_check(where)
        la = LengthAnalysis(self.__dict__)
        op = la.len_summary(where=where, save_txt=save_txt)
        return op
    
    def q_score_analysis(self, loc=None, save_txt=False):
        
        if loc is not None:
            self._where_check('dna')
            self._loc_check(loc, self.D_design)
            
        QA = QScoreAnalysis(self.__dict__)
        op = QA.q_summary(loc=loc, save_txt=save_txt)
        return op
    
    def sequence_convergence_analysis(self, where=None, alphabet=None):

        self._where_check(where)
        alphabet = self._infer_alphabet(where, alphabet)
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.sequence_level_convergence(where=where, alphabet=alphabet)
        return op
        
    def token_convergence_analysis(self, where=None, 
                                   loc=None, 
                                   alphabet=None,
                                   save_txt=False):
    
        self._where_check(where)
        design = self._infer_design(where)
        alphabet = self._infer_alphabet(where, alphabet)
        if loc is not None:
            self._loc_check(loc, design)
        
        CA = ConvergenceAnalys(self.__dict__)
        op = CA.token_level_convergence(where=where, 
                                        loc=loc, 
                                        alphabet=alphabet,
                                        save_txt=save_txt
                                       )
        return op
            
class LengthAnalysis(Handler):
    
    def __init__(self, *args):
        super(LengthAnalysis, self).__init__(*args)

    def __repr__(self):
        return '<LengthAnalysis object>'    
    
    def len_summary(self, where=None, save_txt=False):
        '''
        For each sample in Data, compute the distribution of peptide/DNA sequence lengths
        (specified by 'where') and plot the resulting histogram in the parser output folder
        as specified by config.py. Optionally, the data can also be written to a txt file.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          						  
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots				
    						  
        Returns:
                Data object (no transformation)
        '''
        
        def length_summary(data):
            
            recast_data = self._recast_data(data, where=where)
            self._prepare_destinations(recast_data, self.dirs.analysis_out)        
            for sample in recast_data:
                
                arr = sample.X
                L = self._L_summary(arr)            
                L, counts = np.unique(L, return_counts=True)
                
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                fname = f'{sample.name}_{where}_L_distribution'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.L_distribution(L, counts, 
                                                      where=where, 
                                                      basename=basename
                                                     )
                
                if save_txt:
                    np.savetxt(basename + '.csv',
                               np.array((L, counts)).T,
                               delimiter=',', 
                               header='Seq length,Count')
            
            return data
        return length_summary    
    
class ConvergenceAnalys(Handler):
    
    def __init__(self, *args):
        super(ConvergenceAnalys, self).__init__(*args)
        
        from clibas.misc import (
                                 shannon_entropy,
                                 get_freqs,
                                 positional_conservation
                                )

        self.shannon_entropy = shannon_entropy
        self.get_freqs = get_freqs
        self.positional_conservation = positional_conservation
        return

    def __repr__(self):
        return '<ConvergenceAnalys object>'

    def sequence_level_convergence(self, where=None, alphabet=None):
        '''
        For each sample in Data, perform basic library convergence analysis on 
        a sequence level. Computes normalized Shannon entropy and postition-wise 
        sequence conservation. Plots the results in the parser output folder as 
        specified by config.py.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.
                          
                alphabet: token alphabet for the datasets to be analyzed.
                          will be automatically inferred if 'where' is 
                          specified. otherwise, dtype: any of (list, tuple, ndarray)
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        
        def sequence_level_convergence_summary(data):
    
            recast_data = self._recast_data(data, where)
            self._prepare_destinations(recast_data, self.dirs.analysis_out)
            
            for sample in recast_data:
    
                self._transform_check(sample, inspect.stack()[0][3])                
                arr = sample.X
                
                shannon, counts = self.shannon_entropy(arr, norm=True)
                freq = self.get_freqs(arr, alphabet)
                seq_conservation = self.positional_conservation(freq)
                
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                fname = f'{sample.name}_{where}_library_convergence'
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.dataset_convergence(counts, shannon, where, basename)                
                
                fname = f'{sample.name}_{where}_sequence_conservation'
                basename = os.path.join(destination, fname)                
                Plotter.SequencingData.conservation(seq_conservation, where, basename)
                
            return data
        return sequence_level_convergence_summary


    def token_level_convergence(self, where=None, loc=None, alphabet=None, save_txt=False):
        '''
        Perform basic library convergence analysis at a token level. For each sample in Data, 
        computes the frequency of each token in the dataset. Plots the results in the parser 
        output folder as specified by config.py. Optionally, the data can also be written to
        a txt file.
    
        Parameters:
                   where: 'dna' or 'pep' to specify which dataset the op 
                          should work on.

                alphabet: token alphabet for the datasets to be analyzed.
                          will be automatically inferred if 'where' is 
                          specified. otherwise, dtype: any of (list, tuple, ndarray)
    						  
                     loc: if not None: a list of ints to specify regions to be 
                          analyzed; in this case, the op will collapse sample's
                          internal state (see explanation for Data objects)
                          
                          if left None: get the same statistics over the entire
                          sequence; in this case, the op will NOT collapse  
                          sample's internal state
    
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)
        '''
        def token_level_convergence_analysis(data):
            
            #TODO: see whether this code can be simplified and generalized
            if loc is None:
                recast_data = self._recast_data(data, where)
            else:
                recast_data = data

            self._prepare_destinations(recast_data, self.dirs.analysis_out)            
            
            for sample in recast_data:
                
                if loc is None:
                    arr = sample.X
                    freq = self.get_freqs(arr, alphabet)
                    nloc = 'overall'
                    fname = f'{sample.name}_{where}_tokenwise_frequency'

                else:     
                    arr = sample[where]
                    design = self._infer_design(where)
                    #array internal state has to be collapsed for this calculation
                    if not sample._is_collapsed:
                        msg = f"<frequency_summary> op will collapse sample {sample.name}'s internal state"
                        self.logger.warning(msg)
                        sample._collapse_internal_state()
                    
                    #initialize the frequency array: 3D array to be reduced 
                    #along axis 0 at the end
                    maxlen = self._find_max_len(design, loc)
                    freq = np.zeros((len(design), len(alphabet), maxlen),
                                     dtype=np.float32)
                    
                    for i,template in enumerate(design):                    
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)
    
                        #calculated weighed contributions of each design
                        #to the overall frequency array
                        norm = np.divide(np.sum(row_mask), arr.shape[0])
                        freq[i,:,:len(col_mask)] = norm * np.nan_to_num(
                            
                            self.get_freqs(arr[row_mask][:,col_mask], 
                                           alphabet) 
                        )
    
                    #reduce back to a 2D array and plot/save
                    freq = np.sum(freq, axis=0)
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_{where}_reg_{nloc}_tokenwise_frequency'
                
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.tokenwise_frequency(freq, alphabet, where,
                                                           nloc, basename
                                                          )
                if save_txt:          
                    np.savetxt(basename + '.csv',
                               freq,
                               delimiter=',')
                
            return data
        return token_level_convergence_analysis
    
class QScoreAnalysis(Handler):
    
    def __init__(self, *args):
        super(QScoreAnalysis, self).__init__(*args)

    def __repr__(self):
        return '<QScoreAnalysis object>'

    def q_summary(self, loc=None, save_txt=False):
        '''
        For each sample in Data, compute some basic Q score statistics.
    	For each position in regions specified by 'loc', computes the mean and standard deviation
        of Q scores. Plots the results in the parser output folder as specified by config.py.
        Optionally, the data can also be written to a txt file.
        	    	
        Parameters:					  
                     loc: if not None: a list of ints to specify regions to be 
                          analyzed; in this case, the op will collapse sample's
                          internal state (see explanation for Data objects)
                          
                          if left None: get the same statistics over the entire
                          sequence; in this case, the op will NOT collapse  
                          sample's internal state
                          
                save_txt: if True, the data will be written to a txt file saved
                          in the same folder as the .png and .svg plots						 
                          						  							  
        Returns:
                Data object (no transformation)	
        '''            
        def q_score_summary(data):
            
            if loc is None:
                recast_data = self._recast_data(data, 'Q')
            else:
                recast_data = data        
                
            self._prepare_destinations(recast_data, self.dirs.analysis_out)
            
            for sample in recast_data:
            
                if loc is None:
                    relevant_arr = sample.X.astype(np.float32)
                    nloc = 'overall'
                    fname = f'{sample.name}_q_score_summary'
                    
                else:
                    arr = sample.Q
                    if not sample._is_collapsed:
                        msg = f"<q_score_summary> op will collapse sample {sample.name}'s internal state"
                        self.logger.info(msg)
                        sample._collapse_internal_state()
                                    
                    maxlen = self._find_max_len(self.D_design, loc)
                    #iterate over templates and append all of the relevant arr views to this array
                    #relevant view: masked (row/columnwise) arr
                    relevant_arr = []
                    
                    for i,template in enumerate(self.D_design):
                        
                        row_mask = sample._internal_state[:,i]
                        col_mask = template(loc, return_mask=True)                
                        
                        arr_view = np.zeros((np.sum(row_mask), maxlen), dtype=np.float32)
                        arr_view[:,:len(col_mask)] = arr[row_mask][:,col_mask]
                        relevant_arr.append(arr_view)
    
                    #assemble into a single array
                    relevant_arr = np.vstack(relevant_arr)
                    
                    nloc =  ', '.join(str(x + 1) for x in loc)
                    fname = f'{sample.name}_reg{nloc}_q_score_summary'
                    
                #mask out pads (0) as nans for nanmean/nanstd statistics
                relevant_arr[relevant_arr == 0] = np.nan
                
                #get the stats; plot
                q_mean = np.nanmean(relevant_arr, axis=0)
                q_std = np.nanstd(relevant_arr, axis=0)
    
                destination = os.path.join(self.dirs.analysis_out, sample.name)
                basename = os.path.join(destination, fname)
                Plotter.SequencingData.Q_score_summary(q_mean, q_std, nloc, basename)  
    
                if save_txt:
                    q = np.vstack((q_mean, q_std))
                    np.savetxt(basename + '.csv',
                               q.T,
                               delimiter=',',
                               header='Q mean, Q std')  
                
            return data
        return q_score_summary
           
    