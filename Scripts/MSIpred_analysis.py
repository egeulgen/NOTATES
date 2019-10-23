#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' MSIpred Analyses
Perform MSI status prediction using MSIpred
for the given sample
'''

##################################################
## Project: NOTATES
## Script purpose: Perform MSIpred Analysis
## Date: Oct 22, 2019
## Author: Ege Ulgen
##################################################
import sys
import MSIpred as mp
# import pandas as pd
# pd.set_option('display.max_row', 5) 

def main():
    data_sources_dir = sys.argv[1]
    exome_size = float(sys.argv[2])

    ## Initialization of the MAF file object
    maf_obj = mp.Raw_Maf(maf_path='MSIpred/somatic_maf.maf')

    ## Generate an annotated MAF file for further analysis 
    ## by adding one extra column called “In_repeats”, which indicates whether 
    ## mutation events happen in simple repeats region, to the original MAF file 
    ## given the 'simpleRepeat.txt' file
    maf_obj.create_tagged_maf(ref_repeats_file = data_sources_dir + '/simpleRepeat.txt', 
                              tagged_maf_file  = 'MSIpred/tagged_somatic.maf')

    ## Initialization of an annotated MAF file object using the annotated MAF file 
    tagged_maf_obj = mp.Tagged_Maf(tagged_maf_path = 'MSIpred/tagged_somatic.maf')

    ## Create a feature dataframe used for further MSI prediction given size (Mb) 
    ## of captured exome sequence from which MAF file is generated.
    maf_features = tagged_maf_obj.make_feature_table(exome_size=exome_size)
    maf_features.to_csv('MSIpred/maf_feats.csv')

    ## Predict tumor MSI status (MSS or MSI-H) using a SVM classifier given the 
    ## feature dataframe obtained in last step. A pandas dataframe containing 
    ## predicted MSI status for all tumors in the very beginning MAF file will be 
    ## obtained. If not specified (svm_model=None), the default svm model will be 
    ## used for prediction
    predicted_MSI = mp.msi_prediction(feature_table=maf_features, svm_model=None)

    ## add POLE def.
    predicted_MSI['Likely_POLE_deficiency'] = (maf_features.iloc[0]['SNP'] > 60) & (maf_features.iloc[0]['INDEL_R'] < 0.18)
    
    predicted_MSI.to_csv('MSIpred/MSIpred_prediction.csv')

if __name__ == '__main__':
    main()
    exit(0)