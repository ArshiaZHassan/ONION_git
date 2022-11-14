# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 23:32:21 2020

@author: Arshia Zernab Hassan
"""

# Imports packages
import sys
import argparse
import pandas as pd
from sklearn.preprocessing import StandardScaler

# Create parser for user input
def create_arg_parser():
    parser = argparse.ArgumentParser(prog='pre_process.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file",
                        default=argparse.SUPPRESS,
                        type=str,
                        required=True,
                        help="Path to data to be visualized [required]")
    parser.add_argument("-o", "--output_file",
                        default=argparse.SUPPRESS,
                        type=str,
                        required=True,
                        help="Path to output file [required]")
    parser.add_argument("-rst", "--row_standardize",
                        action="store_true",
                        help="Row standardize data [optional]")
    parser.add_argument("-cst", "--col_standardize",
                        action="store_true",
                        help="Column standardize data [optional]")
    parser.add_argument("-mms", "--min_max_scale",
                        action="store_true",
                        help="min-max scaling of data (-1 to 1) [optional]")
    parser.add_argument("-clp", "--clip",
                        action="store_true",
                        help="clip data (-4 to 4) [optional]")
    return parser

'''
Function:
    row_standardize()
Arguments:
    data : dataset to apply row-wise standardization. type: pandas dataframe (required)
Return: 
    row-wise standardized dataset. type: pandas dataframe
Description:
    - transpose data
    - apply standard scaler to data to standardize columns
    - re-transpose data
#Standardize features by removing the mean and scaling to unit variance
#The standard score(z-score) of a sample x is calculated as: z = (x - u) / s
#https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
''' 
def row_standardize(data):
    data_T = pd.DataFrame(data.transpose())
    scalerT = StandardScaler()
    scalerT.fit(data_T)
    data_T_c_st = pd.DataFrame(scalerT.transform(data_T), index = data_T.index, columns = data_T.columns)
    data_r_st = pd.DataFrame(data_T_c_st.transpose())
    print('Row standardized')    
    return data_r_st
    
'''
Function:
    column_standardize()
Arguments:
    data : dataset for applying column-wise standardization. type: pandas dataframe (required)
Return: 
    column-wise standardized dataset. type: pandas dataframe
Description:
    - apply standard scaler to data to standardize columns
#Standardize features by removing the mean and scaling to unit variance
#The standard score(z-score) of a sample x is calculated as: z = (x - u) / s
#https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
''' 
def column_standardize(data):
    scaler = StandardScaler()
    scaler.fit(data)
    data_c_st = pd.DataFrame(scaler.transform(data), index = data.index, columns = data.columns)
    print('Column standardized')
    return data_c_st

'''
Function:
    min_max_scale_neg1_to_1()
Arguments:
    data : dataset for applying min-max scaling. type: pandas dataframe (required)
Return: 
    min-max scaled dataset. type: pandas dataframe
Description:
    -apply min-max scaling to entire data (range [-1,1])
    
                x - minimum(data)
      x'= 2 * --------------------------- - 1
              maximum(data)-minimum(data)
'''
def min_max_scale_neg1_to_1(data):
    data_mms = pd.DataFrame()
    mx= data.values.max()    
    mn= data.values.min()    
    data_mms = (((data - mn)*2) / (mx - mn))-1
    print('Min-max scaled to range [-1,1]')
    return data_mms

'''
Function:
    clip_range()
Arguments:
    data : dataset for applying clipping clipping. type: pandas dataframe (required)
    lb : lower bound for clipping. type: float (required)
    ub : upper bound for clipping. type: float (required)
Return: 
    clipped dataset
Description:
    -clip data to range [lb,ub]
'''  
def clip_range(data, lb, ub):
    data_clip = data.clip(lb, ub)
    print('Clipped to range [' + str(lb) +',' + str(ub) +']')
    return data_clip


'''
Function:
    pre_process()
Arguments:
    input_file : path to the data file. type: str (required)
    output_file : path to output file. type: str (required)
    pp_list : sequential list of pre-processing steps. type: list of strings
Return: 
    none
Description:
    -load data as pandas dataframe from input_file (should have no header or index), 
    -replace NAN values with 0 
    -run pre-processing steps in the order of cmd. line arguments 
    (column-wise stardardization or row-wise stardardization or min-max scale or clip), 
    -save pre-processed data in output_file
'''  
def pre_process(input_file,output_file, pp_list):
    
    #works on data created by this code 
    data = pd.read_csv(input_file,sep = "\t",index_col=0)
    
    #works on original data
    #data = pd.read_csv(input_file,sep = "\t")
    
    data_ = data.fillna(0)
    
    for process in pp_list:
    
        if process == '-rst':
            data_ = row_standardize(data_)
        elif process == '-cst':
            data_ = column_standardize(data_)
        elif process == '-mms':
            data_ = min_max_scale_neg1_to_1(data_)
        elif process == '-clp':
            data_ = clip_range(data_, -4, 4)
    
    data_.to_csv(output_file, sep='\t')
    
    return

'''
-Parse command line arguments list and create (argument name: value) dictionary
-extract arguments for pre-processing steps in the given order and save them in list
-call pre_process function with input file path, output file path and pre-processing steps list as arguments
'''
if __name__ == '__main__':
    parser = create_arg_parser()
    args = parser.parse_args()
    args = vars(args)
    
    arg_list = sys.argv    
    pp_list = [i for i in arg_list if i in ['-mms', '-rst', '-cst', '-clp']]
    print('Pre-processing commands:' + str(pp_list))
    
    pre_process(args['input_file'], args['output_file'], pp_list)

'''
run instruction:
    pre_process_2.py -i <input file path> -o <output file path> <pre-processing step arguments in the sequence of execution>
run example:
    python3 pre_process_2.py -i input/depmap_q3_2019.tsv -o output/depmap_q3_2019_rst_clp_mms.tsv -rst -clp -mms  
'''