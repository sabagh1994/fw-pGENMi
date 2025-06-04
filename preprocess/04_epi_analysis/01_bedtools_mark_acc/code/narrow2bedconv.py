import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
dir_ = f'{SCRIPT_DIR}/..' # 01_bedtools_mark_acc

path= f'{dir_}/final/narrowPeak/all/'
outdir= f'{dir_}/final/bed5/all/'
os.makedirs(outdir, exist_ok= True)
for filename in os.listdir(path):
    if not filename.startswith('.'):
        df= pd.read_csv(path+filename, sep= '\t', header= None)
        df.drop([4,5,6,7,8], axis= 1, inplace= True)
        df.to_csv(outdir+filename.split('.')[0]+'.bed', sep= '\t', index= False, header= False)
