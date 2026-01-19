import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
import argparse
import scipy.io
import pandas as pd
import itertools

def call_rand(ar1,ar2):
    return adjusted_mutual_info_score(ar1,ar2)

def parse_args():
    parser=argparse.ArgumentParser(description='Permutation for connectivity')
    parser.add_argument('--matfile', type=str)
    parser.add_argument('--flag', type=int)
    parser.add_argument('--suffix', type=str)
    args=parser.parse_args()
    return args

if __name__== '__main__':

    output_dir='/data/chamal/projects/natvik/knox_qc_full_06232025/derivatives/community_louvain/'
    inputs=parse_args()
    mat = scipy.io.loadmat(inputs.matfile)
    if inputs.flag==0:
        rands=dict()
        communities=mat['citemp']
        means=[]
        stds=[]
        for col in range(communities.shape[1]+1):
            print('this is col', col)
            r=[]
            
            if col == 0:
                continue
            elif col==1:
                continue
            else:
                for subset in itertools.combinations(range(col),2):
                # print(subset[0],'and',subset[1])
                    r.append(call_rand(communities[subset[0]],communities[subset[1]]))
                # print('mean',np.average(r))
                means.append(np.average(r))
                print(len(means))
                # print('std',np.std(r))
                stds.append(np.std(r))
            # rands[f'{col}']=r
        
        scipy.io.savemat(output_dir+'rand_output.mat',mdict={'means':means,'stds':stds})



    if inputs.flag==1:
        communities=mat['ciu']
        rands=dict()
        for col in range(communities.shape[1]):
            r=[]
            ar1=communities[:,col]
            for col2 in range(communities.shape[1]):
                
                ar2=communities[:,col2]
                r.append(call_rand(ar1,ar2))
            rands[f'{col}']=r
        rands_df=pd.DataFrame(rands)
        scipy.io.savemat(output_dir+'rand_output_'+inputs.suffix+'.mat',mdict={'a':pd.DataFrame(rands).to_numpy()})
    # print(rands_df)