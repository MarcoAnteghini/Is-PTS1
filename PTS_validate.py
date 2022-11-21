#!/usr/bin/env python

from time import time
import pickle
import os
import numpy as np
from Bio import SeqIO
import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import tensorflow as tf

start_time=time()

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


tf.set_random_seed(42)
np.random.seed(42)

fastafile=sys.argv[1]

# create DATA folder
os.makedirs(fastafile[:-6] + "/DATA/", exist_ok=True)


from unirep import babbler1900 as babbler

# Where model weights are stored. 
# Here we consider the directory where the script is executed. Insert the proper path below if needed.
MODEL_WEIGHT_PATH = "."
b = babbler('DATA/1900_weights')

print('UniRep Encoding:')
with open(fastafile, "r") as handle:
    if is_fasta(fastafile):
        fasta = SeqIO.parse(handle, "fasta")
        d={}
        for record in fasta:
            try:
                sequence=record.seq
                ids=record.id
                d[ids.split('|')[1]] = sequence
                    
            except:
                print('Something wrong','\n','fasta file should start with >sp|ID|ORGANISM')
                pass

for k in d.keys():
    if len(d[k])>1200:
        print('\n'*10+k+' Seqeunce too long, please remove sequences longer then 1000')
        exit()

c=0
for keys in d.keys():
    try:
        ur = b.get_rep(d[keys])
        tosave1 = np.asarray(ur[0])
        #        tosave2 = np.asarray(ur[1])
        #        tosave3 = np.asarray(ur[2])
# We here save just one of the 3 arrays that Unirep produces.
        np.save(fastafile[:-6] + "/DATA/" + keys+'_UniRep1', tosave1)
        #        np.save(ids.split('|')[1]+'_UniRep2', tosave2)
        #        np.save(ids.split('|')[1]+'_UniRep3', tosave3)
        c=c+1
        print('ID:',keys,' '*20,c,'/',len(d))
    except:
        pass
#        print('Not encoded:',to_check)
        
from subprocess import Popen, PIPE
print('\n'*2)
print('SeqVec encoding:','\n')

p1 = Popen(["seqvec", "-i", fastafile, "-o", fastafile[:-6]+str('_seqvec.npz'),"--protein"], stdout=PIPE)

p1.communicate()

#comment/uncomment the model you want to use
filename='Is-PTS1_model.sav'
#filename='SVM_model.sav'


model = pickle.load(open(filename, 'rb'))

unireps = {}
for filenames in os.listdir(fastafile[:-6] + '/DATA'):
    if filenames.endswith('UniRep1.npy'):
        unireps[filenames.split('_')[0]] = np.load(fastafile[:-6] + "/DATA/" + filenames)

seqvecs = {}
seqvec_archive = np.load(fastafile[:-6]+str('_seqvec.npz'), allow_pickle=True)
for uniprotid in seqvec_archive.files:
    seqvecs[uniprotid] = seqvec_archive[uniprotid]

final_d={}
for keys in seqvecs.keys():
    for k in unireps.keys():
        if k==keys:
            final_d[k]= np.concatenate([unireps[k], seqvecs[k]])

true_pero, false_pero = [],[]
for uniprotid in final_d:

    pred = model.predict(final_d[uniprotid].reshape(1, -1))[0]
    if pred==1:
        true_pero.append(uniprotid)
    else:
        false_pero.append(uniprotid)
print('\n'*2)
print('True peroxisomal PTS: ',true_pero)
print('False peroxisomal PTS: ',false_pero)



notenc=set(seqvecs.keys())-set(unireps.keys())

    



output = open(fastafile[:-6]+str('_output.txt'), 'w')
output.write("%s\n" % 'True pero PTS:')
for e in true_pero:
    output.write("%s\n" % e)
output.write("%s\n" % 'False pero PTS:')
for e in false_pero:
    output.write("%s\n" % e)
output.write("%s\n" % 'not encoded:')
for e in notenc:
    output.write("%s\n" % e)
output.close()

print((time()-start_time)/60)
