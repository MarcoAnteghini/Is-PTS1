#!/usr/bin/env python

import sys
import re
from Bio import SeqIO
from time import time
import pickle
import os
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf


# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print("Performs a PTS search in the C-terminal region of a fasta sequence")
		print("By Marco Anteghini\n")
		print("Usage: " + sys.argv[0] + " <sequences.fasta>")
		print("Examples: " + sys.argv[0] + " mySeqs.fasta")
		exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#===========================================================================================================
# Main program code:

start_time=time()

# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.

# Stores file one for input checking.
fastafile=sys.argv[1] #file.fasta

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)



#Known PTS signal
#[SACHEQ]-[KRH]-[LAF] yeast
#[ASCNPHTG]-[RKHQNSL]-[LMIVF] euK

species='other'
fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
leng=40
ms={}
non_PTS=[]
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    if '|' in name:
        name=name.split('|')[1]
    if species=='other':
        #print(sequence+'\n')
        #print(sequence[-int(cterm_len):])
        match = re.search(r'[ASCNPHTG][RKHQNSL][LMIVF]$', sequence[-int(leng):])
    if species=='yeast':
        match = re.search(r'[SACHEQ][KRH][LAF]$', sequence[-int(leng):])
    if match:
        ms[name]=match.group()
    if not match:
        non_PTS.append(name)
#print(ms)
with open(fastafile, "r") as handle:
    if is_fasta(fastafile):
        fasta = SeqIO.parse(handle, "fasta")
        d={}
        for record in fasta:
            #print('GUARDA QUIIIIIIIIIIIIIIIIIIII')
            #print(record.id.split('|')[1])
            if record.id.split('|')[1] in ms.keys():
             #   print('AAAAAAAAAAAAAAAAAAAAAAAASFOJEAPOFJOJFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
             #   print(k,record.id.splip('|')[1])
                try:
                    sequence=record.seq
                    ids=record.id
             #       print(ids+'\n'+sequence)
                    d[ids.split('|')[1]] = sequence
                except:
                    print('Something wrong','\n','fasta file should start with >sp|ID|ORGANISM')
                    pass

print(d)
#d should contain only the matching fasta files

#with open(fastafile) as fasta, open(fastafile[:-6]+'_matches.txt','w') as result:
#    c=0
#    if species=='other':
#        for line in fasta:
#            if line[0] == '>':
#                idx=line.rstrip()[1:]   
#            else:
#                match = re.search(r'[ASCNPHTG][RKHQNSL][LMIVF]$', line[-int(cterm_len):])
#                if match:
#                    result.write(idx+' '+match[-int(cterm_len):])
#                    c+=1
#    if species=='yeast':
#        for line in fasta:
#            if line[0] == '>':
#                idx=line.rstrip()[1:]   
#            else:
#                match = re.search(r'[SACHEQ][KRH][LAF]$', line[-int(cterm_len):])
#                if match:
#                    result.write(idx+' '+match[-int(cterm_len):])
#                    c+=1    
with open(fastafile[:-6]+'_matches.txt','w') as out:
    for k,v in ms.items():
        print(k,v)
        out.write(k+'\t'+v+'\n')


print('_'*15)
print('\n')
print('FOUND:'+' '+str(len(ms))+' PTSs')
print('\n'*2)
print('The output is stored in: '+str(fastafile[:-6])+'_matches.txt')


# create DATA folder
os.makedirs("DATA/", exist_ok=True)

tf.set_random_seed(42)
np.random.seed(42)


from unirep import babbler1900 as babbler

# Where model weights are stored.
# Here we consider the directory where the script is executed. Insert the proper path below if needed.
MODEL_WEIGHT_PATH = "."
b = babbler('1900_weights')

print('UniRep Encoding:')
for k in d.keys():
    if len(d[k])>1200:
        print('\n'*10+k+' Seqeunce too long, please remove sequences longer then 1200')
        exit()

c=0
print(d.keys())
for keys in d.keys():
    try:
        ur = b.get_rep(d[keys])
        tosave1 = np.asarray(ur[0])
        #        tosave2 = np.asarray(ur[1])
        #        tosave3 = np.asarray(ur[2])
# We here save just one of the 3 arrays that Unirep produces.
        np.save("DATA/" + keys+'_UniRep1', tosave1)
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

p1 = Popen(["seqvec", "-i", fastafile, "-o", fastafile[:-6]+str('_seqvec.npz'),"--protein","True"], stdout=PIPE)

p1.communicate()

#comment/uncomment the model you want to use
filename='Is-PTS1_model.sav'
#filename='SVM_model.sav'


model = pickle.load(open(filename, 'rb'))

unireps = {}
for filenames in os.listdir('DATA'):
    if filenames.endswith('UniRep1.npy'):
        unireps[filenames.split('_')[0]] = np.load("DATA/" + filenames)

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
print('PTS1 not found:',non_PTS)


notenc=set(seqvecs.keys())-set(unireps.keys())

    



output = open(fastafile[:-6]+str('_PTS_output.txt'), 'w')
output.write("%s\n" % 'True pero PTS:')
for e in true_pero:
    output.write("%s\n" % e)
output.write("%s\n" % 'False pero PTS:')
for e in false_pero:
    output.write("%s\n" % e)
output.write("%s\n" % 'PTS1 not detected:')
for e in non_PTS:
    output.write("%s\n" % e)
output.write("%s\n" % 'not encoded:')
for e in notenc:
    output.write("%s\n" % e)
output.close()

print((time()-start_time)/60)
