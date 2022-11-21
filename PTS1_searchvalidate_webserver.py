#!/usr/bin/env python

import sys
import re
from Bio import SeqIO
from time import time
import pickle
import os
import numpy as np
import pandas as pd
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf



import re    
#check if sample_str is a number or string
def is_number_or_float(sample_str):
    sample_str = str(sample_str)
    print ("is_number_or_float::: type(sample_str)", type(sample_str))

    ''' Returns True if the string contains only
        number or float '''
    result = True
    if re.search("[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$", sample_str) is None:
        result = False
    return result



def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

        
print(isfloat('s12'))
print(isfloat('1.123'))



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

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


### Added after acceptence
def check_fasta(filename):
    with open(filename) as f:
        fasta = SeqIO.parse(f,"fasta")
        for record in fasta:
            if 'sp' in record.id:
                #print(record.id)
                return True
        else:
            return False
       
def generate_fasta(filename, newfasta):
    with open(newfasta,'w') as nfasta, open(filename) as f:
        fasta = SeqIO.parse(f,"fasta")
        for record in fasta:
            nfasta.write('>sp|'+str(record.id)+'|ORGANISM'+'\n'+str(record.seq)+'\n')




# Stores file one for input checking.
fastafile=sys.argv[1]
#newfasta=fastafile.rsplit('/', 1)[0] + '/newfastafile.fasta'
newfasta=fastafile[:-6] + '_newfastafile.fasta'


print ("fastafile =================== ", fastafile)
print ("newfasta =================== ", newfasta)


## check fasta
if is_fasta(fastafile)==False:
    print ("Please submit a correct fasta file")
    output = open(sys.argv[1][:-6]+str('_output.txt'), 'w')
    output.write("Please submit a correct fasta file\n")
    output.close()
    exit()

if check_fasta(fastafile)==False:
    generate_fasta(fastafile, newfasta)
    fastafile=newfasta




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
        match = re.search(r'[ASCNPHTGEQ][RKHQNSL][LAMIVF]$', sequence[-int(leng):])
    if match:
        ms[name]=match.group()
    if not match:
        non_PTS.append(name)

with open(fastafile, "r") as handle:
    if is_fasta(fastafile):
        fasta = SeqIO.parse(handle, "fasta")
        d={}
        for k,record in zip(ms.keys(),fasta):
            #print(k,record.id.split('|')[1])
            if k==record.id.split('|')[1]:
                try:
                    sequence=record.seq
                    ids=record.id
                    d[ids.split('|')[1]] = sequence
                except:
                    print('Something wrong','\n','fasta file should start with >sp|ID|ORGANISM')
                    pass

print ("()"*50)
print ("()"*50)
print(d)
print ("()"*50)
print ("()"*50)
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
os.makedirs(fastafile[:-6] + "/DATA/", exist_ok=True)

tf.set_random_seed(42)
np.random.seed(42)


from unirep import babbler1900 as babbler

# Where model weights are stored.
# Here we consider the directory where the script is executed. Insert the proper path below if needed.
MODEL_WEIGHT_PATH = "."
b = babbler('/data/p290092/PTS_validate/1900_weights')

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
                     
        c=0
        for keys in d:
            try:
                ur = b.get_rep(d[keys])
                tosave1 = np.asarray(ur[0])
        #        tosave2 = np.asarray(ur[1])
        #        tosave3 = np.asarray(ur[2])

                np.save(fastafile[:-6] + "/DATA/" + keys+'_UniRep1', tosave1)
        #        np.save(ids.split('|')[1]+'_UniRep2', tosave2)
        #        np.save(ids.split('|')[1]+'_UniRep3', tosave3)
                c=c+1
                print('ID:',keys,' '*20,c,'/',len(d))
            except:
                pass
#        print('Not encoded:',to_check)
                       
    else:
        print('YOU HAVE TO INPUT A CORRECT FASTA FILE')
        
from subprocess import Popen, PIPE
print('\n'*2)
print('SeqVec encoding:','\n')

p1 = Popen(["seqvec", "-i", fastafile, "-o", fastafile[:-6]+str('_seqvec.npz'),"--protein"], stdout=PIPE)

p1.communicate()

#comment/uncomment the model you want to use
filename='test_pts_pred_LR_all_features_model.sav'
#filename='SVM_model.sav'


LR_model = pickle.load(open("/data/p290092/PTS_validate/" + filename, 'rb'))

unireps = {}
for filenames in os.listdir(fastafile[:-6] + '/DATA'):
    if filenames.endswith('UniRep1.npy'):
        unireps[filenames.split('_')[0]] = np.load(fastafile[:-6] + "/DATA/" + filenames)

seqvecs = {}
seqvec_archive = np.load(fastafile[:-6]+str('_seqvec.npz'), allow_pickle=True)
for uniprotid in seqvec_archive.files:
    seqvecs[uniprotid] = seqvec_archive[uniprotid]


print ("seqvecs == ", seqvecs)
print ("unireps == ", unireps)

print ("ms === ", ms)
final_d={}
for keys in seqvecs.keys():
    print ("keys ==== ", keys)
    for k in unireps.keys():
        print ("k ====== ", k)
        if k==keys and k in ms:
            final_d[k]= np.concatenate([unireps[k], seqvecs[k]])

print ("final_d ========= ", final_d)

true_pero, false_pero, protID, probs, true_pero2, false_pero2  = [],[],[],[],[],[]
for uniprotid in final_d:

    print ("uniprotid === ", uniprotid)
    protID.append(uniprotid)
    pred = LR_model.predict(final_d[uniprotid].reshape(1, -1))[0]
    proba = LR_model.predict_proba(final_d[uniprotid].reshape(1, -1))[0]
    probs.append(proba)
    if pred==1:
        true_pero.append((uniprotid,float("%.3f" % proba[1])))
    else:
        false_pero.append((uniprotid,float("%.3f" % proba[0])))
for e in probs:
    true_pero2.append(e[1])
    false_pero2.append(e[0])
head=['ProteinID','True PTS1','False PTS1']
df_output=pd.DataFrame(columns = head)
df_output['ProteinID']=protID
df_output['True PTS1']=true_pero2
df_output['False PTS1']=false_pero2

print('\n'*2)
print('True peroxisomal PTS: ',true_pero)
print('False peroxisomal PTS: ',false_pero)
print('PTS1 not found:',non_PTS)

print ("*" * 50)
print ("true_pero === ", true_pero)
print ("false_pero === ", false_pero)
print ("protID === ", protID)
print ("probs === ", probs)
print ("true_pero2 === ", true_pero2)
print ("false_pero2 === ", false_pero2)
print ("*" * 50)
print ("*" * 50)


notenc=set(seqvecs.keys())-set(unireps.keys())

print(notenc)


df_output.to_csv(sys.argv[1][:-6]+str('_output.csv'), index=False)




print ("*" * 50)
print ("*" * 50)
output = open(sys.argv[1][:-6]+str('_output.txt'), 'w')
output.write("%s\n" % 'True pero PTS:')
for e in true_pero:
    #print (e, is_number_or_float(e))
    #print (str(e), str(e).isdigit(), isinstance(str(e), float), isinstance(e, float))
    #print ("True pero PTS:", type(e), e, e.isdigit())
    #print (e, type(e), is_number_or_float(e), isfloat(e))
    #if isfloat(e):
    #    output.write("%.3f\n" % e)
    #else:
    #    output.write("%s\n" % e)
    print ("e === ", e)
    output.write("%s\n" % str(e))

output.write("%s\n" % 'False pero PTS:')
for e in false_pero:
    #print (str(e), str(e).isdigit(), isinstance(str(e), float))
    #print ("False pero PTS:", type(e), e, e.isdigit())
    #print (e, type(e), is_number_or_float(e), isfloat(e))
    #if isfloat(e):
    #    output.write("%.3f\n" % e)
    #else:
    #    output.write("%s\n" % e)
    print ("e === ", e)
    output.write("%s\n" % str(e))

output.write("%s\n" % 'PTS1 not detected:')
for e in non_PTS:
    #print (e.isnumeric())
    #print (str(e), str(e).isdigit(), isinstance(str(e), float))
    #print ("PTS1 not detected:", type(e), e, e.isdigit())
    #print (e, type(e), is_number_or_float(e), isfloat(e))
    #if isfloat(e):
    #    output.write("%.3f\n" % e)
    #else:
    #    output.write("%s\n" % e)
    print ("e === ", e)
    output.write("%s\n" % str(e))

output.write("%s\n" % 'not encoded:')
for e in notenc:
    #print (str(e), str(e).isdigit(), isinstance(str(e), float))
    #print ("not encoded:", type(e), e, e.isdigit())
    #print (e, type(e), is_number_or_float(e), isfloat(e))
    #if isfloat(e):
    #    output.write("%.3f\n" % e)
    #else:
    #    output.write("%s\n" % e)
    print ("e === ", e)
    output.write("%s\n" % str(e))


output.close()


print((time()-start_time)/60)


