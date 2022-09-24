#!/usr/bin/env python

import sys
import re
from Bio import SeqIO

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print("Performs a PTS search in the C-terminal region of a fasta sequence")
		print("By Marco Anteghini\n")
		print("Usage: " + sys.argv[0] + " <sequences.fasta>" + ' species' + ' lenght_of_Cterminal_to_keep')
		print("Examples: " + sys.argv[0] + " mySeqs.fasta" + ' yeast' + ' 70' )
		exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(4) # Checks if the number of arguments are correct.

# Stores file one for input checking.
fastafile=argv = sys.argv[1] #file.fasta
species= sys.argv[2] # either choose between yeast or other
cterm_len=sys.argv[3] #the lenght of the Cterminal part the user want to keep

# Possible inputs as 3rd argument:
# yeast
# other


#Known PTS signal
#[SACHEQ]-[KRH]-[LAF] yeast
#[ASCNPHTG]-[RKHQNSL]-[LMIVF] euK

fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
leng=70
ms={}
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    if '|' in name:
        name=name.split('|')[1]
    if species=='other':
        #print(sequence+'\n')
        #print(sequence[-int(cterm_len):])
        match = re.search(r'[ASCNPHTG][RKHQNSL][LMIVF]$', sequence[-int(cterm_len):])
    if species=='yeast':
        match = re.search(r'[SACHEQ][KRH][LAF]$', sequence[-int(cterm_len):])
    if match:
        ms[name]=match.group()

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
print('\n')

