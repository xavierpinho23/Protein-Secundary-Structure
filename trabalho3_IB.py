# -*- coding: utf-8 -*-
"""
Protein Secondary Structure
Goals:
    - Secondary Structure Forecast
    - Implementation of the algorithm proposed by Nussinov

Xavier Pinho & Jorge Melo - Introduction to Bioinformatics, University of Coimbra, 2018/2019
"""

import pandas as pd
import numpy as np
import urllib
import os 

print("------- Importation Complete ----------")
home = os.getcwd()
os.chdir(home)

def delta(i,j):
    """function that analyse the complementarity
    return 1 if i and j are complementary
    return 0 if not"""
    if i=='A'   and j=='U':
        return 1;
    elif i=='U' and j=='A':
        return 1;
    elif i=='G' and j=='C':
        return 1;
    elif i=='C' and j=='G':
        return 1;
    else:
        return 0;

def DP(sequence):
    """Function to build the DP matrix """
    size = len(sequence)
    matrix = np.zeros((size,size))
    scores = np.zeros((size,size))
    for n in range(1,size):
        for j in range(n,size):
            i = j-n
            case1 = matrix[i+1,j-1] + delta(sequence[i],sequence[j])
            case2 = matrix[i+1,j]
            case3 = matrix[i,j-1]
            if i + 3 <= j:
                tmp = []
                for k in range(i+1,j):
                    tmp.append(matrix[i,k] + matrix[k+1,j])
                case4 = max(tmp)
                matrix[i,j] = max(case1,case2,case3,case4)
            else:
                matrix[i,j] = max(case1,case2,case3)
            if matrix[i,j]==case1:
                scores[i,j]=1
            elif matrix[i,j]==case2:
                scores[i, j]=2
            elif matrix[i,j]==case3:
                scores[i, j]=3
            elif matrix[i,j]==case4:
                scores[i, j]=4

    #np.save(home+"/DP_matrix",matrix)
    return matrix,scores
 
def traceback(matrix,sequence,i,j,pair):
    """Function that compute the traceback"""
    if i <j:
        # case1
        if matrix[i,j] == (matrix[i+1,j-1] + delta(sequence[i],sequence[j])):
            pair.append([i,j])
            traceback(matrix,sequence, i+1,j-1,pair)    
        # case2
        elif matrix[i,j] == matrix[i+1,j]:
            traceback(matrix,sequence, i+1,j,pair)
        # case3
        elif matrix[i,j] == matrix[i,j-1]: 
            traceback(matrix,sequence, i,j-1,pair)
        # case4
        else:
            for k in range(i+1,j):
                if matrix[i,j] == (matrix[i,k] + matrix[k+1,j]):
                    traceback(matrix,sequence, i,k,pair)
                    traceback(matrix,sequence, k+1,j,pair)
                    break
    return pair

##Let's test for 4 different cases
#case 1 - tRNA
tRNA = "GGTTCCATGGTGTAGTGGTTATCACATCTGCTTTACACGCAGAAGGTCCTGGGTTCAAGCCCCAGTGGAACCA"
matrix,key_scores = DP(tRNA)
pairs = traceback(matrix,tRNA,0,len(tRNA)-1,[])
out = ["."]*len(tRNA)
for i in range(len(pairs)):
    out[pairs[i][0]] = "("
    out[pairs[i][1]] = ")"
print("tRNA: \n", tRNA) 
print("sequence: \n", out)
print("max # of folding pairs: ",len(pairs),"\n")
print("Pairs: ", pairs,"\n\n")

#case 2
file = open(home+"/sequence.txt")
seq1 = file.read()
seq1 = seq1.replace("\n","")
matrix1,key_scores1 = DP(seq1)
pairs1 = traceback(matrix1,seq1,0,len(seq1)-1,[])
out1 = ["."]*len(seq1)
for i in range(len(pairs1)):
    out1[pairs1[i][0]] = "("
    out1[pairs1[i][1]] = ")"
print("AB063105: \n", seq1 ) 
print("sequence: \n", out1)
print("max # of folding pairs: ",len(pairs1),"\n")
print("Pairs: ",pairs1,"\n\n")              

#case 3
seq_Q3YRK7 = open(home+'/seq_Q3YRK7.txt')
seq2 = seq_Q3YRK7.read()
seq2 = seq2.replace("\n","")
matrix2,key_scores2 = DP(seq2)
pairs2 = traceback(matrix2,seq2,0,len(seq2)-1,[])
out2 = ["."]*len(seq2)
for i in range(len(pairs2)):
    out2[pairs2[i][0]] = "("
    out2[pairs2[i][1]] = ")"
print("Q3YRK7: \n",   seq2)   
print("sequence: \n", out2)
print ("max # of folding pairs: ",len(pairs2), "\n")
print('pairs: ',pairs2, "\n\n")

#case 4
seq_Q4A597 = open(home+'/seq_Q4A597.txt')
seq3 = seq_Q4A597.read()
seq3 = seq3.replace("\n","")
matrix3,key_scores3 = DP(seq3)
pairs3=traceback(matrix3,seq3,0,len(seq3)-1,[])
out3 = ["."]*len(seq3)
for i in range(len(pairs3)):
    out3[pairs3[i][0]] = "("
    out3[pairs3[i][1]] = ")"
print("Q4A597: \n",   seq3)
print("sequence: \n", out2)
print ("max # of folding pairs: ",len(pairs3), "\n")
print('pairs: ',pairs3, "\n\n")


