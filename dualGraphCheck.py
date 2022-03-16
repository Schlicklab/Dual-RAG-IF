#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 22:31:46 2020

@author: qiyaozhu
"""
import random as RANDOM
import os
import os.path
import sys
import time
from functools import partial
import multiprocessing
from ClassesFunctions import *
from dualGraphs import *


def graphFinder(jobID):
    
    RNA = getCTInfo("tmpRNAfold"+jobID+".ct")
    os.system("rm -rf tmpRNAfold"+jobID+".ct")
    
    countHelices(RNA) 
    changeHelices(RNA)
    RNA.makeMatrices()
    connectHelices(RNA)
    for i in range(0,len(RNA.adjMatrix)): # S.J. 07/11/2018 - to keep track of vertexOrder
        vertexOrder.append(0)
        
    success, graph = calcEigen(RNA)
    correctHNumbers(RNA)
    if len(RNA.adjMatrix)==1 or len(RNA.adjMatrix)>9:
        print ("No matching graph exists because vertex number is either 1 or greater than 10.")
        return None
    elif success == 0: # no graph ID was assigned as eigen values not in the library S.J. 11/09/2017
        print ("No matching graph exists (even if the vertex number is between 2 and 9).")
        return None
    else:
        return graph
    

def nupackCheck(survivor):
    
    seq = survivor[0]
    jobID = survivor[1]
    
    with open("tmpRNAfold"+jobID+".in","w") as f:
        f.write(">seq\n")
        f.write(seq)
    
    os.system("mfe -pseudo -material rna tmpRNAfold"+jobID+" 2>/dev/null ")
    with open("tmpRNAfold"+jobID+".mfe", 'r') as f:
        fold = f.readlines()[16]
    with open("tmpRNAfold"+jobID+".mfe", 'w') as f:
        f.write(">seq\n")
        f.write(seq + "\n")
        f.write(fold)
    os.system("dot2ct tmpRNAfold"+jobID+".mfe tmpRNAfold"+jobID+".ct")    
    os.system("ct2dot tmpRNAfold"+jobID+".ct 1 tmpRNAfold"+jobID+".out")
    with open("tmpRNAfold"+jobID+".out") as f:
        lines = f.readlines()
        nupackFold = lines[2].split('\n')[0]
    os.system("rm -rf tmpRNAfold"+jobID+".in tmpRNAfold"+jobID+".mfe tmpRNAfold"+jobID+".out")
    
    return graphFinder(jobID), nupackFold


def getSeqJob(heaven):
    
    Survivors = []
    
    with open(heaven, 'r') as f:
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i][0] != '>':
            continue
        else:
            seq = lines[i+1].split(" ")[0]
            fold = lines[i+2].split("\n")[0]            
            jobID = str(RANDOM.randint(10000,99999))
            jobID = jobID+str(time.time()).split('.')[1]
            list1= ['a','b','c','d','e','f','g','h','i','j']
            list2= [1,2,3,4,5,6,7,8,9,0]
            jobID = jobID + RANDOM.choice(list1)+ str(RANDOM.choice(list2))
            jobID = jobID + RANDOM.choice(list1)+ str(RANDOM.choice(list2))
            jobID = jobID + RANDOM.choice(list1)+ str(RANDOM.choice(list2))
            Survivors.append([seq, jobID, fold, heaven])
            
    return Survivors


def IPknotvsNUPACK(survivor):
    
    seq = survivor[0]
    heaven = survivor[3]
    
    name = heaven.split('h')[0]
    filename = name + 'Sequences.txt'
    design = str(name.split('_')[0]) + '_' + str(name.split('_')[1])
    
    nupackGraph, nupackFold = nupackCheck(survivor)
    if nupackGraph == design:
        with open(filename, 'a+') as f:
            f.write('>\n')
            f.write(seq + '\n')
            f.write(nupackFold + ' NUPACK\n')
        print("Both right")
    else:
        print("IPknot right, NUPACK graph: " + nupackGraph + ", NUPACK fold:\n" + seq + "\n" + nupackFold)
            
           
def doubleCheck(heaven):
    
    if not os.path.isfile(heaven):
      print("Heaven file not exist...")
      sys.exit()
    
    Survivors = getSeqJob(heaven)
    pool = multiprocessing.Pool(4)
    pool.map( IPknotvsNUPACK, Survivors )
    pool.close()
    pool.join()
    
            
            
    
    