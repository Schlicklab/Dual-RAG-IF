#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 01:04:19 2020

@author: qiyaozhu
"""

import os
import os.path
import sys


# @ Sequences file contains all the mutated candidate sequences
# @ origf file contains the original sequence
# Organize the genetic algorithm candidates in ascending order of mutation numbers compared to the original sequence
def minCount(Sequences, origf):
    
    if not os.path.isfile(Sequences):
      print("Sequence file not exist...")
      sys.exit()
    
    if not os.path.isfile(origf):
      print("Original sequence file not exist...")
      sys.exit()
    
    dic = {}
    num = 0
    
    with open(Sequences, 'r') as f:
        lines = f.readlines()
    
    with open(origf, 'r') as f:
        f.readline()
        origseq = f.readline().split('\n')[0]    
        
    for i in range(len(lines)):
        if lines[i][0] == '>':
            num = num + 1
            count = 0
            seq = lines[i+1].split('\n')[0]
            mutseq = ''
            for j in range(len(seq)):
                if seq[j] == origseq[j]:
                    mutseq = mutseq + '.'
                else:
                    mutseq = mutseq + seq[j]
                    count = count + 1
            dic[mutseq] = count           
        else:
            continue
        
    dic = sorted(dic, key=dic.get)
    
    file = Sequences.split('S')[0] + 'minimalMut.txt'
    with open(file, 'w') as f:
        for i in range(len(dic)):
            f.write(dic[i]+'\n')
        f.write('\nNumber of sequences: ' + str(num))
    print(num)
    
    