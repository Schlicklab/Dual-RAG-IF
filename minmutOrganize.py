#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 01:04:19 2020

@author: qiyaozhu
"""

import re


# @ Sequences file contains all optimized mutation candidates
# @ origf file contains the original sequence
# Organize the optimized mutations in ascending order
def minmutOrganize(Sequences, origf):
    
    dic = {}
    
    with open(origf, 'r') as f:
        f.readline()
        tmpseq = f.readline().split('\n')[0]  
        tmpseq = list(tmpseq)
    
    muts = []
    with open(Sequences, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i][0] != 'u':
                continue
            else:
                muts = lines[i+1:len(lines)]
                break
    
    for mut in muts:
        MUT = '['
        m = re.findall('\d+-\w', mut)
        for i in range(len(m)):
            x = m[i]
            num = x.split('-')[0]
            res = x.split('-')[1]
            ori = tmpseq[int(num)-1]
            if i != len(m)-1:
                MUT = MUT + num + ori + '-' + res + ', '
            else:
                MUT = MUT + num + ori + '-' + res + ']\n'
        dic[MUT] = len(m)
    
    dic = sorted(dic, key=dic.get)
    with open(Sequences.split('m')[0]+'min_mut', 'w') as f:
        for i in range(len(dic)):
            f.write(dic[i])  
    
    