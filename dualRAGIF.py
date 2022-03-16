#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:37:15 2020

@author: qiyaozhu
"""
from copy import *
from dualGraphs import *
from ClassesFunctions import *
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from dualGA import *
from gaif import *
from dualGraphCheck import *
from minimalCount import *
from mutationOptimization import *
from minmutOrganize import *
import sys
import os.path
import re


# from the adjacency matrix of a graph, get all possible vertex orders from 5' to 3' end
# row 0 of adjacency matrix is for vertex 0, etc.
def adjToSequence(AdjM):
    
    Seqs = []
    Mats = []
    successDic = {}
    
    # find one sequence end
    begin = 0    
    while begin < len(AdjM):
        if sum(AdjM[begin]) < 4:
            break
        else:
            begin = begin+1
    Seqs.append(str(begin)+'_')
    Mats.append(deepcopy(AdjM))
    
    # complete every sequence and check if it is possible
    s = 0
    while s < len(Seqs):
        seq = Seqs[s]
        adjC = deepcopy(Mats[s])
        end = False
        old_pos = -1
        
        while not end:
            seq = Seqs[s]
            adjC = deepcopy(Mats[s])
            
            # check if every stem number has been written twice, i.e. end
            end = True
            for n in range(len(adjC)):
                if seq.count(str(n)) < 2:
                    end = False
                        
            if end:
#                print('\nseq %d end check'%s)
#                print(Seqs)
#                print(Mats)
                
                success = True
                for n in range(len(adjC)):
                    if seq.count(str(n)) > 2:
                        success = False
                        successDic[seq] = success
                        break
                # get rid of remaining '_'                       
                if success:
                    pos = seq.find('_')
                    while pos != -1:
                        if pos == len(seq)-1:
                            seq = seq[:pos]
                            pos = seq.find('_')
                        else:
                            n = int(seq[pos-1:pos])
                            m = int(seq[pos+1:pos+2])
                            if adjC[n][m] > 0:
                                seq = seq[:pos] + seq[pos+1:]
                                adjC[n][m] = adjC[n][m] - 1
                                adjC[m][n] = adjC[m][n] - 1
                                Mats[s] = deepcopy(adjC)
                                pos = seq.find('_')
                            else:
                                success = False
                                break
                # check if every edge has been made
                if success:
                    for i in range(len(adjC)):
                        for j in range(len(adjC)):
                            if adjC[i][j] != 0:
                                success = False
                successDic[seq] = success
                s = s+1
                
            else:
                # find first '_' and fulfill it
                pos = seq.find('_')
                if pos == old_pos:
                    end = True
                    successDic[seq] = False
                    s = s+1
                else:
                    old_pos = pos
                
                    x = seq[pos-1]
                    
                    # self-loop, x_ become xx_
                    if adjC[int(x)][int(x)] == 2:
#                        print('\nself-loop')
#                        print('Before')
#                        print(Seqs)
#                        print(Mats)
                        
                        Seqs[s] = seq[:pos] + x + seq[pos:]
                        adjC[int(x)][int(x)] -= 2
                        Mats[s] = deepcopy(adjC)
                        
#                        print('After')
#                        print(Seqs)
#                        print(Mats)
                    
                    # 3, pseudoknot, x_ become xnxn_
                    elif 3 in adjC[int(x)]:
#                        print('\npseudoknot')
#                        print('Before')
#                        print(Seqs)
#                        print(Mats)
                        
                        n = 0
                        while n < len(adjC[int(x)]):
                            if adjC[int(x)][n] == 3:
                                break
                            else:
                                n = n+1
                        if seq.count(x) == 1 and seq.count(str(n)) == 0:
                            Seqs[s] = seq[:pos] + str(n) + x + str(n) + seq[pos:]
                            adjC[int(x)][n] -= 3
                            adjC[n][int(x)] -= 3
                            Mats[s] = deepcopy(adjC)
                        else:
                            end = True
                            successDic[seq] = False
                            s = s+1
                            
#                        print('After')
#                        print(Seqs)
#                        print(Mats)
                    
                    # 2
                    elif 2 in adjC[int(x)]:
#                        print('\n2 in adj')
#                        print('Before')
#                        print(Seqs)
#                        print(Mats)
                        
                        # if x only appears once
                        if seq.count(x) == 1:
                            n = 0
                            while n < len(adjC[int(x)]):
                                if adjC[int(x)][n] == 2:
                                    break
                                else:
                                    n = n+1
                            # xnnx_
                            if adjC[n][n] == 2:
                                if seq.count(str(n)) == 0:
                                    Seqs[s] = seq[:pos] + str(n) + str(n) + x + seq[pos:]
                                    adjC[int(x)][n] -= 2
                                    adjC[n][int(x)] -= 2
                                    adjC[n][n] -= 2
                                    Mats[s] = deepcopy(adjC)
                                else:
                                    end = True
                                    successDic[seq] = False
                                    s = s+1
                            else:
                                # n not written yet, all possibilities: xn_nx_, xn_xn_, x_nxn_, xnx_
                                if seq.count(str(n)) == 0:
                                    Seqs[s] = seq[:pos] + str(n) + '_' + str(n) + x + seq[pos:]                           
                                    Seqs.append(seq[:pos] + str(n) + '_' + x + str(n) + seq[pos:])
                                    Seqs.append(seq[:pos] + '_' + str(n) + x + str(n) + seq[pos:])
                                    Seqs.append(seq[:pos] + str(n) + x + seq[pos:])
                                    adjC[int(x)][n] -= 2
                                    adjC[n][int(x)] -= 2
                                    Mats[s] = deepcopy(adjC)
                                    Mats.append(deepcopy(adjC))
                                    Mats.append(deepcopy(adjC))
                                    Mats.append(deepcopy(adjC))
                                # n written once, it has some other stem around it
                                elif seq.count(str(n)) == 1:
                                    w = seq.find(str(n))
                                    # bring x and n together is separated by '_' only
                                    if w == pos+1:
                                        Seqs[s] = seq[:pos] + seq[pos+1:]
                                        adjC[int(x)][n] -= 1
                                        adjC[n][int(x)] -= 1
                                        Mats[s] = deepcopy(adjC)
                                    # xn_...
                                    else:
                                        Seqs[s] = seq[:pos] + str(n) + seq[pos:]
                                        adjC[int(x)][n] -= 1
                                        adjC[n][int(x)] -= 1
                                        Mats[s] = deepcopy(adjC)
                                # n written twice, they have some other stems around, so forming xnx there not possible
                                elif seq.count(str(n)) == 2:
                                    w1 = seq.find(str(n))
                                    w2 = seq.find(str(n),w1+1)
                                    # nxn forming possible only if n_n exists
                                    if w2-w1 == 2 and seq[w1+1:w1+2] == '_':
                                        Seqs.append(seq[:w1+1] + x + seq[w2:])
                                        tempAdj = deepcopy(adjC)
                                        tempAdj[int(x)][n] -= 2
                                        tempAdj[n][int(x)] -= 2
                                        Mats.append(tempAdj)
                                    # one n needs to be next to x
                                    if w1 == pos+1:
                                        Seqs[s] = seq[:pos] + seq[pos+1:]
                                        adjC[int(x)][n] -= 1
                                        adjC[n][int(x)] -= 1
                                        Mats[s] = deepcopy(adjC)                                   
                                    else:
                                        end = True
                                        successDic[seq] = False
                                        s = s+1
                                else:
                                    end = True
                                    successDic[seq] = False
                                    s = s+1
                                    
                        # x appears twice, two n must be on two available sides of two x        
                        else:
                            p1 = seq.find(x)
                            p2 = seq.find(x,p1+1)
                            n = 0
                            while n < len(adjC[int(x)]):
                                if adjC[int(x)][n] == 2:
                                    break
                                else:
                                    n = n+1
                            if seq.count(str(n)) == 0:
                                # x_x becomes xnnx
                                if adjC[n][n] == 2:
                                    if p2-p1 == 2 and seq[p1+1:p2] == '_':
                                        Seqs[s] = seq[:p1+1] + str(n) + str(n) + seq[p2:]
                                        adjC[int(x)][n] -= 2
                                        adjC[n][int(x)] -= 2
                                        adjC[n][n] -= 2
                                        Mats[s] = deepcopy(adjC)
                                    else:
                                        end = True
                                        successDic[seq] = False
                                        s = s+1
                                # two n must be on two available sides of two x
                                else:
                                    if seq[p2-1:p2] == '_':
                                        w2 = p2-1
                                    elif seq[p2+1:p2+2] == '_':
                                        w2 = p2+1
                                    else:
                                        end = True
                                        successDic[seq] = False
                                        s = s+1
                                        
                                    if not end:
                                        # xn_..._nx
                                        if w2 < p2: 
                                            Seqs[s] = seq[:p1+1] + str(n) + seq[p1+1:w2+1] + str(n) + seq[w2+1:]
                                            adjC[int(x)][n] -= 2
                                            adjC[n][int(x)] -= 2
                                            Mats[s] = deepcopy(adjC)
                                            # xnx
                                            if pos == w2:
                                                Seqs.append(seq[:p1+1] + str(n) + seq[p2:])
                                                Mats.append(deepcopy(adjC))
                                        # xn_...xn_
                                        else:
                                            Seqs[s] = seq[:p1+1] + str(n) + seq[p1+1:p2+1] + str(n) + seq[p2+1:]
                                            adjC[int(x)][n] -= 2
                                            adjC[n][int(x)] -= 2
                                            Mats[s] = deepcopy(adjC)
                                        
                            # x appears twice, two n must be on two available sides of two x, so connect the n with its neighbor x           
                            elif seq.count(str(n)) == 1 or seq.count(str(n)) == 2:
                                w = seq.find(str(n))
                                # n next to first x
                                if w == p1 + 2:
                                    Seqs[s] = seq[:p1+1] + seq[p1+2:]
                                    adjC[int(x)][n] -= 1
                                    adjC[n][int(x)] -= 1
                                    Mats[s] = deepcopy(adjC)
                                    
                                # n before second x
                                elif w == p2-2:
                                    Seqs[s] = seq[:p2-1] + seq[p2:]
                                    adjC[int(x)][n] -= 1
                                    adjC[n][int(x)] -= 1
                                    Mats[s] = deepcopy(adjC)
                                    
                                # n after second x
                                elif w == p2+2:
                                    Seqs[s] = seq[:p2+1] + seq[p2+2:]
                                    adjC[int(x)][n] -= 1
                                    adjC[n][int(x)] -= 1
                                    Mats[s] = deepcopy(adjC)
                                        
                            else:
                                end = True
                                successDic[seq] = False
                                s = s+1
                                                        
#                        print('After')
#                        print(Seqs)
#                        print(Mats)
                    
                    # 1              
                    elif 1 in adjC[int(x)]:
#                        print('\n1 in adj')
#                        print('Before')
#                        print(Seqs)
#                        print(Mats)
                        
                        # find all n that has 1 edge with x
                        l = []
                        n = 0
                        while n < len(adjC[int(x)]):
                            if adjC[int(x)][n] == 1:
                                l.append(n)
                            n = n+1
                            
                        # One of the n must be next to the current x. Change seq for the first possibility, add other possibilities.
                        first = True
                        for n in l:
                            tmpAdj = deepcopy(adjC)
                            
                            if seq.count(str(n)) < 2:
                                if first:
                                    Seqs[s] = seq[:pos] + str(n) + seq[pos:]
                                    tmpAdj[int(x)][n] -= 1
                                    tmpAdj[n][int(x)] -= 1
                                    Mats[s] = deepcopy(tmpAdj)
                                    first = False
                                    if pos < len(seq)-1:
                                        m = int(seq[pos+1])
                                        if n == m:
                                            Seqs.append(seq[:pos] + seq[pos+1:])
                                            Mats.append(deepcopy(tmpAdj))
                                        if tmpAdj[n][m] > 0:
                                            Seqs.append(seq[:pos] + str(n) + seq[pos+1:])
                                            tmpAdj[n][m] -= 1
                                            tmpAdj[m][n] -= 1
                                            Mats.append(deepcopy(tmpAdj))
                                else:
                                    Seqs.append(seq[:pos] + str(n) + seq[pos:])
                                    tmpAdj[int(x)][n] -= 1
                                    tmpAdj[n][int(x)] -= 1
                                    Mats.append(deepcopy(tmpAdj))
                                    if pos < len(seq)-1 and n == int(seq[pos+1]):
                                        Seqs.append(seq[:pos] + seq[pos+1:])
                                        Mats.append(deepcopy(tmpAdj))                               
                                    if pos < len(seq)-1:
                                        m = int(seq[pos+1])
                                        if n == m:
                                            Seqs.append(seq[:pos] + seq[pos+1:])
                                            Mats.append(deepcopy(tmpAdj))
                                        if tmpAdj[n][m] > 0:
                                            Seqs.append(seq[:pos] + str(n) + seq[pos+1:])
                                            tmpAdj[n][m] -= 1
                                            tmpAdj[m][n] -= 1
                                            Mats.append(deepcopy(tmpAdj))
                                            
                            elif seq.count(str(n)) == 2 and pos < len(seq)-1:
                                m = int(seq[pos+1])
                                if n == m:
                                    if first:
                                        Seqs[s] = seq[:pos] + seq[pos+1:]
                                        tmpAdj[int(x)][n] -= 1
                                        tmpAdj[n][int(x)] -= 1
                                        Mats[s] = deepcopy(tmpAdj)
                                        first = False
                                    else:
                                        Seqs.append(seq[:pos] + seq[pos+1:])
                                        tmpAdj[n][m] -= 1
                                        tmpAdj[m][n] -= 1
                                        Mats.append(deepcopy(tmpAdj))
                                
                            else:
                                first = False
                                end = True
                                successDic[seq] = False
                                s = s+1
                        
#                        print('After')
#                        print(Seqs)
#                        print(Mats)
                        
                    else:
#                        print('\n0 in adj')
#                        print(Seqs)
#                        print(Mats)
                        
                        end = True
                        successDic[seq] = False
                        s = s+1
                    
                    
    successSeqs = []
    Dic = {}
    for s in successDic:
        l = []
        for n in list(s):
            if n == '_':
                l.append('_')
            else:
                l.append(str(int(n)+1))
        l = ''.join(l)
        Dic[l] = successDic[s]
        if Dic[l] == True:
            successSeqs.append(l)
            
    return successSeqs, Dic



# sort the vertex orders so that they are in ascending order from 5' to 3' end
def orderSequence(successSeqs):
    orderSeqs = []  
    # number vertices in ascending order from begin to end
    for s in successSeqs:
        seq = list(s)
        orderseq = []
        perm = {}
        pn = 1
        for n in seq:
            if n in perm:
                pass
            else:
                perm[n] = pn
                pn += 1
        for n in seq:
            orderseq.append(str(perm[n]))
        if ''.join(orderseq) not in orderSeqs:
            orderSeqs.append(''.join(orderseq))
    # number vertices in ascending order from end to begin, to get the reversed sequence
    for s in successSeqs:
        seq = list(s)
        orderseq = []
        perm = {}
        pn = 1
        for i in range(len(seq)-1,-1,-1):
            n = seq[i]
            if n in perm:
                pass
            else:
                perm[n] = pn
                pn += 1
        for i in range(len(seq)-1,-1,-1):
            n = seq[i]
            orderseq.append(str(perm[n]))
        if ''.join(orderseq) not in orderSeqs:
            orderSeqs.append(''.join(orderseq))
    return orderSeqs



# get simplified dot-bracket notation, bracket for vertex, may insert dot to represent loops
def sequenceToDB(orderSeqs):
    DBs = []
    for s in orderSeqs:
        seq = s
        DB = []
        perm = {}
        takenList = []
        for n in seq:
            if n in perm:
                pass
            else:
                perm[n] = 'H'
                w1 = seq.find(n)
                w2 = seq.find(n, w1+1)
                for i in range(w1+1, w2):
                    if int(seq[i:i+1]) < int(n):
                        if perm[seq[i]] == 'PK':
                            perm[n] = 'NPK'
                            break
                        else:
                            perm[n] = 'PK'
        for n in seq:
            if perm[n] == 'H':
                if n not in takenList:
                    DB.append('(')
                    takenList.append(n)
                else:
                    DB.append(')')
            elif perm[n] == 'PK':
                if n not in takenList:
                    DB.append('<')
                    takenList.append(n)
                else:
                    DB.append('>')
            else:
                if n not in takenList:
                    DB.append('[')
                    takenList.append(n)
                else:
                    DB.append(']')
#        DBs.append('.'+'.'.join(DB)+'.')
        DBs.append(''.join(DB))
    return DBs
    


def helixOrder(RNA):
    dic = {}
    for i in range(1,len(RNA.Helices)):
        dic[RNA.Helices[i].start] = i
        dic[RNA.Bases[RNA.Helices[i].end].indexBP] = i
    order = sorted(dic.keys())
    seq = []
    for o in order:
        seq.append(str(dic[o]))
    return seq



def ctToSequence(arg):
    RNA = getCTInfo(arg)
    countHelices(RNA)
    changeHelices(RNA)
    order = helixOrder(RNA)
    return RNA, ''.join(order)
    


# Use Smith-Waterman local pairwise sequence alignment to find optimal mutation regions
def mutationRegion(RNA, ori_order, seq):
    mut_region = []
        
    alignments = pairwise2.align.globalms(ori_order, seq, 2, -1, -1, -1)
    align = format_alignment(*alignments[0])
    align = align.split('\n')
    ori_seq = align[0]
    alignment = list(align[1])
    tar_seq = align[2]
    
    for i in range(len(alignment)):
        # mismatch, get all residues in that helix strand
        if alignment[i] == '.':
            n = ori_seq[i]
            first = True
            for j in range(i):
                if ori_seq[j] == n:
                    first = False
            if first:
                start = RNA.Helices[int(n)].start
                end = RNA.Helices[int(n)].end
                for r in range(start, end+1):
                    mut_region.append(r)
            else:
                start = RNA.Bases[RNA.Helices[int(n)].end].indexBP
                end = RNA.Bases[RNA.Helices[int(n)].start].indexBP
                for r in range(start, end+1):
                    mut_region.append(r)
                    
        # gap
        elif alignment[i] == ' ':
            # need to make a stem use loop region if at least 6 residues, or add in neighboring residues
            if ori_seq[i] == '-':                   
                if i == 0:
                    loop = RNA.Helices[1].start
                    if loop > 6:
                       for r in range(1, loop):
                           mut_region.append(r)
                    else:
                        for r in range(1, 7):
                           mut_region.append(r)
                           
                elif ori_seq[i-1] != '-':
                    if i == len(ori_seq)-1:
                        loop = RNA.Bases[RNA.Helices[int(ori_seq[i-1])].start].indexBP
                        if len(RNA.Bases)-1-loop > 6:
                           for r in range(loop+1, len(RNA.Bases)):
                               mut_region.append(r)
                        else:
                            for r in range(len(RNA.Bases)-6, len(RNA.Bases)):
                               mut_region.append(r)
                               
                    else:
                        b = ori_seq[i-1]
                        first = True
                        for j in range(i-1):
                            if ori_seq[j] == b:
                                first = False
                        if first:
                            loop = RNA.Helices[int(b)].end
                        else:
                            loop = RNA.Bases[RNA.Helices[int(b)].start].indexBP    
                        k = loop + 1
                        num = 0
                        while k < len(RNA.Bases) and RNA.Bases[k].helixNumber == 0:
                            mut_region.append(k)
                            num += 1
                            k += 1                              
                        while num < 6:
                            mut_region.append(loop)
                            loop -= 1
                            if k < len(RNA.Bases):
                                mut_region.append(k)
                                k += 1
                                num += 2
                            else:
                                num += 1
                        
            # need to take away a stem strand
            else:
                n = ori_seq[i]
                first = True
                for j in range(i):
                    if ori_seq[j] == n:
                        first = False
                if first:
                    start = RNA.Helices[int(n)].start
                    end = RNA.Helices[int(n)].end
                    for r in range(start, end+1):
                        mut_region.append(r)
                else:
                    start = RNA.Bases[RNA.Helices[int(n)].end].indexBP
                    end = RNA.Bases[RNA.Helices[int(n)].start].indexBP
                    for r in range(start, end+1):
                        mut_region.append(r)
                        
    mut_region = sorted(list(set(mut_region)))
    return mut_region



# Run the entire dual RAG-IF to transform the folding of an RNA sequence into a target graph
# @ arg: a .ct file containing the original sequence's 2D structure in ct format
# @ designM: specifies the design method, 1 for a target dual graph ID, 2 for a target 2D structure file
# @ target: if designM=1, a target dual graph ID; if designM=2, a design file containing a target 2D structure in dot-bracket format and a sequence specifying the mutation regions in 'N'
# @ kwargs: optional arguments in dictionary format, tmpf=a file containing a template sequence, k=folding prediction program (1 for PKNOTS, 2 for NUPACK, 3 for IPknot)
def main(arg, designM, target, kwargs):
    
    # User only gives the target dual graph, need to find good mutation regions
    if designM == '1':
        # get original vertex sequence from the ct file
        RNA, ori_order = ctToSequence(arg)
        
        # get target graph adjacency matrix
        Graphs = []
        n = int(target.split('_')[0])
        eigen_file = "%dEigen"%n
        adj_file = "V%dAdjDG"%n
        loadEigenvalues(Graphs,n,eigen_file)
        loadAdjMatrices(Graphs,n,adj_file)
        AdjM = []
        for g in Graphs:
            if g.graphID == target:
                AdjM = g.adjMatrix
                break
        
        # get all possible target vertex sequences    
        successSeqs, successDic = adjToSequence(AdjM)
        orderSeqs = orderSequence(successSeqs)
        
        # sort the possible vertex sequences by alignment score with original vertex sequence
        # alignment, match score 2, mismatch score -1, gap opening score -1, gap extension score -1
        dic = {}
        for seq in orderSeqs:
           alignments = pairwise2.align.globalms(ori_order, seq, 2, -1, -1, -1)
           dic[seq] = int(alignments[0][2])
           print(format_alignment(*alignments[0]))
        Seqs = [k for k,v in sorted(dic.items(), key=lambda x: x[1], reverse=True)]
        
        # identify mutation regions by alignment, prepare input file for GA
        for i in range(len(Seqs)):
            mut_region = mutationRegion(RNA, ori_order, Seqs[i])        
            with open(target+'_'+str(i+1)+'inpf', 'w') as f:
                f.write(Seqs[i]+'\n')
                ss = []
                for j in range(1,len(RNA.Bases)):
                    if j not in mut_region:
                        ss.append(RNA.Bases[j].nt)
                    else:
                        ss.append('N')
                f.write(''.join(ss))
            print(mut_region)
        
       
            # run dualGA
            runGA_graph(target+'_'+str(i+1)+'inpf', kwargs)
            
            # run dualGraphCheck
            doubleCheck(target+'_'+str(i+1)+'heaven.txt')
            
            # run minimalCount
            os.system("ct2dot "+arg+" 1 "+arg.split('.')[0]+".out")
            minCount(target+'_'+str(i+1)+'Sequences.txt', arg.split('.')[0]+".out")
            
            # run mutationOptimization
            optimization(target+'_'+str(i+1)+'Sequences.txt', arg.split('.')[0]+".out")
            
            # organize optimization results (minmutOrganize)
            minmutOrganize(target+'_'+str(i+1)+'min_mut_analysis', arg.split('.')[0]+".out")
    
    
    # User specifies the target 2D structure and the mutation regions
    elif designM == '2':
        runGA(target, kwargs)                
        
        # run dualGraphCheck
        doubleCheck(target.split('inpf')[0]+'heaven.txt')
        
        # run minimalCount
        os.system("ct2dot "+arg+" 1 "+arg.split('.')[0]+".out")
        minCount(target.split('inpf')[0]+'Sequences.txt', arg.split('.')[0]+".out")
        
        # run mutationOptimization
        optimization(target.split('inpf')[0]+'Sequences.txt', arg.split('.')[0]+".out")
        
        # organize optimization results (minmutOrganize)
        minmutOrganize(target.split('inpf')[0]+'min_mut_analysis', arg.split('.')[0]+".out")
                


# Run the entire dual RAG-IF to transform the folding of an RNA sequence into a target graph
# @ arg: a .ct file containing the original sequence's 2D structure in ct format
# @ designM: specifies the design method, 1 for a target dual graph ID, 2 for a target 2D structure file
# @ target: if designM=1, a target dual graph ID; if designM=2, a design file containing a target 2D structure in dot-bracket format and a sequence specifying the mutation regions in 'N'
# @ kwargs: optional arguments, tmpf=a file containing a template sequence, k=folding prediction program (1 for PKNOTS, 2 for NUPACK, 3 for IPknot)     
if __name__== "__main__":
    
    if len(sys.argv) < 4:
        print("Missing input: \n* a ct file containing the original sequence's 2D structure \n* design method, 1 for a target dual graph ID, 2 for a target 2D structure file \n* if designM=1, a target dual graph ID; if designM=2, a design file containing a target 2D structure in dot-bracket format and a sequence specifying the mutation regions in 'N'")
        sys.exit()
    
    arg = sys.argv[1]
    if not os.path.isfile(arg):
        print("input ct file not exist...")
        sys.exit()
    
    designM = sys.argv[2]   
    target = sys.argv[3]
    
    if designM == '1': # check if the target dual graph is valid
        dualList = []
        with open('dualGraphList.txt', 'r') as f:
            lines = f.readlines()
        for l in lines:
            dualList.append(l.split('\n')[0])
        if target not in dualList:
            print('Please enter a valid target dual graph ID...')
            sys.exit()
            
    elif designM == '2':
        if not os.path.isfile(target):
            print("target design file not exist...")
            sys.exit()    
        
    else:
        print('Please enter a valid design method... 1 for a target dual graph ID, 2 for a target 2D structure file')
        sys.exit()

    
    kwargs = {}
        
    for i in range(4, len(sys.argv)):
        s = sys.argv[i].split('=')
        if len(s) == 2:
            if s[0] == 'tmpf':
                if not os.path.isfile(s[1]):
                    print("template sequence file not exist...")
                    sys.exit()
                else:
                    kwargs['tmpf'] = s[1]
            elif s[0] == 'k':
                if s[1] != '1' and s[1] != '2' and s[1] != '3':
                    print("engine selection invalid, 1 for PKNOTS, 2 for NUPACK, 3 for IPknot...")
                    sys.exit()
                else:
                    kwargs['k'] = int(s[1])
            else:
                print('Invalid optional arguments format: tmpf=a file containing a template sequence, k=folding prediction program (1 for PKNOTS, 2 for NUPACK, 3 for IPknot)')
                sys.exit()
        else:
            print('Invalid optional arguments format: tmpf=a file containing a template sequence, k=folding prediction program (1 for PKNOTS, 2 for NUPACK, 3 for IPknot)')
            sys.exit()

    main(arg, designM, target, kwargs)
  
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                        