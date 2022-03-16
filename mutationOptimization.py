# This script is designed to order the final sequences in:
# 1. heaven.txt
# 2. *.survivors   
# according to the number of point mutations 
#
# also it should print out the mutation pattern 
#
# to be added: automatically remove redundant point mutation(s) and keep only
#              necessary point mutations
#
# * 2018-09-05
# issue: deleting useless mutations one by one has no problem, but sometimes if they are deleted  
#        together, then the secondary structure will be changed. 
#        This implies that there exists interactions between these mutations.  
#        --> if these mutations form a set, then we need to check its all possible subsets and find the largest   
#            set that can keep the secondary structure 
# another issue: above remedy might give a set which is not minimal !  
#
# * 2018-09-06
# issue: when analyzing the results, not all sequences in the '.survivors' are included   
#        --> comment out 
#
#
# issue: try to remove pairs introduced by two mutations, but also test whether the secondary structure 
#        stays the same or not     
#       
#
#
# * 2018-09-11
# issue 1: check inclusion among the unique mutation sets. 
#          If found, only keep the shorter mutation set.    
#
# issue 2: we did not check the topology by RNAfold for those unique mutation sets.
#
#
# * 2018-09-13
# issue: inclusion check not complete 
#
# * 2018-09-20
# issue: the "target_topo" should be specified via the command line with 'sys.argv' 
#
#

import os
import sys
from itertools import chain, combinations
import numpy
from ClassesFunctions import *
from dualGraphs import *


# ------------------- S.J. 03/19/2019 - function to get the NUPACK topology
def getTopo_Nupack(ex_seq):

    f1 = open("tmpNupack.in","w")
    f1.write(ex_seq)
    f1.close()

    os.system("mfe -pseudo -material rna tmpNupack")       
    with open("tmpNupack.mfe", 'r') as f:
        for l in f.readlines():
            if l[0] != '.' and l[0] != '(' and l[0] != ')' and l[0] != '{' and l[0] != '}':
                continue
            else:
                fold = l                  
    with open("tmpNupack.mfe", 'w') as f:
        f.write(">seq\n")
        f.write(ex_seq + "\n")
        f.write(fold)
    # export DATAPATH=/Users/qiyaozhu/Downloads/RNAstructure/data_tables/
    os.system("dot2ct tmpNupack.mfe tmpNupack.ct")
    
    RNA = getCTInfo("tmpNupack.ct")
    os.system("rm -rf tmpNupack.in tmpNupack.ct tmpNupack.mfe")
    
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


# ----------------- S.J. 03/20/2019
def check_inclusion(pattern, unique_pattern): # to check if this pattern includes any of the unique patterns already selected

    for pat in unique_pattern: # for every pattern in unique pattern
        
        flag = True # assume inclusion
        for mut in pat: # for every mutation in this unique pattern
            if mut not in pattern: # if this mutation is not in the pattern
                flag = False
                break

        if flag: # if all mutations in pat are in pattern, therefore this is redundant
            return flag

    # code will come here if none of the unique patterns are included in the pattern
    return False


# -----------------
def getSS(ex_seq): # based on NUPACK 

       f1 = open("tmpNupack.in","w")
       f1.write(ex_seq)
       f1.close()

       os.system("mfe -pseudo -material rna tmpNupack")       
       with open("tmpNupack.mfe", 'r') as f:
           for l in f.readlines():
               if l[0] != '.' and l[0] != '(' and l[0] != ')' and l[0] != '{' and l[0] != '}':
                   continue
               else:
                   fold = l                  
       with open("tmpNupack.mfe", 'w') as f:
           f.write(">seq\n")
           f.write(ex_seq + "\n")
           f.write(fold)
       # export DATAPATH=/Users/qiyaozhu/Downloads/RNAstructure/data_tables/
       os.system("dot2ct tmpNupack.mfe tmpNupack.ct")
       os.system("ct2dot tmpNupack.ct 1 tmpNupack.out")
       with open("tmpNupack.out", 'r') as f:
            target = f.readlines()[2]
       os.system("rm -rf tmpNupack.in tmpNupack.ct tmpNupack.mfe tmpNupack.out")
       
       return target


#--------------
def inputDotBracket(array): # find pairs from dot-bracket notations 

    bracketPositions = []
    bracket2Positions = []
    results = []
    for i, item in enumerate(array):
        if i == 0 and (item == ')' or item == '>'):
            pass
            #print("Non sense ! Exit")
            break

        if item == '(':
            bracketPositions.append(i)
        elif item ==')':
            if len(bracketPositions) > 0:
                openingPosition = bracketPositions.pop()
                results.append([openingPosition+1,i+1 ]  )
                #print(openingPosition+1, '-->', i+1)
                
        elif item == '<':
            bracket2Positions.append(i)
        elif item =='>':
            if len(bracket2Positions) > 0:
                openingPosition = bracket2Positions.pop()
                results.append([openingPosition+1,i+1 ]  )
        else:
            pass
                #print('ERROR: Not a bracket. Word is: %s.' % item)
    return results


#--------------- S.J. 03/20/2019 - changing the name of the parameters
def get_seq(seq1, res_num, seq2 ):
    
    list1 = list(seq1)
    list2 = list(seq2)
    
    #print "Replacing " + str(res_num) + " from " + list1[ res_num-1 ] + " to " + list2[ res_num-1 ]

    list1[ res_num-1 ] = list2[ res_num-1 ]
    test_seq = ''.join(list1)

    return test_seq


#--------------
def all_subsets(ss): # find all subsets except the null set   
    # ss is also a list 
    a = chain(*map(lambda x: combinations(ss, x), range(1, len(ss)+1)))
    list_a = []
    for k in a:
        list_a.append( list(k) )
    
    return list_a


#---------------
# to add topology check[getTopo_RNAfold] - 2018/09/12    
def kick_one( pattern, test_seq, seq_org, ss_current, target_topo ): 

    print("Kicking...")

    istop = 0
    to_kick = []

    for i in range(len(pattern)):
        res_num = int( pattern[i].split('-')[0] )

        new_seq = get_seq( test_seq, res_num, seq_org )
        #SJ commented new_ss = getSS( new_seq )
        new_ss = getTopo_Nupack( new_seq )
        new_topo = getTopo_IPknot( new_seq )

        #SJ commented if new_ss == ss_current and new_topo == target_topo:# modified
        if new_ss == target_topo and new_topo == target_topo:
           to_kick.append( pattern[i] )
           out_seq = new_seq+'.'[:-1]
           break

        if i == (len(pattern)-1):
           out_seq = test_seq+'.'[:-1] 
           istop = 1

    if len(to_kick) != 0:
       pattern.remove( to_kick[0] )

    return pattern, istop, out_seq


#---------------
def getTopo_IPknot( ex_seq ):

    with open("tmpIPknot.in","w") as f:
        f.write(">seq\n")
        f.write(ex_seq)
    
    os.system("ipknot -g 2 -g 16 -e CONTRAfold -r 1 tmpIPknot.in > tmpIPknot.txt")
    with open("tmpIPknot.txt", 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i][0] == '>':
                fold = lines[i+2]
    with open("tmpIPknot.txt", 'w') as f:
        f.write(">seq\n")
        f.write(ex_seq + "\n")
        f.write(fold)
    os.system("dot2ct tmpIPknot.txt tmpIPknot.ct")
    
    RNA = getCTInfo("tmpIPknot.ct")
    os.system("rm -rf tmpIPknot.in tmpIPknot.txt tmpIPknot.ct")
    
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

# -------------------------
# we need to add topology check for RNAfold secondary structure - 2018/09/12
def pat_analyzer( pattern, ss_current, ss_org, seq_sur, seq_org, target_topo, resultfile ):

    pairs_current = inputDotBracket( ss_current ) # find pairs from dot-brackets
    pairs_org     = inputDotBracket( ss_org )

    pattern_2 = []
    info_2 = []

    useful_i = []
    
    print(pattern, "before")
    with open(resultfile, 'a+') as f:
        f.write('[ ')
        for k in pattern:
            f.write(k + ' ')
        f.write('] before\n')


    info_1 = []
    # Analysis 1: find paired partner
    for i in range(len( pattern )):
       
        #print pattern 
        res_num = int( pattern[i].split('-')[0] )

        pair_num_org = 0 # the residue number of its paired partner
        pair_num_current = 0

        for j in range(len(pairs_org)):
            if pairs_org[j][0] == res_num:
               pair_num_org = pairs_org[j][1]
            if pairs_org[j][1] == res_num:
               pair_num_org = pairs_org[j][0]

        for j in range(len(pairs_current)):
            if pairs_current[j][0] == res_num:
               pair_num_current = pairs_current[j][1]
            if pairs_current[j][1] == res_num:
               pair_num_current = pairs_current[j][0]
        print(res_num, pair_num_org, pair_num_current)

        info_1.append( [res_num, pair_num_org, pair_num_current] )
        

    # Analysis 2: try to find the mininal set of mutation 'pattern' 
    #             we need to loop over all mutations in this pattern, remove it then check  
    #             whether the secondary structure will change or not. If SS is not changed, 
    #             then safely remove this point mutation. 
    #             This should repeat until all point mutations will lead to change in SS. 
    #               
    # input: pattern, ss_current, seq_sur, seq_org

    allchange = 0
    min_pattern = pattern[:] # need to be updated
    new_seq_sur = (seq_sur + '.')[:-1] # need to be updated
    while allchange == 0:
          delList = []
          keep_num_list = []
          for i in range(len( min_pattern )):
              res_num = int( min_pattern[i].split('-')[0] )
              #print "Removing the mutation"
              test_seq = get_seq( new_seq_sur, res_num, seq_org ) # delete this point mutation
              # SJ commented - test_ss = getSS( test_seq ) # test the secondary structure afterwards
              test_ss = getTopo_Nupack(test_seq) # S.J. 03/19/2018 - NUPACK topology
              
              test_topo = getTopo_IPknot( test_seq ) # get the RNAfold topology
              # SJ commented if test_ss == ss_current and test_topo == target_topo:
              if test_ss == target_topo and test_topo == target_topo:
                 delList.append( min_pattern[i] )# collects the mutation that can be safely deleted
              else:
                 keep_num_list.append( res_num ) # necessary point mutations
          if len(delList) == 0:
             allchange = 1
             break # break the while loop

          # update the 'new_seq_sur'
          new_seq_sur = (seq_org+'.')[:-1]
          for i in range(len(keep_num_list)):
              #print "Adding back the mutation"
              new_seq_sur = get_seq( new_seq_sur, keep_num_list[i], seq_sur)
          
          # update the 'min_pattern'
          for i in range(len(delList)):
              min_pattern.remove( delList[i] )

    # - - - - - -  - - - - - - -  - - - - - -  - - - -  -  - -

    # double check the secondary structure of 'new_seq_sur' 
    
    print(new_seq_sur)
    for s in range(0,len(seq_org)):
        if seq_org[s] != new_seq_sur[s]:
            print(s+1,seq_org[s],"-",new_seq_sur[s]) # this should be exactly the same as the min pattern being printed next
    
    # SJ commented check_ss = getSS( new_seq_sur )
    check_ss = getTopo_Nupack( new_seq_sur )
    # and topology...
    check_topo = getTopo_IPknot( new_seq_sur )
    print(check_ss, " ", check_topo)

    #SJ commented if check_ss == ss_current and check_topo == target_topo:
    if check_ss == target_topo and check_topo == target_topo:
        print("remove",str(len(pattern)-len(min_pattern)), "point mutations")
        print(min_pattern) #, "b"
        pattern_2 = min_pattern
        #print ""
    else:
        print("try to remove ",str(len(pattern)-len(min_pattern)), "point mutations")
        print(min_pattern)
        print("failed")
        
        # remedy solution -> add the smallest subset of 'deleted' into min_pattern  
        deleted = list( set(pattern) - set(min_pattern) )
        subsets_del = all_subsets( deleted )
        test_min_pattern = min_pattern[:] # added on 2018-09-22
        for i in range(len(subsets_del)):
            test_min_pattern = min_pattern + subsets_del[i]
            test_seq = (seq_org+'.')[:-1]
            for j in range(len(test_min_pattern)):
                res_num = int( test_min_pattern[j].split('-')[0] )
                test_seq = get_seq( test_seq, res_num, seq_sur )
             # SJ commented test_ss = getSS( test_seq )
            test_ss = getTopo_Nupack(test_seq)
            test_topo = getTopo_IPknot( test_seq )
            # SJ commented if test_ss == ss_current and test_topo == target_topo:
            if test_ss == target_topo and test_topo == target_topo:
               print(test_min_pattern) #, "b"
               #pattern_2 = test_min_pattern # wait a minute...
               break # the for loop of 'i'
    
        # to solve another issue... 
        #    input: test_min_pattern 
        #
        #
        istop = 0
        while istop == 0:
            # kick_one has been modified with regard to topology
            # issue: 2018-09-22 
            # 'test_min_pattern' referenced before assignment 
            test_min_pattern, istop, test_seq = kick_one( test_min_pattern, test_seq, seq_org, ss_current, target_topo )

        print(test_min_pattern) #, 'b'
        pattern_2 = test_min_pattern


    # Analysis 3: find pairs introduced by two mutations, and test whether removing them will change  
    #             the secondary structure or not.  
    # 
    # input: info_1, pattern_2, ss_current, ... 
    #
    #
    pairs_int = []
    for i in range(len(info_1)):
        first_n = info_1[i][0]
        third_n = info_1[i][2]
        for j in range(len(info_1)):
            if j > i:
               first_m = info_1[j][0]
               third_m = info_1[j][2]
               if first_n == third_m and first_m == third_n:
                  print("found pair!")
                  pairs_int.append( [first_n, first_m] )

    num_list = []
    for i in range(len(pattern_2)):
        res_num = int( pattern_2[i].split('-')[0] )
        num_list.append( res_num )

    for i in range(len(pairs_int)):
        a = pairs_int[i][0]
        b = pairs_int[i][1]
 
        if a in num_list and b in num_list:
             print("pairs found in pattern_2")
             
             # try to remove these two mutations and test the secondary structure 
             test_seq2 = (seq_org+".")[:-1]
             for j in range(len(pattern_2)):
                 res_num = int( pattern_2[j].split('-')[0] )
                 if ( res_num != a ) and ( res_num != b ):
                    test_seq2 = get_seq( test_seq2, res_num, seq_sur )

             # SJ commented test_ss2 = getSS( test_seq2 )
             test_ss2 = getTopo_Nupack( test_seq2 )

             # get topology 
             test_topo2 = getTopo_IPknot( test_seq2 )

             # SJ commented - if test_ss2 == ss_current and test_topo2 == target_topo:
             if test_ss2 == target_topo and test_topo2 == target_topo:
                to_remove = []
                for k in range(len(pattern_2)):
                    res_num = int( pattern_2[k].split('-')[0] )
                    if res_num == a or res_num == b:
                       to_remove.append( pattern_2[k] )
                for k in to_remove:
                    pattern_2.remove( k )
                #print pattern_2, 'b'
             # update pattern_2
    print(pattern_2, 'b')
    
    with open(resultfile, 'a+') as f:
        f.write("remove " + str(len(pattern)-len(pattern_2)) + " point mutations\n")
        f.write('[ ')
        for k in pattern_2:
            f.write(k + ' ')
        f.write('] after\n')

    return pattern_2



###########
#  MAIN START #
###########
# @ inpf is the same design file provided to genetic algorithm
# @ origf contains the original sequence and its 2D structure in dot-bracket format
def optimization(inpf, origf):
    
    if not os.path.isfile(inpf):
        print("input sequences file not exist...")
        sys.exit()
    
    resultfile = inpf.split('S')[0] + 'min_mut_analysis'
    target_topo = inpf.split('_')[0] + '_' + inpf.split('_')[1].split('S')[0]
    
    with open(origf, 'r') as f:
        lines = f.readlines()
        seq_org = lines[1].split('\n')[0]
        ss_org = lines[2].split('\n')[0]
        
    # order "Sequences.txt" first
    seq_sur = []
    ss_sur = [] # folds by NUPACK
    dif_sur = []
    
    with open(inpf, 'r') as f: # read in sequences in ".survivors" and calculate the number of different residues  
        lines = f.readlines()
        for l in range(len(lines)):
            if lines[l][0] == '>':
                temp_seq = lines[l+1].split('\n')[0]
                count = 0
                seq_sur.append(temp_seq)
                
                temp_ss = lines[l+2].split(' ')[0]
                ss_sur.append(temp_ss)
                
                for x, y in zip( seq_org, temp_seq ):
                    if x != y:
                        count = count + 1
                dif_sur.append( count )
                
            else:
                continue
    
    order = list(numpy.argsort(dif_sur))
    ord_seq_sur = []
    ord_ss_sur = []
    ord_dif_sur = []
    for i in order:
        ord_seq_sur.append(seq_sur[i])
        ord_ss_sur.append(ss_sur[i])
        ord_dif_sur.append(dif_sur[i])
    
    
    # S.J. 03/20/2019 - new code to check patterns and then take unique sequences - based on old code
    unique_pattern = []
    for i in range(len(ord_seq_sur)):
    
        cur_pat = []
        cur_seq = ord_seq_sur[i]
        cur_ss = ord_ss_sur[i]
        for j in range(len(seq_org)):
            if seq_org[j] != cur_seq[j]:
                cur_pat.append( str(j+1)+"-"+cur_seq[j] ) # calculating the pattern
        
        flag = check_inclusion(cur_pat, unique_pattern) # checking if this pattern includes any of the other unique patterns already selected
        if flag: # this includes one of the already selected unique patterns, so this is redundant, so skip this sequence
            continue
        
        with open(resultfile, 'a+') as f:
            f.write("\nOptimizing Sequence " + str(i+1) + '\n')        
        print("\nOptimizing Sequence " + str(i+1))
        # come here if this sequence needs to be analyzed
        opt_pat = []
        opt_pat = pat_analyzer(cur_pat, cur_ss, ss_org, cur_seq, seq_org, target_topo, resultfile)
        opt_pat = sorted(opt_pat)
        unique_pattern.append(opt_pat) # add it to the unique patterns, it will be different from others already selected, otherwise this sequnce would have been skipped anyways
    
    with open(resultfile, 'a+') as f:
        f.write("\nunique patterns: --b\n")
        for k in unique_pattern:
            f.write(str(k) + ' b\n')



