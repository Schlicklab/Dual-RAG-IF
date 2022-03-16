# -*- coding: UTF-8 -*-

# Genetic Alogorithm Inverse Folding 



import random
import os
import os.path
import sys
import time
from functools import partial
import multiprocessing
import statistics
import subprocess


# Chrom contains the points of mutations that can lead to target folding
class Chrom:
    chrom = []
    fitness = 0.0
    folding = ""

    def __init__(self, ngene):
        self.chrom = ['N']*ngene 

    def assign(self):
        lib = ['A','U','C','G']
        ngene = len(self.chrom)
        for i in range(ngene):
            self.chrom[i] = random.choice( lib )

# This function mutates the seq to residues given in chrom, and calculates the fitness by comparing folding to ss_target.
# @ seq is our current sequence with points of mutations written as 'N'
# @ ss_target is in dot format
# @ Nindex is a list of the mutation positions
# @ chrom gives the residues for the mutation positions
# @ k decides the engine of folding to use, we use pknots now for k=1
#                                           use NUPACK for k=2
#                                           use IPknot for k=3
def eachFit(chrom, seq, ss_target, Nindex,k):
    
    fullseq = list(seq)
    for i in range(len(Nindex)):
        fullseq[ Nindex[i] ] = chrom[i]
    fullseq = ''.join(fullseq)

    jobID = str(random.randint(10000,99999))
    jobID = jobID+str(time.time()).split('.')[1]
    list1= ['a','b','c','d','e','f','g','h','i','j']
    list2= [1,2,3,4,5,6,7,8,9,0]
    jobID = jobID + random.choice(list1)+ str(random.choice(list2))
    jobID = jobID + random.choice(list1)+ str(random.choice(list2))
    jobID = jobID + random.choice(list1)+ str(random.choice(list2))
 
    f1 = open("tmpRNAfold"+jobID+".in","w")
    f1.write(">seq\n")
    f1.write(fullseq)
    f1.close()
       
    if k == 1:
       # -k allow pseudoknots, too slow
       os.system("pknots -g tmpRNAfold"+jobID+".in tmpRNAfold"+jobID+".ct 2>/dev/null " )
       #delete the first 4 lines and add "84 MT246482.1 0 13405 84" as the first line for ct file
       with open("tmpRNAfold"+jobID+".ct", "r+") as f:
           lines = f.readlines()
       with open("tmpRNAfold"+jobID+".ct", "w") as f:
           f.write("seq\n")
           for i in range(4, len(lines)):
               f.write(lines[i])
       # export DATAPATH=/Users/qiyaozhu/Downloads/RNAstructure/data_tables/"
       os.system("ct2dot tmpRNAfold"+jobID+".ct 1 tmpRNAfold"+jobID+".out 2>/dev/null" )
       os.system("rm -rf tmpRNAfold"+jobID+".in tmpRNAfold"+jobID+".ct")
            
    elif k == 2:
       os.system("mfe -pseudo -material rna tmpRNAfold"+jobID+" 2>/dev/null ")
       with open("tmpRNAfold"+jobID+".mfe", 'r') as f:
           fold = f.readlines()[16]
       with open("tmpRNAfold"+jobID+".mfe", 'w') as f:
           f.write(">seq\n")
           f.write(fullseq + "\n")
           f.write(fold)
       os.system("dot2ct tmpRNAfold"+jobID+".mfe tmpRNAfold"+jobID+".ct")
       os.system("ct2dot tmpRNAfold"+jobID+".ct 1 tmpRNAfold"+jobID+".out")
       os.system("rm -rf tmpRNAfold"+jobID+".in tmpRNAfold"+jobID+".ct tmpRNAfold"+jobID+".mfe")
       
    elif k == 3:
       os.system("ipknot -g 2 -g 16 -e CONTRAfold -r 1 tmpRNAfold"+jobID+".in > tmpRNAfold"+jobID+"IPknot.txt")
       with open("tmpRNAfold"+jobID+"IPknot.txt", 'r') as f:
           lines = f.readlines()
           for i in range(len(lines)):
               if lines[i][0] == '>':
                   fold = lines[i+2]
       with open("tmpRNAfold"+jobID+"IPknot.txt", 'w') as f:
           f.write(">seq\n")
           f.write(fullseq + "\n")
           f.write(fold)
       os.system("dot2ct tmpRNAfold"+jobID+"IPknot.txt tmpRNAfold"+jobID+".ct")
       os.system("ct2dot tmpRNAfold"+jobID+".ct 1 tmpRNAfold"+jobID+".out")
       os.system("rm -rf tmpRNAfold"+jobID+".in tmpRNAfold"+jobID+".ct tmpRNAfold"+jobID+"IPknot.txt")
       
    passMe = os.path.isfile("tmpRNAfold"+jobID+".out")
    if passMe:
        with open("tmpRNAfold"+jobID+".out") as f:
            lines = f.readlines()
            ss_new = lines[2]
        os.system("rm -rf tmpRNAfold"+jobID+".out")
    else:
        print("no .out\n")
        ss_new = '.'*len(ss_target)

    # get fitness 
    fitness = 0
    for x,y in zip( ss_new, ss_target ):
        if x == y:
            fitness = fitness + 1
    
    print(fitness)
    print(ss_new)
    return fitness, ss_new

# Calculates the fitness for mutation chrom idx in pop
# @ pop is a list of chroms for mutations
def fit4par(idx, pop, seq, ss, Nindex, k):
    pop[idx].fitness, pop[idx].folding = eachFit(pop[idx].chrom, seq, ss, Nindex, k)
    return pop[idx] 

# Calculates the fitness of a population of mutation chroms
# @ nproc is number of CPU processors
def calcFit(pop,seq,ss,Nindex,k,nproc):   
    partial_fit4par = partial(fit4par, pop=pop, seq=seq, ss=ss, Nindex=Nindex, k=k)
    
    if nproc <= 1:
        for idx in range(len(pop)):
            pop[idx] = partial_fit4par(idx)
        return pop    
    else:
        whole = range(len(pop))
        pool = multiprocessing.Pool(nproc)
        result = pool.map( partial_fit4par, whole )
        pool.close()
        pool.join()
        return result

# Calculates the mean fitness of a pop
def meanFit( pop ):
    Fitness = [pop[i].fitness for i in range(len(pop))]
    mean = statistics.mean(Fitness)
    return mean

# Get order of the pop fitness, from highest fitness to lowest
def getOrder( pop ):
    Fitness = [pop[i].fitness for i in range(len(pop))]
    dic = dict(zip(range(len(pop)), Fitness))
    sortPos = sorted(dic, key=dic.get, reverse=True)
    return sortPos

# Finds the position of the mutation chrom with highest fitness
def findBest( pop ):   
    order = getOrder(pop)
    best = [order[0], pop[order[0]].fitness]
    return best

# Replace 50 sequences with the lowest fitness with 50 with the highest fitness
def select( pop, nreplace ):
    
    order = getOrder(pop)
    bestPos = order[0:nreplace]
    worstPos = order[len(pop)-nreplace:len(pop)]

    # replace worst chromosomes with best chromosomes
    for i in range(nreplace):
        pop[ worstPos[i] ].chrom = pop[ bestPos[i] ].chrom[:]
        pop[ worstPos[i] ].fitness = pop[ bestPos[i] ].fitness
        pop[ worstPos[i] ].folding = pop[ bestPos[i] ].folding
                   
    return pop

# Waive the top nwaive chroms in the pop, get list of their positions
def getWaive( pop, nwaive ):
    order = getOrder(pop)
    waiveList = order[0:nwaive]
    return waiveList


# For each pair of sequences, probability of acr of swapping identities of nx random mutable residues
# @ nresi is number of mutable residues
# @ acr is probability of cross-over
# @ nx is number of cross-over residues
def xover1( pop, nresi, acr, nx, waiveList ):

    if nx > 0.5*nresi:
       nx = int(0.1*nresi)
    if nx == 0:
       nx = 1

    for i in range(len(pop)):
        for j in range(len(pop)):
            if i != j:
              if acr > random.random():  
                for k in range(nx):
                   
                   acrNode = random.randrange(0,nresi)
                   if (i not in waiveList) and (j not in waiveList): 
                      tmp = pop[i].chrom[acrNode]
                      pop[i].chrom[acrNode] = pop[j].chrom[acrNode]
                      pop[j].chrom[acrNode] = tmp
                      
                   if (i  in waiveList) and (j not in waiveList): 
                      pop[j].chrom[acrNode] = pop[i].chrom[acrNode]
                      
                   if (i not in waiveList) and (j in waiveList): 
                      pop[i].chrom[acrNode] = pop[j].chrom[acrNode]                   
    return pop

# Generate a mutation residue other than r
def mutIt( r ):
    lib = ['A','U','C','G']    
    lib.remove( r )
    s = random.choice( lib )
    return s 


# Mutate mutable residues
def mutation1( pop, nresi, order, mut_good, mut_ave, mut_bad, waiveList):
    # Fitness level:            - H(20%)    - M(50%)   - L(30%) 
    # Mutate residue number:    - 1         - 2        - 1
    # Mutate probability:       - mut_good  - mut_ave  - mut_bad 
    
    N = len(pop)
    highN = int(0.2*N)
    midN = int(0.5*N)
    
    highList = order[0:highN]
    midList = order[highN:highN+midN]
    lowList = order[highN+midN:N]

    nloop = 1
    for ic in range(len(pop)):
        if ic in highList:
           mut = mut_good
           nloop = 1
        if ic in midList:
           mut = mut_ave
           nloop = 2
        if ic in lowList:
           mut = mut_bad
           nloop = 1
        if ic in waiveList:
           continue
        if mut > random.random():
           # do mutation
           for loop in range(nloop):
               mutNode = random.randrange(0, nresi)
               pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] ) 
    return pop


# Nomination of winning chroms, after being recorded in the heavenList, the chrom has probability of heaven_rate getting 1 residue mutated
# @ if a chromosome has a distance to target less than 'nstepheaven', then pick it out
def heaven( pop, fullfit, nstepheaven, heaven_rate, heavenList ):

    #   give 50% chance to use a lower standard (by 2)  
    #        40%                                 by 4
    #        20%                                 by 6
    #        15%                                 by 8
    # should give more chances to use lower standards (by 4, 6 and 8) 
    
    for ic in range(len(pop)):
        diff = fullfit - pop[ic].fitness

        not_good = False
        if diff > nstepheaven:
           if diff > (nstepheaven+2):
              if diff > (nstepheaven+4):
                 if diff > (nstepheaven+6):
                    if diff > (nstepheaven+8): 
                       not_good = True
                    else: # diff <= (nstepheaven+8)
                       if 0.15 > random.random():
                          not_good = False
                       else:
                          not_good = True 

                 else: # diff <= (nstepheaven+6)
                    if 0.2 > random.random():
                       not_good = False     
                    else:
                       not_good = True

              else: # diff <= (nstepheaven+4)
                 if 0.4 > random.random(): 
                    not_good = False
                 else:
                    not_good = True

           else: # diff <= (nstepheaven+2)
              if 0.5 > random.random():
                 not_good = False
              else:
                 not_good = True 

        if not_good:
           continue
        else:
           # must enter heavenList 
           str1 = ''.join( pop[ic].chrom )
           if [ str1, pop[ic].fitness, pop[ic].folding ] not in heavenList:
              heavenList.append( [ str1, pop[ic].fitness, pop[ic].folding ] )
           # but not necessarily to be mutated    
           if heaven_rate > random.random(): # mutation
              # mutate this chromosome once  
              mutNode = random.randrange(0, len(pop[ic].chrom))
              pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] )

    return pop, heavenList


# Mutate the waiveList chroms if no change in heavenList length for a while
def stuckMutation( pop, waiveList ):
    nprotect = 2
    iprotect = 0
    
    for ic in waiveList:
        iprotect = iprotect + 1  # protect 'nprotect' seq. out of the waiveList from mutation 
        if iprotect <= nprotect:
           continue           
        
        for loop in range(30): 
            mutNode = random.randrange(0, len(pop[ic].chrom) )
            pop[ic].chrom[mutNode] = mutIt( pop[ic].chrom[mutNode] ) 
    return pop


# Get mutation template from seq_tem
# seq_tmp gives suggestions for points of mutation, other points should be the same
def chkSeq(seq, seq_tmp, Nindex):
    
    # check seq length 
    if len(seq) != len(seq_tmp):
       print('Template sequence has inconsistent length.')
       sys.exit()

    # check info
    for i in range(len(seq)):
        if i not in Nindex:
           if seq[i] != seq_tmp[i]:
              print('Template sequence has inconsistent info.')
              sys.exit()
    
    # get initial guess in correct format
    chrom = []
    for i in range(len(Nindex)):
        chrom.append( seq_tmp[Nindex[i]] )

    return chrom



# @ inpf, first line be the target dot format,
#         second line be the actual AUCG sequence we have now and needed to be mutated, points of mutations are written as 'N'
# @ tmpf, template sequence that suggests points of mutations
# @ nseq is population size
# @ nreplace is number of worst chroms being replaced by best chroms
# @ nwaive is waiveList size
# @ niter is number of GA iterations
# @ k is engine, 1 for pknots, 2 for RNAfold
# @ nproc is number of CPU processors
# @ acr_ave is probability of cross-over
# @ nx is number of cross-over residues
# @ mut_good, mut_ave, mut_bad: mutation probabilities for different fitness level
# @ nstepheave, if a chromosome has a distance to target less than 'nstepheaven', then pick it out
# @ heaven_rate, after being nominated in the heavenList, the chrom has probability of heaven_rate getting 1 residue mutated
# @ nprintheaven, print out heaven results every 'nprintheaven' steps if any
# @ nsurvivors is number of nominations
# @ nstill_0, it has been nstill_0 steps that heavenList has no sequence, triggers population re-initialization
# @ nstill_1, it has been nstill_1 steps that heavenList has no change in length, triggers waiveList mutation or population re-initialization
def main(inpf,tmpf,nseq,nreplace,nwaive,niter,k,nproc,acr_ave,nx,mut_good,mut_ave,mut_bad,nstepheaven,heaven_rate,nprintheaven,nsurvivors,nstill_0,nstill_1,design):

   # set up timer 
   timer_a = time.time()
   total_time = 12*60*60 # 30 min is the upper limit of execution 

   # Read in info. of input file
   with open( inpf ) as f:
        ss = f.readline().strip()
        seq = f.readline().strip().upper()

   print("Length of Seq.", len(seq))

   Nindex = [ ] # collect the index of unknown residues 
   for r in range(len(seq)):
       if seq[r] == "N":
          Nindex.append( r )

   # Construct chromosomes
   population = []
   for i in range(nseq):
       population.append( Chrom(len(Nindex)) )
       population[i].assign()

   # Use template [optional] 
   if tmpf != None: # use template 
      with open(tmpf) as f:
         seq_tmp = f.readline().strip()
      for i in range(nwaive): # make the first 'nwaive' chromosomes to take the initial guess
         chrom_tmp = chkSeq(seq, seq_tmp, Nindex)
         chrom_fitness, chrom_folding = eachFit(chrom_tmp, seq, ss, Nindex,k)
         population[i].chrom = chrom_tmp

   # Get initial fitness, waiveList, order
   population = calcFit( population, seq, ss, Nindex, k, nproc )     
   waiveList = getWaive( population, nwaive)
   orderC = getOrder( population )
   heavenList = []


   # Start Iteration
   t = 0
   nochange = 0 # no change in heavenList length
   nboost = 200 # a boost in iteration number t

   while t <= niter:

       # check timer 
       timer_b = time.time()
       if timer_b-timer_a > total_time:
          print("\nTime is up.")
          sys.exit()

       # if not sufficient...
       if t == niter and len(heavenList) < nsurvivors:
          t = 0
          nstepheaven = nstepheaven + 2 # lower the standard for nomination

       t = t + 1
       
       # 1. Mutation
       population = mutation1( population, len(Nindex), orderC, mut_good, mut_ave, mut_bad, waiveList)
       
       # 2. Cross-over
       population = xover1( population, len(Nindex), acr_ave, nx, waiveList)
       population = calcFit( population, seq, ss, Nindex, k, nproc )

       # 3. Selection
       population = select( population, nreplace)
        
       if tmpf != None:
          population[0].chrom = chrom_tmp
          population[0].fitness = chrom_fitness
          population[0].folding = chrom_folding

       waiveList = getWaive( population, nwaive) 
       orderC = getOrder( population )

       # 4. Nomination
       population, heavenList = heaven( population, len(seq), nstepheaven, heaven_rate, heavenList )


       # mutate chromosomes in waiveList if the heavenList does not change for a while (stuck somewhere)
       # or re-initialize the population with a boost in iteration number t
       if t == 1: # first iteration 
          old_heavenlen = len(heavenList)

       if t > 1:
          if old_heavenlen != len(heavenList):
             old_heavenlen = len(heavenList) 
             nochange = 0 # go back to 0 in time, it means 'heavenList' has been changed
             
          else: # no change in length
              
             if old_heavenlen == 0: # nothing in the 'heavenList'
                nochange = nochange + 1 # set up a counter  
                if nochange >= nstill_0: # was set to 150 steps 
                   # re-initialize the whole population (99%)
                   for tmpl in range(1,nseq):
                       if 0.01 < random.random(): 
                          population[tmpl].assign()
                   nochange = 0
                   t = nboost
                   nboost = nboost + 100 # was 100
                   print("\n0 nomination. Population re-initialized", 't =',t, 'nboost =', nboost)


             if old_heavenlen > 0: # there are some seq. in the 'heavenList'
                nochange = nochange + 1
                # if nochange exceeds a limit, then trigger mutation, also set nochange back to 0
                if nochange >= nstill_1: # was set to 100 steps 
                   # mutation
                   if len(heavenList) > 10:
                      population = stuckMutation( population, waiveList )
                      print("\nstuckMutation activated.")
                   else:                     
                      for tmpl in range(1,nseq):
                          if 0.01 < random.random():  
                             population[tmpl].assign()
                      t = nboost
                      nboost = nboost + 100 # was 100
                      print("\npopulation re-initialized", 't =',t, 'nboost =', nboost)
                   nochange = 0
                                             

       # print out heaven results
       if len(heavenList) != 0:
          if (t+1)%nprintheaven == 0: 
             f1 = open(design+'heaven.txt','w')
             
             fullseq = list(seq)
             iseq = ''.join(fullseq)
             f1.write( "Inquiry sequence:\n" + iseq +"\n")
             
             tfold = ''.join(ss)
             f1.write( "Target folding:\n" + tfold +"\n")
             
             if tmpf != None:
                 tseq = ''.join(seq_tmp)
                 f1.write( "Template sequence:\n" + tseq +"\n")
             
             for ih in range(len(heavenList)):
                 newseq = []
                 for j in range(len(fullseq)):
                     if j not in Nindex:
                        newseq.append( fullseq[j] )
                     else:
                        newseq.append( heavenList[ih][0][Nindex.index(j)] )
                 newseq = ''.join(newseq)
                 f1.write(">\n" + newseq)
                 f1.write(" ")
                 f1.write( str(heavenList[ih][1]) + "/" + str(len(fullseq)) +"\n" )
                 f1.write( heavenList[ih][2] +"\n" )
             f1.close()
             pass
             # stop running the program
             if len(heavenList) >= nsurvivors:
                print("\nSufficient chromosomes in heaven:",len(heavenList))
                with open(design+'heaven.txt','a+') as f:
                    f.write("\nNumber of chromosomes in heaven: " + str(len(heavenList)))
                break

   return 0


# Run genetic algorithm
def runGA(inpf, kwargs):
    
   if not os.path.isfile(inpf):
       print("input file not exist...")
       sys.exit()
   
   k = 3 # default
   if 'k' in kwargs:
       k = kwargs['k']
       if k != 1 and k != 2 and k != 3:
           print("engine selection invalid, 1 for PKNOTS, 2 for NUPACK, 3 for IPknot...")
           sys.exit()
   
   tmpf = None 
   if 'tmpf' in kwargs:
       tmpf = kwargs['tmpf']
       if not os.path.isfile(tmpf):
           print("template file not exist...")
           sys.exit()

   design = inpf.split("i")[0]

   
   
   # set up
   nseq = 500 # number of starting candidates (population size)
   nreplace = 50 # kill the last 'nreplace', replace them with the first 'nreplace'
   nwaive = 10 # waive any changes for 'nwaive' best chromosomes  

   nstepheaven = 0 # if a chromosome has a distance to target less than 'nstepheaven', then pick it out
   heaven_rate = 0.75 # probability to be mutated after being lifted to heaven
   nprintheaven = 1 # print out heaven results every 'nprintheaven' steps if any
   nstill_1 = 100 # if the heavenList is not changed for up to 'nstill' iterations, mutate chromosomes in waiveList
   nstill_0 = 150 

   nsurvivors = 500 # 500 # when the heaven has 'nsurvivors' chromosomes, stop the program 

   niter = 500 # number of iteration

   nproc = 4 # number of CPU processors 

   mut_good = 0.30 # probability of mutation
   mut_bad =  0.75
   mut_ave =  0.50


   acr_good = 0.25 # probability of cross-over
   acr_bad =  0.25
   acr_ave =  0.25 
   nx = 4 # number of cross over genes   

   main(inpf,tmpf,nseq,nreplace,nwaive,niter,k,nproc,acr_ave,nx,mut_good,mut_ave,mut_bad,nstepheaven,heaven_rate,nprintheaven,nsurvivors,nstill_0,nstill_1,design)





