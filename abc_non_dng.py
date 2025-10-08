#################
   # SET UP #
#################

#   /\_/\  
#  ( o.o ) 
#   > ^ < 

#SET UP
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression

#SETUP


import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import copy
import numpy as np
import sys

from scipy import stats 
from scipy.stats import norm
from scipy.stats import qmc

import itertools

import csv
import pandas as pdszzzzzoewrk,
import gc
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in divide")


#######################
#TOMIMOTO INITIALISE
#######################

# Elongation
# if st_d = num_stem --> Structured Elongation
# if st_d = 0 --> Stochastic Elongation

def mutInStemCells(num_stem, t, mu_0, Genom, st_d):
    # If st_d is a generator or has randomness, ensure it's a float for proportional control
    if isinstance(st_d, np.random.Generator):
        st_d = st_d.random()  # If it's a generator, draw a value
    
    st_d = float(st_d)  # Ensure `st_d` is treated as a float
    
    tCells = []  # List for history of stem cells in meristem
    stemCells = [[i, []] for i in range(num_stem)]  # [(ID number), [list of mutations]]

    tCells = copy.deepcopy(stemCells)

    for k in range(t):  # One cycle of division
        # One daughter cells
        ccells = [[0, 0] for i in range(num_stem)]
        for i in range(num_stem):
            m = random.random()
            if m > mu_0:  # No mutation
                ccells[i] = [stemCells[i][0], stemCells[i][1]]
            else:  # Mutation occurs
                lst = copy.deepcopy(stemCells[i][1])
                lst.append(random.randint(1, Genom))  # Site of Genom in which mutation occurs
                ccells[i] = [stemCells[i][0], lst]

        # The other daughter cells
        ccells2 = [[0, 0] for i in range(num_stem)]
        for i in range(num_stem):
            m = random.random()
            if m > mu_0:  # No mutation
                ccells2[i] = [stemCells[i][0], stemCells[i][1]]
            else:  # Mutation occurs
                lst = copy.deepcopy(stemCells[i][1])
                lst.append(random.randint(1, Genom))  # Site of Genom in which mutation occurs
                ccells2[i] = [stemCells[i][0], lst]

        # Step1 + Step2: choose cells for the next meristem based on `st_d`
        num_structured_cells = int(st_d)  # Integer part: structured selection
        num_random_cells = num_stem - num_structured_cells  # Remaining are selected randomly

        # Select `num_structured_cells` from `ccells2` deterministically (if `st_d > 0`)
        if num_structured_cells > 0:
            structured_cells = ccells2[:num_structured_cells]
        else:
            structured_cells = []

        # Select the remaining cells randomly from both daughter sets (`ccells` and `ccells2`)
        random_cells = random.sample(ccells + ccells2, num_random_cells)

        # Combine structured and random cells to get the stem cells for the next cycle
        stemCells = structured_cells + random_cells
        tCells = tCells + copy.deepcopy(stemCells)

    return tCells

def mutInBrStemCells(num_stem, stemCells, t, mu_0, Genom, st_d):
    st_d = float(st_d)  # Ensure `st_d` is treated as a float

    tCells = copy.deepcopy(stemCells)
    cells = copy.deepcopy(stemCells)  # Prevent reference of stemCells in for loop

    ccells = [[0, 0] for _ in range(num_stem)]
    ccells2 = [[0, 0] for _ in range(num_stem)]

    for k in range(t):  # One cycle of division
        # One daughter cell
        for i in range(num_stem):
            m = random.random()
            if m > mu_0:  # No mutation
                ccells[i] = [cells[i][0], cells[i][1]]
            else:  # Mutation occurs
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1, Genom))  # Site of Genom in which mutation occurs
                ccells[i] = [cells[i][0], lst]

        # The other daughter cell
        for i in range(num_stem):
            m = random.random()
            if m > mu_0:  # No mutation
                ccells2[i] = [cells[i][0], cells[i][1]]
            else:  # Mutation occurs
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1, Genom))  # Site of Genom in which mutation occurs
                ccells2[i] = [cells[i][0], lst]

        # Step1 + Step2: choose cells for the next meristem based on `st_d`
        num_structured_cells = int(st_d)  # Integer part: structured selection
        num_random_cells = num_stem - num_structured_cells  # Remaining are selected randomly

        # Select `num_structured_cells` from `ccells2` deterministically (if `st_d > 0`)
        if num_structured_cells > 0:
            structured_cells = ccells2[:num_structured_cells]
        else:
            structured_cells = []

        # Select the remaining cells randomly from both daughter sets (`ccells` and `ccells2`)
        random_cells = random.sample(ccells + ccells2, num_random_cells)

        # Combine structured and random cells to get the cells for the next cycle
        cells = structured_cells + random_cells
        tCells = tCells + copy.deepcopy(cells)

    return tCells


# Branching
# bias_d = 0.5 ->UNBIASED
# bias_d = 10 -> BIASED

#for small sigma
def WNdist_f(x_value, mu_value, sig_value, n_value): # x, mean, var,

    f_func = 1/(sig_value*math.sqrt(2*math.pi))*np.sum([math.exp(-(x_value + 2*math.pi*k - mu_value)**2/(2*sig_value**2)) for k in range(-n_value,n_value+1)]) #norm.pdf(x, loc, scale)

    return f_func

#for large sigma
def WNdist_g(x_value, mu_value, sig_value, n_value): # x, mean, var,

    g_func = 1/(2*math.pi)*(1 + 2*np.sum([math.exp(-1/2*sig_value**2)**(k**2)*math.cos(k*(x_value-mu_value)) for k in range(1,n_value+1)])) #norm.pdf(x, loc, scale)

    return g_func


def WNdist(x_value, mu_value, sig_value):

    n_value = 5

    if sig_value < 3: #for small sigma
        dist_func = WNdist_f(x_value, mu_value, sig_value, n_value)

    elif sig_value >= 3: #for large sigma
        dist_func = WNdist_g(x_value, mu_value, sig_value, n_value)

    return dist_func


def weightList(num_stem, num_div, bias_d): # Num of stem cells, Num of division, sigma

    cell_dt = 2*math.pi/(num_stem*(2**num_div)) #distance between cells
    mean_num = random.uniform(0,2*math.pi) #center of dist

    wList = np.array([WNdist(i*cell_dt, mean_num, bias_d) for i in range(num_stem*(2**num_div))])

    return wList/np.sum(wList)

#Sampling of stem cells in Branching

def sample_mutations(stem_cells, num_stem, num_div, mu_0, Genom, bias_d): #　cells in apical meristem, Num of stem cells, Num of cell divisions in branching, Mutation rate, Genome size, sigma

    wList = weightList(num_stem, num_div, bias_d)

    cells = copy.deepcopy(stem_cells[-num_stem:])

    for k in range(0,num_div):

        ccells = [[0,0] for i in range(len(cells))]
        ccells2 = [[0,0] for i in range(len(cells))]

        #one daughter cells
        for i in range(len(cells)):
            m = random.random()
            if m > mu_0: #non mutaion
                ccells[i] = [copy.deepcopy(cells[i][0]), copy.deepcopy(cells[i][1])]
            else: #mutation
                lst = ''
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1,Genom)) #site of Genom in which mut occure
                ccells[i] = [copy.deepcopy(cells[i][0]), lst]

        #the other daughter cells
        for i in range(len(cells)):
            m = random.random()
            if m > mu_0: #non mutaion
                ccells2[i] = [copy.deepcopy(cells[i][0]), copy.deepcopy(cells[i][1])]
            else: #mutaionI 
                lst = copy.deepcopy(cells[i][1])
                lst.append(random.randint(1,Genom)) #site of Genom in which mut occure
                ccells2[i] = [copy.deepcopy(cells[i][0]), lst]

        lst2 = [] #list for daughter cells
        for p,q in zip(ccells,ccells2):
            lst2.append(p)
            lst2.append(q)

        cells = copy.deepcopy(lst2) #daughter cells arranged in genealogical order.

        #check
        ##print(cells)

    sampleNum = np.random.choice([i for i in range(len(cells))], size=num_stem, replace=False, p=wList) #sampling of cells for axillary meristem
    sample = [cells[i] for i in sampleNum]

    return sample
def pickUpMut(allmutList):
    mutList = []
    for i in allmutList:
        for j in i[1]:
            mutList.append(j)

    reAllmutList = list(set(mutList)) #eliminate overlapping of mutation
    reAllmutList.sort()

    return reAllmutList

#making mutation matrix from history of meristems
def makeMutMatrix(List, mutList, num_stem): # List_of_meristems、list_of_all_mutated_site、Num_of_stemcells
    preMutMatrix = [[] for i in range(len(List))]
    for i,j in enumerate(List):
        for k in mutList:
            if k in j[1]:
                preMutMatrix[i].append(1)
            else:
                preMutMatrix[i].append(0)

    #check
    for i,j in zip(List, preMutMatrix):
        if len(i[1]) != sum(j):
            print("mutation occured at the same site")
            break

    #sum into meristems
    np_mer_mutMatrix = np.array([0 for i in range(int(len(List)/num_stem))])
    for i in range(int(len(preMutMatrix)/num_stem)): #int(len(preMutMatrix)/num_stem) = number of history of merstems
        mer_mutVector = np.zeros(len(mutList),dtype = int)
        for j in range(num_stem):
            arraylist = np.array(preMutMatrix[i*num_stem +j])
            mer_mutVector = mer_mutVector + arraylist

        if i == 0:
            np_mer_mutMatrix = mer_mutVector
        else:
            np_mer_mutMatrix = np.vstack((np_mer_mutMatrix, mer_mutVector))

    return np_mer_mutMatrix

#null regression
def reg1dim(x, y):
    a = np.dot(x, y)/ (x ** 2).sum()
    return a

#_for saving the data

#for saving chage the shape of list
def add_None_func(ori_List):
    max_length = max([len(i[0]) for i in ori_List]) #scaling of maximum length
    none_List = np.array([[None] for _ in range(len(ori_List[0]))]) #len(ori_List[0]): len(copy_List) due to the type of data

    pre_List = []
    #count = 0 #check
    for j_ in ori_List: #loop by NumTime
        #count += 1 #check
        copy_List = copy.deepcopy(j_)

        while len(copy_List[0]) < max_length: #len(i_[0]): number of mutations
            copy_List = np.hstack((copy_List ,none_List))

        pre_List.append(copy.deepcopy(copy_List))
        #print(count)#check

    return pre_List

#for Refrost data
def non_None(ori_List):
    non_matrix = []
    for i in ori_List:
        pre_matrix = []
        for j in i:
            pre_matrix.append([k for k in j if k is not None])
        non_matrix.append(np.array(copy.deepcopy(pre_matrix)))

    return non_matrix

#######################
#Tree list/dict function
#######################

# Define the configurations for each tree topology
tree_topologies_dict = {
    "test":{
        "numBranch":8,
        "age": 200,
        "s10": 1,
        "b11": 36,"bb11": 130,"b12": 11,"bb12": 103,"b13": 11,"bb13": 98,"b14":93,
        "s40": 36, "b41": 84, "s41": 77, "b42": 80, "s42": 8, "b43": 70, "s44":64
    },
        
    "ubL12": {
        "numBranch": 12,
        "age": 250,
        "s10": 10,
        "b11": 8, "bb11": 16, "b12": 8, "bb12": 16, "b13": 8, "bb13": 16, "b14": 8, "bb14": 16, "b15": 8, "bb15": 16, "b16": 8, "bb16": 16, "b17": 8, "bb17": 16, "b18": 8, "bb18": 16, "b19": 8, "bb19": 16, "b20": 8, "bb20": 16, "b21": 16,
        "s40": 16
    }
}

def create_tree_list_and_dict(tree_topology):
    # Extract the keys from the dictionary
    keys = list(tree_topology.keys())

    # Initialize the tree list
    tree_list = []

    # Initialize the tree dictionary
    tree_dict = {
        "s10": tree_topology["s10"],
        "Ri": [],
        "Rt": [],
        "Li": [],
        "Lt": []
    }

    # Populate the tree list
    for key in keys:
        if key.startswith('s') or key.startswith('b'):
            tree_list.append(tree_topology[key])

    # Find the right-hand branches (before s40)
    s40_index = keys.index('s40')
    right_keys = keys[1:s40_index]  # Exclude s10 and s40
    for i, key in enumerate(right_keys):
        if key.startswith('b') and not key.startswith('bb'):
            if i == len(right_keys) - 1 or not right_keys[i + 1].startswith('bb'):
                tree_dict["Rt"].append(tree_topology[key])  # Last 'b' before s40
            else:
                tree_dict["Ri"].append(tree_topology[key])
        elif key.startswith('bb'):
            tree_dict["Rt"].append(tree_topology[key])

    # Find the left-hand branches (after s40)
    left_keys = keys[s40_index:]  # Include s40
    for i, key in enumerate(left_keys):
        if key.startswith('s') and not key.startswith('b'):
            if i == len(left_keys) - 1 or (i + 1 < len(left_keys) and not left_keys[i + 1].startswith('b')):
                tree_dict["Lt"].append(tree_topology[key])  # Last 'b' on the left side
            else:
                tree_dict["Li"].append(tree_topology[key])
        elif key.startswith('b'):
            tree_dict["Lt"].append(tree_topology[key])

    return tree_list, tree_dict, tree_topology["numBranch"], tree_topology["age"]




def simulate_somatic_mutations(tree_dict, NumStem, NumTime, mu_0, GenSize, StD, NumDiv, nDiv, biasVar):
    job_id = os.getenv('SLURM_ARRAY_TASK_ID', 'default_job')  # Fallback to a default job ID if not set
    unique_filename = f'br_mutmatrix_{job_id}.npy'
    
    NumTime_br_brmutMatrix = []
    #print("mut", mu_0)
    µ_peryear = mu_0 * GenSize
    #print("mut per year", µ_peryear)
    µ14 = (1 - (1 - µ_peryear) ** (1 / nDiv))
    µ13 = (1 - (1 - µ_peryear) ** (1 / nDiv))
    #print("mut13", µ13)

    for num_time in range(NumTime):  # iteration of simulation
        # Simulation begin
        List10 = mutInStemCells(NumStem, tree_dict["s10"], µ14 / 2 + µ13 / 2, GenSize, StD)
        sampleCells10 = sample_mutations(List10, NumStem, NumDiv, µ14, GenSize, biasVar)
        print(List10)

        List_br = []
        List_left = []

        # Right branches (List_br)
        if len(tree_dict["Ri"]) == 0:
            for i in range(len(tree_dict["Rt"])):
                List_br.append(mutInBrStemCells(NumStem, sampleCells10, tree_dict["Rt"][i], µ14, GenSize, StD))
        else:
            for i in range(len(tree_dict["Ri"])):
                if i == 0:
                    List_br.append(mutInBrStemCells(NumStem, sampleCells10, tree_dict["Ri"][i], µ14, GenSize, StD))
                else:
                    List_br.append(mutInBrStemCells(NumStem, List_br[-2][-NumStem:], tree_dict["Ri"][i], µ14, GenSize, StD))

                sampleCells_br = sample_mutations(List_br[-1], NumStem, NumDiv, µ14, GenSize, biasVar)
                List_br.append(mutInBrStemCells(NumStem, sampleCells_br, tree_dict["Rt"][i], µ14, GenSize, StD))

            # Adding the last terminal branches for the last internal branch if more than one terminal branch exists
            if len(List_br) > 1:
                List_br.append(mutInBrStemCells(NumStem, List_br[-2][-NumStem:], tree_dict["Rt"][-1], µ14, GenSize, StD))

        #print(f"Right branches List_br (length): {len(List_br)}")
        sampleCells10_2 = sample_mutations(List10, NumStem, NumDiv, µ13, GenSize, biasVar)

        # Left branches (List_left)
        if len(tree_dict["Li"]) == 0:
            for i in range(len(tree_dict["Lt"])):
                List_left.append(mutInBrStemCells(NumStem, sampleCells10_2, tree_dict["Lt"][i], µ13, GenSize, StD))
        else:
            for i in range(len(tree_dict["Li"])):
                if i == 0:
                    List_left.append(mutInBrStemCells(NumStem, sampleCells10_2, tree_dict["Li"][i], µ13, GenSize, StD))
                else:
                    List_left.append(mutInBrStemCells(NumStem, List_left[-2][-NumStem:], tree_dict["Li"][i], µ13, GenSize, StD))

                sampleCells_left = sample_mutations(List_left[-1], NumStem, NumDiv, µ13, GenSize, biasVar)
                List_left.append(mutInBrStemCells(NumStem, sampleCells_left, tree_dict["Lt"][i], µ13, GenSize, StD))

            # Adding the last terminal branches for the last internal branch if more than one terminal branch exists
            if len(List_left) > 1:
                List_left.append(mutInBrStemCells(NumStem, List_left[-2][-NumStem:], tree_dict["Lt"][-1], µ13, GenSize, StD))

       # print(f"Left branches List_left (length): {len(List_left)}")

        # Combine all lists from branches
        combined_list = [List10] + List_br + List_left
        #print(len(combined_list))

        # Flatten the combined list
        flattened_combined_list = [item for sublist in combined_list for item in sublist]

        # Pass the flattened combined list to pickUpMut
        allMutations = pickUpMut(flattened_combined_list)
        #print("All Mutations (count):", len(allMutations))

        # Ensure we do not access out of range
        if len(List_br) > 0 and len(List_left) > 0:
            endallMutations = pickUpMut([item for sublist in List_br[1::2] for item in sublist] +
                                        [item for sublist in List_left[1::2] for item in sublist] +
                                        List_br[-1] + List_left[-1])
        elif len(List_br) > 0:
            endallMutations = pickUpMut([item for sublist in List_br[1::2] for item in sublist] + List_br[-1])
        elif len(List_left) > 0:
            endallMutations = pickUpMut([item for sublist in List_left[1::2] for item in sublist] + List_left[-1])
        else:
            endallMutations = []

        #print(f"End all mutations count: {len(endallMutations)}")

        # Generate mutation matrices for internal and terminal branches
        ListTo_br_mut_matrices = []
        for i in range(0, len(List_br), 2):
            combined_list = List10 + [item for sublist in List_br[:i + 2] for item in sublist]
            ListTo_br_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        # Ensure the last terminal branch is included if more than one terminal branch exists
        if len(List_br) % 2 == 0 and len(List_br) > 1:
            combined_list = List10 + [item for sublist in List_br for item in sublist]
            ListTo_br_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        ListTo_left_mut_matrices = []
        for i in range(0, len(List_left), 2):
            combined_list = List10 + [item for sublist in List_left[:i + 2] for item in sublist]
            ListTo_left_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        # Ensure the last terminal branch is included if more than one terminal branch exists
        if len(List_left) % 2 == 0 and len(List_left) > 1:
            combined_list = List10 + [item for sublist in List_left for item in sublist]
            ListTo_left_mut_matrices.append(makeMutMatrix(combined_list, allMutations, NumStem))

        pre_br_brmutMatrix = np.array([mat[-1, :] for mat in ListTo_br_mut_matrices + ListTo_left_mut_matrices])
        br_brmutMatrix = np.array([i for i in pre_br_brmutMatrix.T if sum(i) != 0]).T
        #print(f"br_brmutMatrix shape: {br_brmutMatrix.shape}")

        NumTime_br_brmutMatrix.append(br_brmutMatrix)

    #print(f"NumTime_br_brmutMatrix length: {len(NumTime_br_brmutMatrix)}")

    if not NumTime_br_brmutMatrix:
        raise ValueError("NumTime_br_brmutMatrix is empty. Please check the simulation logic.")

    # Save the mutation matrix to the unique file
    try:
        np.save(unique_filename, NumTime_br_brmutMatrix)
        print(f"Successfully saved the mutation matrix to {unique_filename}")
    except Exception as e:
        print(f"Error while saving the file: {e}")

    # Load the mutation matrix from the unique file
    try:
        test_matrixsim = np.load(unique_filename, allow_pickle=True)
        print(f"Loaded mutation matrix shape: {np.array(test_matrixsim[0]).shape}")
    except EOFError as eof_error:
        print(f"EOFError encountered: No data left in file while loading {unique_filename}. Please check the file integrity.")
        raise eof_error  # Re-raise the error for further handling
    except Exception as e:
        print(f"An error occurred while loading the file: {e}")
        raise e

    return test_matrixsim


##############################
#Check for Bugs!
##############################
def is_br_mutmatrix_empty(resultList_2):
    for matrix in resultList_2:
        if len(matrix) == 0 or matrix.shape[0] == 0:
            return True
    return False

##############################
#Count Mutations Function!
##############################
def mut_dist_func(resultList_2):
    try:
        numBranch = len(resultList_2[0])
        #print(f"Number of branches (numBranch): {numBranch}")

        mutShapeTemplate = np.array(list(itertools.product([0, 1], repeat=numBranch)))
        #print(f"mutShapeTemplate shape: {mutShapeTemplate.shape}")

        NumStemList = []
        Nummut_freqTimeFreq = []  # for phylogeny with freq
        Numallmut_freqDifMatrix = []  # for scatter plot with freq

        # The pattern of expansion of somatic mutations
        NumAveallMutDist = []
        NumAveSDallMutDist = []  # for SD

        NumAveallMutFreqDist = []  # freq

        for NumStem in [5]:
            NumTime = len(resultList_2)

            for num_time in range(NumTime):
                try:
                    br_brmutMatrix = resultList_2[num_time]
                    #print(f"br_brmutMatrix shape: {br_brmutMatrix.shape}")

                    # The distribution pattern of somatic mutations
                    # focusing on mutations at the tip of branches
                    mutShapeMatrix = [0 for i in range(2 ** numBranch)]
                    br_brmut01Matrix = np.nan_to_num(br_brmutMatrix / br_brmutMatrix)  # convert to 01 (changing Nan, resulting from 0/0, into 0), row: branch, column: mutation

                    for icut, i in enumerate(mutShapeTemplate):  # mutShapeTemplate is a list for all possible combinations
                        for k in range(len(br_brmutMatrix[0])):
                            if np.all(br_brmut01Matrix[:, k] == i):
                                mutShapeMatrix[icut] += 1

                    if num_time == 0:
                        allMutDist = np.array(mutShapeMatrix)
                        totalMutDist = [mutShapeMatrix]  # for SD
                    else:
                        allMutDist = allMutDist + np.array(mutShapeMatrix)
                        totalMutDist.append(mutShapeMatrix)  # for SD

                    # The distribution pattern of somatic mutations
                    # freq_(枝末端間の変異細胞数の average)
                    mutShapeFreqMatrix = [[0 for i in range(sum(j) * NumStem)] for j in mutShapeTemplate]  # for different Ave_mutations, sum(j)*NumStem: number of average variations

                    for icut, i in enumerate(mutShapeTemplate):  # mutShapeTemplate is a list for all possible combinations
                        for j in br_brmutMatrix.T:
                            if sum(j) != 0 and np.all(np.nan_to_num(j / j) == i):  # eliminate [0,0,...,0] all 0 list
                                mutShapeFreqMatrix[icut][sum(j) - 1] += 1  # sum(j): sum of mutated cells, sum(j)-1: index of list

                    if num_time == 0:
                        allMutFreqDist = copy.deepcopy(mutShapeFreqMatrix)
                    else:
                        for icut, i in enumerate(copy.deepcopy(mutShapeFreqMatrix)):
                            for jcut, j in enumerate(i):
                                allMutFreqDist[icut][jcut] += j

                except IndexError:
                    print(f"IndexError encountered at NumTime {num_time}. Skipping this trial.")
                    continue

            # Distribution of mutations
            AveallMutDist = allMutDist / NumTime  # The pattern of expansion of somatic mutations
            NumAveallMutDist.append(AveallMutDist)  # for different NumStem

            # Standard deviation of dist
            sumtotalMutDist = np.array([float(0) for _ in mutShapeTemplate])
            for k in totalMutDist:
                sumtotalMutDist += (np.array(k) - AveallMutDist) ** 2  # branch-wise
            NumAveSDallMutDist.append((sumtotalMutDist / len(totalMutDist)) ** (1 / 2))  # standard deviation for each branch

            # freq_pattern_Distribution of mutations
            AveallMutFreqDist = [[i / NumTime for i in j] for j in allMutFreqDist]  # The pattern of expansion of somatic mutations with freq
            NumAveallMutFreqDist.append(AveallMutFreqDist)

            NumStemList.append(NumStem)

        # Calculate shared, unique, and total mutations
        shared_mutations = {}
        unique_mutations = [0 for _ in range(numBranch)]
        mutations_counts = [0 for _ in range(numBranch)]

        mutDist_TemFreqList = [[i, j] for i, j in zip(mutShapeTemplate.tolist(), NumAveallMutFreqDist[0])]
        mutDist_TemFreqList_sort_key = sorted(mutDist_TemFreqList, key=lambda x: sum(x[0]))  # sort by template

        mutDist_TemFreqList_sort_key_short = [i for i in mutDist_TemFreqList_sort_key if sum(i[1]) >= 0.1]  # removing mutation with freq < 0.1

        for pattern_freq in mutDist_TemFreqList_sort_key_short:
            pattern, freq_list = pattern_freq[0], pattern_freq[1]
            # Identify shared mutations
            if sum(pattern) > 1:  # Mutation is shared between branches
                branches = tuple([i + 1 for i, x in enumerate(pattern) if x == 1])
                shared_mutations[branches] = shared_mutations.get(branches, 0) + sum(freq_list)

            # Count unique mutations
            if sum(pattern) == 1:  # Mutation is unique to one branch
                branch = pattern.index(1)
                unique_mutations[branch] += sum(freq_list)

            # Count total mutations for each branch
            for i, has_mutation in enumerate(pattern):
                if has_mutation == 1:
                    mutations_counts[i] += sum(freq_list)

        #print(f"Unique mutations: {unique_mutations}")
        #print(f"Mutation counts: {mutations_counts}")

        return NumAveallMutFreqDist, NumAveSDallMutDist, shared_mutations, unique_mutations, mutations_counts

    except IndexError:
        print("IndexError encountered in mut_dist_func. Skipping this function.")
        return [], [], {}, [], []

##############################
#Make dist matrices function
##############################

def gen_matrices(tree_dict, unique_mutations, shared_mutations, mutations_counts):
    numBranch = len(mutations_counts)
    #print(f"Number of branches (numBranch): {numBranch}")

    # Create the genetic distance matrix
    genetic_distance_matrix = np.zeros((numBranch, numBranch))

    # Calculate pairwise genetic distances
    for i in range(numBranch):
        for j in range(numBranch):
            if i == j:
                genetic_distance_matrix[i][j] = 0
            else:
                shared = shared_mutations.get((min(i+1, j+1), max(i+1, j+1)), 0)
                genetic_distance_matrix[i][j] = mutations_counts[i] + mutations_counts[j] - shared

    # Extract branch lengths from tree_dict
    root_length = tree_dict['s10']
    right_internal = tree_dict['Ri']
    right_terminal = tree_dict['Rt']
    left_internal = tree_dict['Li']
    left_terminal = tree_dict['Lt']

    # Calculate root-to-tip distances for the right side
    right_root_to_tip = []
    for i in range(len(right_terminal)):
        distance = root_length + sum(right_internal[:i+1]) + right_terminal[i]
        right_root_to_tip.append(distance)

    # Calculate root-to-tip distances for the left side
    left_root_to_tip = []
    for i in range(len(left_terminal)):
        distance = root_length + sum(left_internal[:i+1]) + left_terminal[i]
        left_root_to_tip.append(distance)

    root_to_tip_distances = right_root_to_tip + left_root_to_tip

    # Create the physical distance matrix
    physical_distance_matrix = np.zeros((numBranch, numBranch))

    # Calculate pairwise physical distances
    for i in range(numBranch):
        for j in range(numBranch):
            if i == j:
                physical_distance_matrix[i][j] = 0
            else:
                physical_distance_matrix[i][j] = abs(root_to_tip_distances[i] - root_to_tip_distances[j])

    # Check if lengths match
    if len(root_to_tip_distances) != numBranch:
        raise ValueError("The length of root_to_tip_distances does not match numBranch")

    #print(f"Root to tip distances: {root_to_tip_distances}")
    #print(f"Length of root_to_tip_distances: {len(root_to_tip_distances)}")

    return genetic_distance_matrix, physical_distance_matrix


#################################
#Calc Variants + mut_rate function
#################################

def calc_variants(genetic_matrix, physical_matrix, tree_dict, GenSize):
    physical_distances = []
    genetic_distances = []

    # Extract the upper triangle of the matrices
    for i in range(len(physical_matrix)):
        for j in range(i + 1, len(physical_matrix)):
            physical_distances.append(physical_matrix[i, j])
            genetic_distances.append(genetic_matrix[i, j])

    physical_distances = np.array(physical_distances).reshape(-1, 1)
    genetic_distances = np.array(genetic_distances)

    # Perform linear regression
    regressor = LinearRegression()
    regressor.fit(physical_distances, genetic_distances)
    predicted_genetic_distances = regressor.predict(physical_distances)

    # Extract branch lengths from tree_dict
    root_length = tree_dict['s10']
    right_internal = tree_dict['Ri']
    right_terminal = tree_dict['Rt']
    left_internal = tree_dict['Li']
    left_terminal = tree_dict['Lt']

    # Calculate root-to-tip distances for the right side
    right_root_to_tip = []
    for i in range(len(right_terminal)):
        distance = root_length + sum(right_internal[:i+1]) + right_terminal[i]
        right_root_to_tip.append(distance)

    # Calculate root-to-tip distances for the left side
    left_root_to_tip = []
    for i in range(len(left_terminal)):
        distance = root_length + sum(left_internal[:i+1]) + left_terminal[i]
        left_root_to_tip.append(distance)

    root_to_tip_distances = right_root_to_tip + left_root_to_tip

    # Find the largest root-to-tip distance
    max_root_to_tip_distance = max(root_to_tip_distances)
    sum_root_tip = sum(root_to_tip_distances)
    
    #TOTAL length sampled in m
    sum_samp = root_length + sum(right_internal) + sum(right_terminal) + sum(left_internal) + sum(left_terminal)
    length_samp = sum_samp/10 


    # Calculate the total variants using the regression formula
    #variants = regressor.predict(np.array([[max_root_to_tip_distance]]))[0]
    variants = regressor.predict(np.array([[length_samp]]))[0]

    # Calculate mutation rate per metre
    output_mut_rate = variants / length_samp
  

    # Return the variants, the regression equation, and the mutation rate per nucleotide
    regression_equation = f"y = {regressor.coef_[0]} * x + {regressor.intercept_}"
    return variants, regression_equation, max_root_to_tip_distance, output_mut_rate


###################
#ABC-ACCEPT/REJECT#
###################

# EPSILON = 20

import numpy as np
import os
import pandas as pd
import time

def sample_prior():
    StD = np.random.uniform(0, 5)
    biasVar = np.random.uniform(0.5, 10)
    input_mut = np.random.uniform(5.0E-11, 1.12E-09)
    return StD, biasVar, input_mut

def run_simulation(StD, biasVar, input_mut):
    GenSize = 500000000
    NumStem = 5
    nDiv = 1
    NumDiv = 7
    NumTime = 1

    tree_list, tree_dict, numBranch, age = create_tree_list_and_dict(tree_topologies_dict['test'])

    resultList_2 = simulate_somatic_mutations(tree_dict, NumStem=int(NumStem), NumTime=int(NumTime),
                                              mu_0=input_mut, GenSize=int(GenSize), 
                                              StD=StD, NumDiv=NumDiv, nDiv=nDiv, biasVar=biasVar)

    NumAveallMutFreqDist, NumAveSDallMutDist, shared_mutations, unique_mutations, mutations_counts = mut_dist_func(resultList_2)
    genetic_matrix, physical_matrix = gen_matrices(tree_dict, unique_mutations, shared_mutations, mutations_counts)
    variants, regression_equation, max_root_to_tip_distance, output_mut_rate = calc_variants(genetic_matrix, physical_matrix, tree_dict, GenSize)

    return unique_mutations, output_mut_rate

def calculate_distance(unique_mutations):
    real_distribution = [37.89349405, 37.99117253, 35.91959651, 42.62965378, 
                         38.59110184, 30.84725743, 32.74895451, 43.37876936]
    
    real_distribution_array = np.array(real_distribution)
    simulated_distribution_array = np.array(unique_mutations)
    distribution_distance = np.linalg.norm(simulated_distribution_array - real_distribution_array)

    return distribution_distance

def save_accepted_samples(accepted_samples, job_id, batch_id, epsilon):
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    accepted_file = os.path.join(output_dir, f"accepted_samples_{job_id}_{batch_id}_{timestamp}_epsilon{epsilon}.csv")
    pd.DataFrame(accepted_samples, columns=["StD", "biasVar", "input_mut", "output_mut_rate"]).to_csv(accepted_file, index=False)
    print(f"Accepted samples saved to {accepted_file}")

def abc_rejection(num_samples, epsilon, job_id, batch_id):
    accepted_samples = []
    save_interval = 100  # Save every 100 iterations

    for i in range(num_samples):
        # Sample from the prior
        StD, biasVar, input_mut = sample_prior()
        unique_mutations, output_mut_rate = run_simulation(StD, biasVar, input_mut)
        
        # Calculate distance and accept/reject
        dist = calculate_distance(unique_mutations)
        if dist < epsilon:
            accepted_samples.append([StD, biasVar, input_mut, output_mut_rate])

        # Save accepted samples periodically
        if i > 0 and i % save_interval == 0:
            save_accepted_samples(accepted_samples, job_id, batch_id, epsilon)
            accepted_samples = []  # Clear the list after saving to avoid duplicates

    # Save any remaining accepted samples after the loop
    if accepted_samples:
        save_accepted_samples(accepted_samples, job_id, batch_id, epsilon)

# Main script to run ABC-Rejection and save results
if __name__ == "__main__":
    num_samples = 5000
    epsilon = 20
    job_id = os.getenv('SLURM_ARRAY_TASK_ID', 'default_job')
    batch_id = sys.argv[1] if len(sys.argv) > 1 else 'default_batch'

    abc_rejection(num_samples, epsilon, job_id, batch_id)