# ADDIT-M: Using multi-class multi-output SVM for imputation with small blocks of missing data
#Authors: Olivia Choudhury and Ankush Chakrabarty
#Last updated: April 18, 2017

import numpy as np
from sklearn import svm
import time, argparse, math

# Parse command-line input
parser = argparse.ArgumentParser(description='ADDIT-M : Impute missing genotype in model species')
# Required arguments
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i','--input', help='Input file containing missing genotype data',required=True)
requiredNamed.add_argument('-tr','--train', help='Training data containing reference genotype panel',required=True)
requiredNamed.add_argument('-o','--output', help='Output file containing imputed genotype data',required=True)
# Optional arguments
parser.add_argument('-w','--window', help='Length of sliding window', default='15')

args = vars(parser.parse_args())


# Initializing SVM parameters
C = 0.001       # Regularization term
gamma = 0.8     # RBF kernel parameter

# Selecting window length for training set
N_wind=int(args['window'])
nn=int(math.floor(N_wind/2))
d = int(0.5*(N_wind-1))

start = time.time()

# Loading data'''
l_missing = open(args['input'],'r').readlines()
missing_matrix = []
for count in range(len(l_missing)):
    temp = l_missing[count].replace('N', '-99')
    l_missing[count] = temp
    temp = temp.split()
    missing_matrix.append(np.array(temp, 'int'))
missing_matrix = np.vstack(missing_matrix)

N_cols = np.shape(missing_matrix)[0]

l_training = open(args['train'], 'r').readlines()
training_matrix = []
for count in range(len(l_training)):
    temp = l_training[count].split()
    training_matrix.append(np.array(temp, 'int'))
training_matrix = np.vstack(training_matrix)

end = time.time()
t_elapsed = end-start
print "Time Lapsed in Loading Data: " + str(t_elapsed)+' sec.'

z = 0
print "Data Loaded! Starting ADDIT-M..."
t = time.time()
imputed_matrix = missing_matrix # for storing imputed values

for i in range(0, np.shape(missing_matrix)[0]):
    
    # Finding locations of missing data '''
    missing_pos = np.where(missing_matrix[i,:] == -99)

    for jj in range(0, len(missing_pos[0])):
        
        k = missing_pos[0][jj]
        # Storing truth for error computation '''
        truth_values = training_matrix[:,k]
            
        # Quick imputing if all truth values are identical
        if len(np.unique(truth_values)) == 1:
            imputed_value = truth_values[0]
        
        # Quick imputing if neighbors are identical (refer to ADDIT-NM) 
        elif missing_matrix[i, k-1] ==  missing_matrix[i, k+1]:
            imputed_value = missing_matrix[i, k-1]
                
        # Imputation via multi-class SVM '''
        else:
            
            missing_vect = missing_pos[0]
            num_found_right = 0
            num_found_left = 0
            search_count_left = 0
            search_count_right = 0
            left_index = np.zeros(nn, 'int')
            right_index = np.zeros(nn, 'int')
            
            # Generating training set '''
            while num_found_right < nn and num_found_left < nn:
                search_count_left += 1
                search_count_right += 1
                if np.any(missing_vect == k + search_count_right) == False:
                    right_index[num_found_right] = k + search_count_right
                    num_found_right += 1
                if np.any(missing_vect == k - search_count_left) == False:
                    left_index[num_found_left] = k - search_count_left
                    num_found_left += 1
                    
            training_set = np.append(training_matrix[:,left_index], training_matrix[:,right_index], axis=1)
            
            # Multiple SVM classifiers for training... '''
            SVMImputer = svm.SVC(kernel='rbf', gamma=gamma, C=C).fit(training_set, truth_values)
            col_missing = l_missing[i].split()
            col_missing = np.array(col_missing, 'int')

            temp = np.append(col_missing[left_index], col_missing[right_index], axis = 0)
            temp = temp.reshape(1,-1)
            
            # Impute via MC-SVM '''
            imputed_value = SVMImputer.predict(temp)
        
        imputed_matrix[i,k] = imputed_value

t_elapsed = time.time() - t
np.savetxt(args['output'], np.uint8(imputed_matrix), delimiter = ' ', fmt='%d')

print('Finished ADDIT-M')
