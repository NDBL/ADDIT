# ADDIT
Imputation tool for missing genotype data in model and non-model species


Please cite: O. Choudhury, A. Chakrabarty, S. Emrich. Highly Accurate and Efficient Data-Driven Methods For Genotype Imputation. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 2017



ADDIT-NM Dependencies:

Python 2.7.5 or higher

Running ADDIT-NM:

usage: python Addit_NM.py [-h] -i INPUT -o OUTPUT [-w WINDOW]
                   [-t SIMILARITY_THRESHOLD] [-sh SIMILARITY_HIGH]
                   [-sl SIMILARITY_LOW] [-fh FREQUENCY_HIGH]
                   [-fl FREQUENCY_LOW]

Example: python Addit_NM.py -i Test_Nonmodel_Missing.txt -o Test_Nonmodel_Imputed.txt
				   
				   
				   
ADDIT-M Dependencies:

NumPy (v1.11.3)
Scikit-Learn (v0.18.1)

Running ADDIT-M:
Test_Model_Training.zip has to be unzipped. It contains the txt file.

usage: python Addit_M.py [-h] -i INPUT -tr TRAIN -o OUTPUT [-w WINDOW]

Example: python Addit_M.py -i Test_Model_Missing.txt -tr Test_Model_Training.txt -o Test_Model_Imputed.txt


Contact: semrich@nd.edu

