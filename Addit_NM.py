
#Authors: Olivia Choudhury and Ankush Chakrabarty
#Last updated: April 18, 2017

import os, fileinput
from difflib import SequenceMatcher
import itertools, argparse, math, time

# Instead of using KNNSearch() to find global similarity, use sliding window
# to find k samples which are similar near the missing data region

def hamming1(str1, str2):
  return sum(itertools.imap(str.__ne__, str1, str2))

#------------------------------------------------------------------------------------
# Parse command-line input
parser = argparse.ArgumentParser(description='ADDIT-NM : Impute missing genotype in non-model species')
# Required arguments
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i','--input', help='Input file containing missing genotype data',required=True)
requiredNamed.add_argument('-o','--output', help='Output file containing imputed genotype data',required=True)
# Optional arguments
parser.add_argument('-w','--window', help='Length of sliding window', default='15')
parser.add_argument('-t','--similarity_threshold', help='Similarity threshold. Select it based on maximum likelihood of similarity scores calculated over all candidate windows.', default='7')
parser.add_argument('-sh','--similarity_high', help='cut-off above which normalized similarity score is considered "high"',default='0.5')
parser.add_argument('-sl','--similarity_low', help='cut-off below which normalized similarity score is considered "low"',default='0.1')
parser.add_argument('-fh','--frequency_high', help='cut-off above which normalized allele frequency is considered "high"',default='0.9')
parser.add_argument('-fl','--frequency_low', help='cut-off below which normalized allele frequency is considered "low"',default='0.05')

args = vars(parser.parse_args())

# calculate param: minn, maxn, allele_freq, k_threshold,

#num_Sample=int(args['sample'])
#num_SNP=int(args['SNP'])
k_threshold=int(args['similarity_threshold'])

# Get size of window
N_window=int(args['window'])
# Calculate range of window on the left
minn=int(math.floor(N_window/2))
# Calculate range of window on the right
maxn=N_window-minn

# Get cut-offs for similarity and allele frequency-based normalized weights
s_high=float(args['similarity_high'])
s_low=float(args['similarity_low'])
f_high=float(args['frequency_high'])
f_low=float(args['frequency_low'])

# Get input and output files
f_input=args['input']
f_out=args['output']

#------------------------------------------------------------------------------------
case1=0
case2=0
case3=0
case4=0

# Load input file
start=time.time()
l_input=open(f_input).readlines()
end=time.time()
t_elapsed=end-start
print "Time Lapsed in Loading Data: " + str(t_elapsed)+' sec.'
f_out=open(f_out,'w')

print "Data Loaded! Starting ADDIT-NM..."

# Calculate number of samples (rows) and SNPs (columns)
num_Sample=len(l_input)
col_SNP=l_input[0].split()
num_SNP=len(col_SNP)

for i in range(num_Sample):

	ln1=l_input[i]
	ln1=ln1.rstrip('\n')

	l1=ln1.split()

	l_allele_parent=l1

	# Find indices where genotype is mising, i.e. has value 'N'
	indices=[j for j, x in enumerate(l1) if x == "N"]
	indices.sort()

	# Create a dictionary where key=range of sliding window in format 'min_index_max' and val=allele in sliding window around missing data
        # Stored index will help to check that area in training data
	dict_window={}

	dict_imputatedallele={}

	# Create a dictionary where key=sample_index and value=similar allele adjacent to index, if any
	dict_adjacentallele={}

	for index in indices:

		window_str=''
		window_range=''

		# Check if missing at the beginning
		if (index<1):
			lineno_1=l_allele_parent[0]

			lineno_2=l_allele_parent[1]
			lineno_2=lineno_2.rstrip('\n')

			lineno_3=l_allele_parent[2]
			lineno_3=lineno_3.rstrip('\n')

			# Check for same allele
			if (lineno_2==lineno_3):

				# Key=Sample#_Index
				keyno=str(i)+'_'+str(index)
				dict_adjacentallele[keyno]=lineno_2
				

		# Check if missing at the end
		if (index>(num_SNP-2)):
			lineno_n2=l_allele_parent[num_SNP-3]
			parent_n2=lineno_n2.rstrip('\n')

			lineno_n1=l_allele_parent[num_SNP-2]
			parent_n1=lineno_n1.rstrip('\n')

			lineno_n=l_allele_parent[num_SNP-1]

			# Check for same allele
			if (parent_n2==parent_n1):
				keyno=str(i)+'_'+str(index)
				dict_adjacentallele[keyno]=parent_n2

		# Check if index lies in between
		if (index>=1 and index<(num_SNP-2)):
			lineno_prev=l_allele_parent[index-1]
			parent_prev=lineno_prev.strip('\n')

			lineno_index=l_allele_parent[index]

			lineno_next=l_allele_parent[index+1]
			parent_next=lineno_next.rstrip('\n')

                        # Check for same allele
                        if (parent_prev==parent_next):
                                keyno=str(i)+'_'+str(index)
				# Save the indices which have identical allele at adjacent positions
				# Create a dictionary, where key=sample ID and position, and value=identical allle
                                dict_adjacentallele[keyno]=parent_prev

		# Get content of window, where window ranges from index-minn to index+maxn (missing genotype is at the center of the window)
		if (index>=minn+1 and index<(num_SNP-maxn+1)):
			for k in range(index-minn,index+maxn): 
				# Remove \n from each allele
				temp=l1[k].rstrip('\n')
				window_str=window_str+temp
			window_range=str(index-minn)+'_'+str(index)+'_'+str(index+maxn)
			
		# Get content of window, where missing genotype is at the beginning
		if (index<minn+1):
			for k1 in range(index+maxn):
				# Remove \n from each allele
                                tmp=l1[k1].rstrip('\n')
                                window_str=window_str+tmp
			window_range=str(0)+'_'+str(index)+'_'+str(index+maxn)

		# Get content of window, where missing genotype is at the end
		if (index>(num_SNP-maxn)):
			for k2 in range(index-minn,num_SNP):
				# Remove \n from each allele
                                tmp1=l1[k2].rstrip('\n')
                                window_str=window_str+tmp1
			window_range=str(index-minn)+'_'+str(index)+'_'+str(num_SNP)

		# Create a dictionary, where key=window range, and value=content of the window
		dict_window[window_range]=window_str

	# For each index, find candidate windows which donot have missing data at the center of the window
	# Find similarity score of these candidate windows based on hamming distance with the query window
	for ky in dict_window.keys():#key=index

		# Get min, index, and max from key
                range_col=ky.split('_')
                range_min=int(range_col[0])
		key=int(range_col[1])
                range_max=int(range_col[2])

		# Create a dictionary where key=Training Sample and value=similarity score between candidate and query windows
		dict_similarity={}

		# QUICK IMPUTE
		# Check if the alleles before and after missing position are identical
		# If true, then assign that values to the missing position
		ln_dict=l_input[i]
		ln_dict=ln_dict.rstrip('\n')
		l_dict=ln_dict.split()		
		
		key_check=str(i)+'_'+str(key)

		# Impute if identical adjacent allele
		if key_check in dict_adjacentallele.keys():
			parent_key=dict_adjacentallele[key_check]
			if (parent_key=='1'):
				allele_dict='1'

			else:
                                allele_dict='0'

			dict_imputatedallele[key]=allele_dict

			continue

		dict_hammdist={}

		# Iter over number of training samples 
		for p in range(num_Sample):
		
			# Skip for the query sample
			if p==i:
				continue

			else:
				ln2=l_input[p]
				l2=ln2.split()
			
				train_str=''
				train_index=l2[key].rstrip('\n')

				# Ignore if training data has a missing genotype at the position of 'index'
				# This also ignores training using the same file
				if train_index=='N':
					continue
				else:
					for q in range(range_min,range_max):
					# Remove \n from each allele
						temp1=l2[q].rstrip('\n')
						train_str=train_str+temp1

					# Calculate similarity score based on hamming distance
					similarity_percent=len(train_str)-hamming1(dict_window[ky],train_str)

				# Create a dictionary for each query window, where key=training sample, and value=similarity score wih respect to query window
				dict_similarity[p]=similarity_percent

		# Sort dictionary in descending order of value (similarity score)
		# This dictionary contains the list of training samples in the order of their similarity score
		dict_sorted=sorted(dict_similarity.items(), key=lambda x:x[1],reverse=True)
		
		# Pick any number of similar training samples with similarity score > similarity threshold
		out_dict=dict_sorted

		list_same_kthreshold=[]
		out_dict_kthreshold=[]
		list_same=[]
		new_k_threshold=k_threshold

		# Select adaptive number of trusted candidates based on similarity threshold
		while (len(list_same)==0):
			for elmnt in out_dict:
				elmnt=str(elmnt)
				col1_elmnt=elmnt.split(', ')
				tmp1=col1_elmnt[0]
				sample_threshold=int(tmp1.split('(')[1])
				tmp=col1_elmnt[1]
				col2_elmnt=tmp.split(')')
				hamm_dist=int(col2_elmnt[0])

				# Check threshold
				if (hamm_dist>new_k_threshold):
					out_tuple=(sample_threshold,hamm_dist)
					list_same_kthreshold.append(out_tuple)

			for iter1 in range(len(list_same_kthreshold)):
				list_same.append(list_same_kthreshold[iter1][0])

		# Check if list is empty for a given threshold
			new_k_threshold=new_k_threshold-1
				

		dict_allele_count={}
		dict_weight={}
		dict_allele_weight={}
		dict_allele_avgweight={}
		dict_allele_normavg={}
		dict_allele_avgfreq={}

		for ci in list_same:

			indx=int(ci)
			ln_NN=l_input[indx]
			line_NN=ln_NN.split()

			l_NN=line_NN[key]
			l_NN=l_NN.rstrip('\n')

			# Count allele frequency of candidate windows 
			if l_NN in dict_allele_count:
				dict_allele_count[l_NN]=dict_allele_count[l_NN]+1
				dict_allele_weight[l_NN]=dict_allele_weight[l_NN]+dict_similarity[ci]

			else:
				dict_allele_count[l_NN]=1
				dict_allele_weight[l_NN]=dict_similarity[ci]

		#Assign weight
		l_avg=[]
		for ci_avg in dict_allele_count:
			dict_allele_avgweight[ci_avg]=(dict_allele_weight[ci_avg]/float(dict_allele_count[ci_avg]))/float(N_window-1)
			tmp_wght=dict_allele_avgweight[ci_avg]
			l_avg.append(ci_avg)

		#Assign normalized weight based on similarity score
		for ci_normavg in dict_allele_count.keys():
			if (len(l_avg)>1):
				dict_allele_normavg[ci_normavg]=dict_allele_avgweight[ci_normavg]/float(dict_allele_avgweight[l_avg[0]]+dict_allele_avgweight[l_avg[1]])

			else:
				dict_allele_normavg[ci_normavg]=dict_allele_avgweight[ci_normavg]

		# Assign normalized weight based on allele frequency 
		for ci_freq in dict_allele_count.keys():
			if (len(l_avg)>1):
				dict_allele_avgfreq[ci_freq]=format(dict_allele_count[ci_freq]/float(dict_allele_count[l_avg[0]]+dict_allele_count[l_avg[1]]),'f')
			else:
				dict_allele_avgfreq[ci_freq]=format(dict_allele_count[ci_freq]/float(dict_allele_count[l_avg[0]]),'f')

		# Sort dict_allele_weight based on allele weight when each allele has the same frequency
		dict_allele_weight_sorted=sorted(dict_allele_weight.items(), key=lambda x:x[1],reverse=True)
		dict_sortedkey={}
		dict_sortedkey=sorted(dict_allele_count.items(), key=lambda x:x[0])
		
		# Determine imputed allele based on combination of normalized weights
		l_avg = sorted(l_avg)
		for chk in l_avg:

			s=float(dict_allele_normavg[chk])
			f=float(dict_allele_avgfreq[chk])

			# Check if both wights are "high"
			if (s>=s_high):
				imputed_allele=chk
				dict_imputatedallele[key]=imputed_allele
				case1+=1

			elif ((s_low<=s<s_high) and (f>=f_low)):
				imputed_allele=chk
				dict_imputatedallele[key]=imputed_allele
				case2+=1

			elif ((0<s<s_low) and (f>=f_low)):
				imputed_allele=chk
				dict_imputatedallele[key]=imputed_allele
				case3+=1

			else:
				case4+=1
				if max(s,f)<0.1:
					imputed_allele=chk
					dict_imputatedallele[key]=imputed_allele
				else:
					imputed_allele=l_avg[1]
					dict_imputatedallele[key]=imputed_allele


	# Open input file and substitute missing with imputed genotypes
	ln_missing=l_input[i]
	line_missing=ln_missing.split()

	res=''
	cmiss=0

	for ctr in range(num_SNP):
		if (ctr in indices):
			if (ctr in dict_imputatedallele.keys()):
				res=res+dict_imputatedallele[ctr]+' '
			else:
				cmiss+=1
				res=res+'0 '
		else:
			res=res+line_missing[ctr]+' '

	res=res+'\n'
	f_out.write(res)

f_out.close()


print 'Finished ADDIT-NM'

