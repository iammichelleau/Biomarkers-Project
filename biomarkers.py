import numpy as np
import os
import json
import matplotlib.pyplot as plt
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()

dir_path = '/Users/iammichelleau/Documents/CSE 283/Final Project/KIRP/TCGA'
filenames = os.listdir(dir_path)
filenames = [filename for filename in filenames if filename[-1] == 't']

def group_exp(samples):
	group = []

	med = np.median(samples)
	for sample in samples:
		if sample < med:
			group.append(1)
		elif sample >= med:
			group.append(2)
	
	return group

def biomart():
	print "Reading epigenetic genes..."

	biomart = []
	f = open('mart_export.txt', 'r')
	line = f.readline()
	line = f.readline()
	while line:
		items = line.split()
		biomart.append(items[0])
		line = f.readline()

	return biomart

def epi_genes(biomart_genes):
	print "Reading TCGA data..."

	epi_genes = []
	f = open(dir_path + '/' + filenames[0], 'r')
	line = f.readline()
	while line:
		items = line.split()
		gene = (items[0].split('.'))[0]
		if gene in biomart_genes:
			epi_genes.append(gene)
		line = f.readline()

	return epi_genes

def exp_values(epi_genes):
	print "Reading expression values..."

	values = np.zeros((len(epi_genes), len(filenames)))
	for j in range(0,len(filenames)): 
		i = 0
		f = open(dir_path + '/' + filenames[j], 'r')
		line = f.readline()
		while line:
			items = line.split()
			gene = (items[0].split('.'))[0]
			value = items[1]
			if gene in epi_genes: 
				values[i][j] = value
				i += 1
			line = f.readline()

	return values

def metadata():
	print "Reading metadata..."
	with open('./KIRP/metadata_KIRP.json') as data_file: 
		data = json.load(data_file)

	T = []
	C = []
	file_names = []
	null_cases = []
	for i in range(0,len(filenames)):
		file_names.append(data[i]["file_name"])

		if data[i]["cases"][0]["diagnoses"][0]["days_to_death"] is not None:
			T.append(data[i]["cases"][0]["diagnoses"][0]["days_to_death"])
			if data[i]["cases"][0]["diagnoses"][0]["vital_status"] == 'dead':
				C.append(True)
			else:
				C.append(False)

		elif data[i]["cases"][0]["diagnoses"][0]["days_to_last_follow_up"] is not None: 
			T.append(data[i]["cases"][0]["diagnoses"][0]["days_to_last_follow_up"])
			if data[i]["cases"][0]["diagnoses"][0]["vital_status"] == 'dead':
				C.append(True)
			else:
				C.append(False)

		if (data[i]["cases"][0]["diagnoses"][0]["days_to_death"] is None) and (data[i]["cases"][0]["diagnoses"][0]["days_to_last_follow_up"] is None): 
			T.append(-1)
			C.append(-1)

	T = [t for (f,t) in sorted(zip(file_names,T))][:]
	C = [c for (f,c) in sorted(zip(file_names,C))][:]
	return T, C

def group_by_exp(values):
	print "Separating genes into groups..."

	groups = []
	# f = open('groups.txt','w')
	for i in range(0,len(epi_genes)):
		group = group_exp(values[i][:])
		groups.append(group)
		# print >> f, group

	return groups

def logrank(T, C):
	print "Running logrank test..."

	p_values = []
	for i in range(0,len(groups)):
		T1 = []
		T2 = []
		C1 = []
		C2 = []
		for j in range(0,len(filenames)):
			if T[j] != -1: 
				if groups[i][j] == 1:
					T1.append(T[j])
					C1.append(C[j])
				if groups[i][j] == 2:
					T2.append(T[j])
					C2.append(C[j])

		if len(T1) != 0:
			logrank = logrank_test(T1, T2, C1, C2, alpha=0.99)		
			p_values.append(logrank.p_value)
		else:
			p_values.append(100)

	return p_values

def kp_plot(min_ids, groups, epi_genes, p_values, T, C):
	for min_id in min_ids: 
		T1 = []
		T2 = []
		C1 = []
		C2 = []

		for j in range(0,len(T)):
			if T[j] != -1: 
				if groups[min_id][j] == 1:
					T1.append(T[j])
					C1.append(C[j])
				if groups[min_id][j] == 2:
					T2.append(T[j])
					C2.append(C[j])

		plt.figure()
		ax = plt.subplot(111)
		kmf.fit(T1, C1, label='Low expression')
		ax = kmf.plot(ax=ax)
		kmf.fit(T2, C2, label='High expression')
		ax = kmf.plot(ax=ax)
		plt.title(epi_genes[min_id] + ' (p-value = ' + str(p_values[min_id]) + ')')
		plt.savefig('./KIRP/Plots/kp' + str(min_id) + '.png')


biomart_genes = biomart()
epi_genes = epi_genes(biomart_genes)
values = exp_values(epi_genes)

T, C = metadata()

groups = group_by_exp(values)
# f = open('groups.txt','r')
# groups = np.loadtxt('groups.txt',delimiter=',')

p_values = logrank(T, C)
min_ids = np.argsort(p_values)[:5]
print "Indices of genes with minimum p-values:", min_ids

kp_plot(min_ids, groups, biomart_genes, p_values, T, C)










