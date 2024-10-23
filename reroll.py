#!/usr/local/anaconda/bin/python
import os
import sys
import numpy as np
import pandas as pd
import re
from scipy.stats import gmean

input_file=str(sys.argv[1])
opt_file=str(sys.argv[2])
C2W=str(sys.argv[3])
restriction_list=str(sys.argv[4]).split(" ")
outputdir=str(sys.argv[5])

data1=pd.read_csv(input_file,header=0)
data2=pd.read_csv(opt_file,header=0)

restriction_pattern='|'.join(restriction_list) if len(restriction_list)>0 else ''
restriction_pattern=restriction_pattern.upper()

if len(restriction_pattern)==0:
	print("No restriction_pattern!")
	o=open(outputdir+"/optimized.tsv",'w')
	o.write("ID\tUTR5\toptmized_CDS\tUTR3\taa\tCAI\tcds_GC\tcds_GC3\tcds_T\n")
	for i in range(len(data1)):
		id=str(data1["id"][i])
		utr5=str(data1["utr5"][i])
		utr3=str(data1["utr3"][i])
		cds=str(data2["optimized"][i])
		aa=str(data1["aa"][i])
		cai=data2["CAI"][i]
		gc=data2["GC"][i]
		gc3=data2["GC3"][i]
		t=data2["T"][i]
		o.write(f"{id}\t{utr5}\t{cds}\t{utr3}\t{aa}\t{cai}\t{gc}\t{gc3}\t{t}\n")
	o.close()
	sys.exit()

if C2W.endswith(".csv"):
	df = pd.read_csv(C2W, header=None)
	if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
		df = pd.read_csv(C2W)
		if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U')!= 3:
			raise ValueError("tsv文件内容错误！")
elif C2W.endswith(".tsv"):
	df = pd.read_csv(C2W, sep = "\t", header=None)
	if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
		df = pd.read_csv(C2W, sep="\t")
		if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U')!= 3:
			raise ValueError("tsv文件内容错误！")
else:
	df = pd.read_excel(C2W, header=None)
	if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
		df = pd.read_excel(C2W)
		if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
			raise ValueError("excel文件内容错误！")
C2W={}
for i in range(len(df)):
	codon = df.iloc[i,0].upper().replace('U','T')
	if len(codon) != 3 or codon.count('A')+codon.count('T')+codon.count('C')+codon.count('G') != 3:
		raise ValueError(codon + '不是合法的codon！')
	if codon in C2W:
		raise ValueError(codon+'重复出现！')
	usage = df.iloc[i,1]
	C2W[codon] = usage

codon_dict={
	'G': ['GGT','GGC','GGA','GGG'],
	'A': ['GCT','GCC','GCA','GCG'],
	'V': ['GTT','GTC','GTA','GTG'],
	'L': ['CTT','CTC','CTA','CTG','TTA','TTG'],
	'I': ['ATT','ATC','ATA'],
	'P': ['CCT','CCA','CCG','CCC'],
	'F': ['TTT','TTC'],
	'Y': ['TAT','TAC'],
	'W': ['TGG'],
	'S': ['TCT','TCA','TCC','TCG','AGT','AGC'],
	'T': ['ACT','ACC','ACG','ACA'],
	'C': ['TGT','TGC'],
	'M': ['ATG'],
	'N': ['AAT','AAC'],
	'Q': ['CAA','CAG'],
	'D': ['GAT','GAC'],
	'E': ['GAA','GAG'],
	'K': ['AAA','AAG'],
	'R': ['CGT','CGC','CGG','CGA','AGA','AGG'],
	'H': ['CAT','CAC'],
	'*': ['TAA','TAG','TGA']
}

def custom_exp(x, base):
    return np.exp(x * np.log(base))

def softmax(x,factor=np.exp(1)):
    x_max = np.max(x, axis=-1, keepdims=True)
    e_x = custom_exp(x - x_max,factor)
    return e_x / np.sum(e_x, axis=-1, keepdims=True)

def max_CAI_GC_GC3_minT_v2(aa_seq,C2W=C2W,CAI_strength=101,GC_strength=1,GC3_strength=2,minT_strength=6):
	p_dict={}
	for j in codon_dict.keys():
		p1=np.array(list(map(lambda x:C2W[x],codon_dict[j])))
		p1=softmax(p1,factor=CAI_strength)
		p2=np.array(list(map(lambda x:(x.count('G')+x.count('C')),codon_dict[j])))
		p2=softmax(p2,factor=GC_strength)
		p3=np.array(list(map(lambda x:(x[2].count('G')+x[2].count('C')),codon_dict[j])))
		p3=softmax(p3,factor=GC3_strength)
		p4=np.array(list(map(lambda x:(3-x.count('T')),codon_dict[j])))
		p4=softmax(p4,factor=minT_strength)
		p=p1*p2*p3*p4
		p_dict[j]=p/sum(p)
		if sum(p)==0:
			raise ValueError("密码子%sp值全为0！" %(j))
	res=''.join([np.random.choice(codon_dict[i],p=p_dict[i]) for i in list(aa_seq)])
	return(res)

def CAI(cds_seq, C2W=C2W):
	seq=cds_seq.upper().replace('U','T')
	if len(seq) % 3 != 0:
		raise ValueError("Not a valid coding sequence. Length is not a multiple of 3.")
	w_list = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		if codon not in ['TGG', 'ATG'] and C2W[codon]!= 0:
			w_list.append(C2W[codon])
	return(gmean(w_list))


o=open(outputdir+"/optimized.tsv",'w')
o.write("ID\tUTR5\toptmized_CDS\tUTR3\taa\tCAI\tcds_GC\tcds_GC3\tcds_T\n")
for i in range(len(data1)):
	id=str(data1["id"][i])
	utr5=str(data1["utr5"][i])
	utr3=str(data1["utr3"][i])
	cds=str(data2["optimized"][i])
	aa=str(data1["aa"][i])
	cai=data2["CAI"][i]
	gc=data2["GC"][i]
	gc3=data2["GC3"][i]
	t=data2["T"][i]
	seq=utr5+cds+utr3

	cds_range=(len(utr5),len(utr5)+len(cds))
	match=re.search(restriction_pattern, seq)
	
	k=1
	while match:
		if k>2000:
			raise ValueError("无法roll掉酶切位点！")
		tmp=set(range(*match.span())) & set(range(*cds_range))
		if len(tmp)==0:
			raise ValueError("UTR内存在酶切位点！")
		old_match_span=match.span()
		match_in_cds_range=(min(tmp)-len(utr5),max(tmp)+1-len(utr5))
		alter_codons_range=(int(match_in_cds_range[0]/3),int((match_in_cds_range[1]-1)/3)+1)
		nochange_5_cds=cds[:(3*alter_codons_range[0])]
		nochange_3_cds=cds[(3*alter_codons_range[1]):]
		tmp_aa=aa[alter_codons_range[0]:alter_codons_range[1]]
		kk=1
		while True:
			if kk>1000:
				raise ValueError("无法roll掉酶切位点！")
			new=max_CAI_GC_GC3_minT_v2(tmp_aa)
			new_cds=nochange_5_cds+new+nochange_3_cds
			new_seq=utr5+new_cds+utr3
			match=re.search(restriction_pattern, new_seq)
			if not match:
				cds=new_cds
				seq=new_seq
				break
			if match.span()[0] > old_match_span[0]:
				cds=new_cds
				seq=new_seq
				break
			kk+=1
		k+=1
	cai=CAI(cds)
	gc=(cds.count('G')+cds.count('C'))/len(cds)
	gc3=(cds[2::3].count('G')+cds[2::3].count('C'))/len(cds[2::3])
	t=cds.count('T')/len(cds)
	o.write(f"{id}\t{utr5}\t{cds}\t{utr3}\t{aa}\t{cai}\t{gc}\t{gc3}\t{t}\n")
o.close()
