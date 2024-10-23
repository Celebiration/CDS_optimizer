#!/usr/local/anaconda/bin/python
import sys
import numpy as np
import pandas as pd
import argparse
import re
from scipy.stats import gmean
sys.path.append('/home/fengchr/software/')
sys.path.append('/home/fengchr/software/arnie')
from arnie.bpps import bpps
from arnie.mfe import mfe
from arnie.free_energy import free_energy
from arnie.pfunc import pfunc
import concurrent.futures

# 创建ArgumentParser对象
parser = argparse.ArgumentParser(description="version 2.0: 增加了提供UTR情况下的AUP14优化以及序列报告功能。\n\n本程序用于根据密码子使用频率联合优化CDS序列的CAI/GC/GC3/T。")

# 添加关键字参数
parser.add_argument('-i','--input', dest='input_fa',type=str, help='需要优化的CDS序列，fasta格式', required=False)
parser.add_argument('-o','--output', dest='output', type=str, help='结果路径', required=False)
parser.add_argument('--roll', action="store_true", help='若指定，则输入为DNA序列，且仅需roll掉限制位点。', required=False)
parser.add_argument('--report', dest='report', help='需要报告的DNA序列。包含id, utr5, cds, utr3, restriction sites 5列 (无header)。', required=False)
parser.add_argument('--C2W', type=str, help='指定使用的C2W表文件路径（xlsx/tsv），若不指定，则默认使用human codon', required=False)
parser.add_argument('--CAI_s', type=float, help='指定CAI的优化强度（default: 2）', required=False, default=2)
parser.add_argument('--GC3_s', type=float, help='指定GC3的优化强度（default: 1）', required=False, default=1)
parser.add_argument('--GC_s', type=float, help='指定GC的优化强度（default: 1）', required=False, default=1)
parser.add_argument('--minT_s', type=float, help='指定minT的优化强度（default: 10）', required=False, default=10)
# parser.add_argument('--CAI_p', type=float, help='指定CAI的优化占比（default: 4）', required=False, default=4)
# parser.add_argument('--GC_p', type=float, help='指定GC的优化占比（default: 1）', required=False, default=1)
# parser.add_argument('--GC3_p', type=float, help='指定GC3的优化占比（default: 2）', required=False, default=2)
# parser.add_argument('--minT_p', type=float, help='指定minT的优化占比（default: 1）', required=False, default=1)
parser.add_argument('-r','--restriction_sites', dest='restriction_sites', type=str, help='指定要避免的位点，可指定多个', required=False, default='', nargs='+')
parser.add_argument('--utr5', dest='utr5', type=str, help='若指定5\' UTR，则优化AUP14', required=False)
parser.add_argument('--utr3', dest='utr3', type=str, help='若不优化AUP14，则不需要指定', required=False)
parser.add_argument('-d', '--AUP14_opt_times', dest='times', type=int, help='优化AUP14的随机次数', required=False, default=20)
parser.add_argument('--rep', dest='rep', type=int, help='每条重复优化次数', required=False, default=1)

# 解析参数
args = parser.parse_args()
if args.input_fa:
	input_fa=args.input_fa
else:
	input_fa='tmp.fa'

if args.output:
	output=args.output
elif args.report:
	output=args.report+'.report.csv'
else:
	output=input_fa+'.opt.csv'

use_linear=True

def read_fasta(input):
	with open(input,'r') as f:
		fasta = {}
		for line in f:
			line = line.strip()
			if line[0] == '>':
				header = line[1:]
			elif len(line)==0:
				continue
			else:
				sequence = line.upper()
				fasta[header] = fasta.get(header,'') + sequence
	return fasta

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

#from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
# C2W_human = {'TTT': 0.85185, 'TTC': 1.0, 'TTA': 0.2, 'TTG': 0.325, 'TCT': 0.79167,
# 	   'TCC': 0.91667, 'TCA': 0.625, 'TCG': 0.20833, 'TAT': 0.78571, 'TAC': 1.0,
# 	   'TAA': 0.6383, 'TAG': 0.51064, 'TGT': 0.85185, 'TGC': 1.0, 'TGA': 1.0,
# 	   'TGG': 1.0, 'CTT': 0.325, 'CTC': 0.5, 'CTA': 0.175, 'CTG': 1.0,
# 	   'CCT': 0.90625, 'CCC': 1.0, 'CCA': 0.875, 'CCG': 0.34375, 'CAT': 0.72414,
# 	   'CAC': 1.0, 'CAA': 0.36986, 'CAG': 1.0, 'CGT': 0.38095, 'CGC': 0.85714,
# 	   'CGA': 0.52381, 'CGG': 0.95238, 'ATT': 0.76596, 'ATC': 1.0, 'ATA': 0.3617,
# 	   'ATG': 1.0, 'ACT': 0.69444, 'ACC': 1.0, 'ACA': 0.77778, 'ACG': 0.30556,
# 	   'AAT': 0.88679, 'AAC': 1.0, 'AAA': 0.75439, 'AAG': 1.0, 'AGT': 0.625,
# 	   'AGC': 1.0, 'AGA': 1.0, 'AGG': 1.0, 'GTT': 0.3913, 'GTC': 0.52174,
# 	   'GTA': 0.26087, 'GTG': 1.0, 'GCT': 0.675, 'GCC': 1.0, 'GCA': 0.575,
# 	   'GCG': 0.275, 'GAT': 0.85185, 'GAC': 1.0, 'GAA': 0.72414, 'GAG': 1.0,
# 	   'GGT': 0.47059, 'GGC': 1.0, 'GGA': 0.73529, 'GGG': 0.73529}

#from compodynamics
C2W_human = {
	'GGT': 0.54376, 'GGC': 1.0, 'GGA': 0.86063, 'GGG': 0.7728, 'GCT': 0.72597, 'GCC': 1.0,
	'GCA': 0.65671, 'GCG': 0.2328, 'GTT': 0.44954, 'GTC': 0.52234, 'GTA': 0.29523, 'GTG': 1.0,
	'CTT': 0.38702, 'CTC': 0.49271, 'CTA': 0.20392, 'CTG': 1.0, 'TTA': 0.23835, 'TTG': 0.36726,
	'ATT': 0.8763, 'ATC': 1.0, 'ATA': 0.43002, 'CCT': 1.0, 'CCA': 0.97973, 'CCG': 0.32428,
	'CCC': 0.99933, 'TTT': 0.97134, 'TTC': 1.0, 'TAT': 0.89147, 'TAC': 1.0, 'TGG': 1.0,
	'TCT': 0.85459, 'TCA': 0.71603, 'TCC': 0.88088, 'TCG': 0.20751, 'AGT': 0.7073, 'AGC': 1.0,
	'ACT': 0.7936, 'ACC': 1.0, 'ACG': 0.31518, 'ACA': 0.92128, 'TGT': 0.96879, 'TGC': 1.0,
	'ATG': 1.0, 'AAT': 0.99982, 'AAC': 1.0, 'CAA': 0.39619, 'CAG': 1.0, 'GAT': 0.97872,
	'GAC': 1.0, 'GAA': 0.84468, 'GAG': 1.0, 'AAA': 0.86642, 'AAG': 1.0, 'CGT': 0.344,
	'CGC': 0.66507, 'CGG': 0.80867, 'CGA': 0.48188, 'AGA': 1.0, 'AGG': 0.91389, 'CAT': 0.8084,
	'CAC': 1.0, 'TAA': 0.0, 'TAG': 0.0, 'TGA': 0.0
}

dict2={}
aminos=list(codon_dict.keys())
codons=list(codon_dict.values())
for i in range(len(codon_dict)):
	for j in codons[i]:
		dict2[j]=aminos[i]

def translate(seq):
	seq=seq.upper().replace('U','T')
	if len(seq) % 3 !=0:
		raise ValueError('''cds length isn't multiple of 3''')
	else:
		seq1=[seq[(3*i):(3*i+3)] for i in range((len(seq)+2)//3)]
		return(''.join(list(map(lambda x:dict2[x],seq1))))

def rev_translate(aa_seq):
	return(''.join([np.random.choice(codon_dict[i]) for i in list(aa_seq)]))

def cai(cds_seq, C2W):
	seq=cds_seq.upper().replace('U','T')
	if len(seq) % 3 != 0:
		raise ValueError("Not a valid coding sequence. Length is not a multiple of 3.")
	w_list = []
	# p=[]
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		# p+=3*[1-C2W[codon]]
		# Do not count W or M codon since there is only one that encodes them
		if codon not in ['TGG', 'ATG'] and C2W[codon]!= 0:
			w_list.append(C2W[codon])
	return(gmean(w_list))

def AUP14(seq):
	seq=seq.upper().replace('U','T')
	UP_vec=1-np.sum(bpps(seq, package='vienna',linear=use_linear),axis=0)
	return(sum(UP_vec[:14])/14)

def custom_exp(x, base):
    return np.exp(x * np.log(base))

def softmax(x,factor=np.exp(1)):
    """
    Compute the softmax of vector x.

    Parameters:
    x (numpy.ndarray): Input vector or matrix.

    Returns:
    numpy.ndarray: Softmax applied to the input.
    """
    # Subtract the max value from each element for numerical stability
    x_max = np.max(x, axis=-1, keepdims=True)
    e_x = custom_exp(x - x_max,factor)
    return e_x / np.sum(e_x, axis=-1, keepdims=True)

# def max_CAI_GC_GC3_minT_v1(aa_seq,CAI_strength=2,GC3_strength=1,GC_strength=1,minT_strength=10,CAI_per=2,GC_per=1,GC3_per=3,minT_per=0,C2W=C2W_human):
# 	p_dict={}
# 	for j in codon_dict.keys():
# 		p1=np.array(list(map(lambda x:C2W[x],codon_dict[j])))**CAI_strength
# 		p1=p1/sum(p1) if sum(p1) != 0 else p1
# 		p2=np.array(list(map(lambda x:(x.count('G')+x.count('C')),codon_dict[j])))**GC_strength
# 		p2=p2/sum(p2)
# 		p3=np.array(list(map(lambda x:(x[2].count('G')+x[2].count('C')),codon_dict[j])))**GC3_strength
# 		p3=p3/sum(p3) if sum(p3) != 0 else p3
# 		p4=np.array(list(map(lambda x:(3-x.count('T')),codon_dict[j])))**minT_strength
# 		p4=p4/sum(p4) if sum(p4) != 0 else p4
# 		p=(CAI_per*p1+GC_per*p2+GC3_per*p3+minT_per*p4)
# 		p_dict[j]=p/sum(p)
# 		if sum(p)==0:
# 			raise ValueError("密码子%sp值全为0！" %(j))
# 	res=''.join([np.random.choice(codon_dict[i],p=p_dict[i]) for i in list(aa_seq)])
# 	return(res)

def max_CAI_GC_GC3_minT_v2(aa_seq,CAI_strength=2,GC3_strength=1,GC_strength=1,minT_strength=10,C2W=C2W_human):
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

if args.input_fa:
	fa=read_fasta(input_fa)
elif not args.report:
	fa = {'tmp': sys.stdin.read().strip().upper().replace(' ','')}

if args.report:
	if args.report.endswith(".csv"):
		dr = pd.read_csv(args.report, header=None)
	elif args.report.endswith(".tsv"):
		dr = pd.read_csv(args.report, sep = "\t", header=None)
	elif args.report.endswith(".xlsx"):
		dr = pd.read_excel(args.report, header=None)
	else:
		dr = pd.read_csv(args.report, sep = "\t", header=None)

if not args.C2W:
	C2W_used=C2W_human
else:
	if args.C2W.endswith(".csv"):
		df = pd.read_csv(args.C2W, header=None)
		if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
			df = pd.read_csv(args.C2W)
			if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U')!= 3:
				raise ValueError("tsv文件内容错误！")
	elif args.C2W.endswith(".tsv"):
		df = pd.read_csv(args.C2W, sep = "\t", header=None)
		if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
			df = pd.read_csv(args.C2W, sep="\t")
			if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U')!= 3:
				raise ValueError("tsv文件内容错误！")
	else:
		df = pd.read_excel(args.C2W, header=None)
		if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
			df = pd.read_excel(args.C2W)
			if len(df.iloc[0, 0]) != 3 or df.iloc[0, 0].count('A') + df.iloc[0, 0].count('T') + df.iloc[0, 0].count('C') + df.iloc[0, 0].count('G') + df.iloc[0, 0].count('U') != 3:
				raise ValueError("excel文件内容错误！")
	C2W_used={}
	for i in range(len(df)):
		codon = df.iloc[i,0].upper().replace('U','T')
		if len(codon) != 3 or codon.count('A')+codon.count('T')+codon.count('C')+codon.count('G') != 3:
			raise ValueError(codon + '不是合法的codon！')
		if codon in C2W_used:
			raise ValueError(codon+'重复出现！')
		usage = df.iloc[i,1]
		C2W_used[codon] = usage

if not args.report:
	restriction_pattern=r'A{8,}|T{8,}|C{8,}|G{8,}|U{8,}|TAATACGACTCACTATAAG|CTTATAGTGAGTCGTATTA|'+'|'.join(args.restriction_sites) if len(args.restriction_sites)>0 else r'A{8,}|T{8,}|C{8,}|G{8,}|U{8,}|TAATACGACTCACTATAAG|CTTATAGTGAGTCGTATTA'
	restriction_pattern=restriction_pattern.upper()
	#print('restriction_pattern: '+restriction_pattern)
	ID=[]
	AA=[]
	opts=[]
	CAI=[]
	GC=[]
	GC3=[]
	T=[]
	AUP14s=[]
	if not args.roll:
		if args.utr5:
			utr3 = args.utr3 if args.utr3 else ''
			for cds in fa:
				# 如果是DNA，则先翻译
				if not re.search(r'[^ATCGU]',fa[cds].strip().upper()):
					aa=translate(fa[cds].strip().upper())
				else:
					aa=fa[cds].strip().upper()
				AA.append(aa)
				
				dups=[]
				dup_AUP14s=[]

				def optimize_once(thread_id):
					res = max_CAI_GC_GC3_minT_v2(aa,CAI_strength=args.CAI_s+1,GC3_strength=args.GC3_s+1,GC_strength=args.GC_s+1,minT_strength=args.minT_s+1,C2W=C2W_used)
					# roll掉限制位点
					match=re.search(restriction_pattern, res)
					while match:
						print(res[match.start():match.end()])
						alter_codons=list(range(int(match.start()/3),int((match.end()-1)/3)+1))
						tmp_aa=''.join([dict2[res[(3*i):(3*i+3)]] for i in alter_codons])
						print(tmp_aa)
						while True:
							new=max_CAI_GC_GC3_minT_v2(tmp_aa,CAI_strength=args.CAI_s+1,GC3_strength=args.GC3_s+1,GC_strength=args.GC_s+1,minT_strength=args.minT_s+1,C2W=C2W_used)
							if not re.search(restriction_pattern, new):
								res=res[:(int(match.start()/3)*3)]+new+res[3*(int((match.end()-1)/3)+1):]
								break
						match=re.search(restriction_pattern, res)
					return(res, AUP14(args.utr5+res+utr3))
				
				# 每条重复生成times次，取AUP14最大
				with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
					results = list(executor.map(optimize_once, range(args.times)))
					print(results)# debug
					dups=[result[0] for result in results]
					dup_AUP14s=[result[1] for result in results]
				best=dups[dup_AUP14s.index(max(dup_AUP14s))]
				opts.append(best)
				CAI.append(cai(best,C2W=C2W_used))
				GC.append((best.count('G')+best.count('C'))/len(best))
				tmp=best[2::3]
				GC3.append((tmp.count('G')+tmp.count('C'))/len(tmp))
				T.append(best.count('T')/len(best))
				AUP14s.append(max(dup_AUP14s))

			result=pd.DataFrame({'id': list(fa.keys()), 'aa': AA, 'optimized': opts, 'CAI': CAI, 'GC': GC, 'GC3': GC3, 'T': T, 'AUP14': AUP14s})
			result.to_csv(output,index=False)
		else:
			for cds in fa:
				# 如果是DNA，则先翻译
				if not re.search(r'[^ATCGU]',fa[cds].strip().upper()):
					aa=translate(fa[cds].strip().upper())
				else:
					aa=fa[cds].strip().upper()
				
				def optimize_once(thread_id):
					res = max_CAI_GC_GC3_minT_v2(aa,CAI_strength=args.CAI_s+1,GC3_strength=args.GC3_s+1,GC_strength=args.GC_s+1,minT_strength=args.minT_s+1,C2W=C2W_used)
					# roll掉限制位点
					match=re.search(restriction_pattern, res)
					while match:
						print(res[match.start():match.end()])
						alter_codons=list(range(int(match.start()/3),int((match.end()-1)/3)+1))
						tmp_aa=''.join([dict2[res[(3*i):(3*i+3)]] for i in alter_codons])
						print(tmp_aa)
						while True:
							new=max_CAI_GC_GC3_minT_v2(tmp_aa,CAI_strength=args.CAI_s+1,GC3_strength=args.GC3_s+1,GC_strength=args.GC_s+1,minT_strength=args.minT_s+1,C2W=C2W_used)
							if not re.search(restriction_pattern, new):
								res=res[:(int(match.start()/3)*3)]+new+res[3*(int((match.end()-1)/3)+1):]
								break
						match=re.search(restriction_pattern, res)
					return(res)
				
				with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
					results = list(executor.map(optimize_once, range(args.rep)))
				
				AA+=[aa]*args.rep
				opts+=results
				CAI+=[cai(res,C2W=C2W_used) for res in results]
				GC+=[(res.count('G')+res.count('C'))/len(res) for res in results]
				GC3+=[(res[2::3].count('G')+res[2::3].count('C'))/len(res[2::3]) for res in results]
				T+=[res.count('T')/len(res) for res in results]
				ID+=[cds]*args.rep
			result=pd.DataFrame({'id': ID, 'aa': AA, 'optimized': opts, 'CAI': CAI, 'GC': GC, 'GC3': GC3, 'T': T})
			result.to_csv(output,index=False)
	else:# roll
		for cds in fa:
			# 如果不是DNA，报错
			if re.search(r'[^ATCGU]',fa[cds].strip().upper()):
				raise ValueError('当指定--roll选项时，输入必须为DNA序列！')
			else:
				aa=translate(fa[cds].strip().upper())
			AA.append(aa)
			res = fa[cds].strip().upper().replace('U','T')
			# roll掉限制位点
			match=re.search(restriction_pattern, res)
			while match:
				#print(res[match.start():match.end()])
				alter_codons=list(range(int(match.start()/3),int((match.end()-1)/3)+1))
				tmp_aa=''.join([dict2[res[(3*i):(3*i+3)]] for i in alter_codons])
				while True:
					new=max_CAI_GC_GC3_minT_v2(tmp_aa,CAI_strength=args.CAI_s+1,GC3_strength=args.GC3_s+1,GC_strength=args.GC_s+1,minT_strength=args.minT_s+1,C2W=C2W_used)
					if not re.search(restriction_pattern, new):
						res=res[:(int(match.start()/3)*3)]+new+res[3*(int((match.end()-1)/3)+1):]
						break
				match=re.search(restriction_pattern, res)
			opts.append(res)
			CAI.append(cai(res,C2W=C2W_used))
			GC.append((res.count('G')+res.count('C'))/len(res))
			tmp=res[2::3]
			GC3.append((tmp.count('G')+tmp.count('C'))/len(tmp))
			T.append(res.count('T')/len(res))
		print("sequence after removing restriction sites:\n"+res)
		#print(res)
		mutated=sum(1 for char1, char2 in zip(res, fa[cds].strip().upper().replace('U','T')) if char1 != char2)
		print(f"number of mutated nucleotides: {mutated} ({100*mutated/len(fa[cds]):.3f}%)")
		result=pd.DataFrame({'id': list(fa.keys()), 'aa': AA, 'optimized': opts, 'CAI': CAI, 'GC': GC, 'GC3': GC3, 'T': T})
		if args.input_fa or args.output:
			result.to_csv(output,index=False)
else:
	use_linear=True
	def report_one(i):
		id=str(dr.iloc[i][0]).strip()
		utr5=dr.iloc[i][1].strip().upper().replace('U','T') if isinstance(dr.iloc[i][1],str) else ''
		cds=dr.iloc[i][2].strip().upper().replace('U','T') if isinstance(dr.iloc[i][2],str) else ''
		utr3=dr.iloc[i][3].strip().upper().replace('U','T') if isinstance(dr.iloc[i][3],str) else ''
		seq=utr5+cds+utr3
		tmp=dr.iloc[i][4].strip().upper().replace('U','T') if (dr.shape[1]>=5 and isinstance(dr.iloc[i][4],str)) else ''
		restriction_pattern=r'A{8,}|T{8,}|C{8,}|G{8,}|U{8,}|TAATACGACTCACTATAAG|CTTATAGTGAGTCGTATTA|'+'|'.join([j.strip() for j in re.split('[ ,]+',tmp)]) if len(tmp)>0 else r'A{8,}|T{8,}|C{8,}|G{8,}|U{8,}|TAATACGACTCACTATAAG|CTTATAGTGAGTCGTATTA'
		restriction_pattern=restriction_pattern.upper()
		# 如果不是DNA，报错
		if re.search(r'[^ATCGU]',cds):
			print(cds)
			raise ValueError('第3列需为DNA序列！')
		else:
			aa=translate(cds)
		# 检测限制位点
		it = list(re.finditer(restriction_pattern,cds))
		if len(it)==0:
			# restriction_sites.append('')
			restriction_site=''
		else:
			# restriction_sites.append(', '.join([r"({match.span()[0]},{match.group()})" for match in it]))
			restriction_site=', '.join([f"({match.span()[0]},{match.group()})" for match in it])

		tmp=cds[2::3]
		UP_vec=1-np.sum(bpps(seq, package='vienna',linear=use_linear),axis=0)
		structure,dG_MFE=mfe(seq,package='vienna',linear=use_linear,return_dG_MFE=True)

		return([id, utr5, cds, utr3, aa, restriction_site, cai(cds,C2W=C2W_used), (cds.count('G')+cds.count('C'))/len(cds), (tmp.count('G')+tmp.count('C'))/len(tmp), cds.count('T')/len(cds), free_energy(seq,package='vienna',linear=True), dG_MFE, sum(UP_vec[0:14])/14, sum(UP_vec[max(len(utr5)-16,0):(len(utr5)+21)])/(len(utr5)+21-max(len(utr5)-16,0))])

	with concurrent.futures.ThreadPoolExecutor(max_workers=150) as executor:
		results = list(executor.map(report_one, range(len(dr))))
		result = list(zip(*results))
		result = pd.DataFrame({'id': result[0], 'UTR5': result[1], 'CDS': result[2], 'UTR3': result[3], 'aa': result[4], 'restriction_sites': result[5], 'CAI': result[6], 'cds_GC': result[7], 'cds_GC3': result[8], 'cds_T': result[9], 'free_energy': result[10], 'dG_MFE':result[11], 'AUP14':result[12], 'AUP_ATG':result[13]})
		result.to_csv(output,index=False)
