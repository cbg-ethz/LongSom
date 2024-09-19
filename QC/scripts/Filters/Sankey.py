import pandas as pd
import numpy as np
import plotly.graph_objects as go

def MultiAllelic_filtering(ALT, FILTER, CTYPES, BC, CC, VAF, CCF):

	if len(ALT.split('|'))>1:
		if len(CTYPES.split(','))>1:
			if len(BC.split(',')[1].split('|'))>1:
				BCS = [int(i) for i in BC.split(',')[1].split('|')]
				MAX = max(BCS)
				index = np.argmax(BCS)
				BCS[index] = 0 # removing max to select next "max"
				MAX2 = max(BCS)
				if not(MAX2/MAX<0.02): # one alt 50x larger than the other)
					return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
			else:
				return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
				
			ALT_Cancer = ALT.split(',')[1].split('|')[index]
			BC_Cancer = BC.split(',')[1].split('|')[index]
			CC_Cancer = CC.split(',')[1].split('|')[index]
			VAF_Cancer = VAF.split(',')[1].split('|')[index]
			CCF_Cancer = CCF.split(',')[1].split('|')[index]
			try:
				ALT_NonCancer = ALT.split(',')[0].split('|')[index]
				BC_NonCancer = BC.split(',')[0].split('|')[index]
				CC_NonCancer = CC.split(',')[0].split('|')[index]
				VAF_NonCancer = VAF.split(',')[0].split('|')[index]
				CCF_NonCancer = CCF.split(',')[0].split('|')[index]
			except:
				ALT_NonCancer = ALT.split(',')[0].split('|')[0]
				BC_NonCancer = BC.split(',')[0].split('|')[0]
				CC_NonCancer = CC.split(',')[0].split('|')[0]
				VAF_NonCancer = VAF.split(',')[0].split('|')[0]
				CCF_NonCancer = CCF.split(',')[0].split('|')[0]
			ALT = ','.join([ALT_NonCancer,ALT_Cancer])
			BC = ','.join([BC_NonCancer,BC_Cancer])
			CC = ','.join([CC_NonCancer,CC_Cancer])
			VAF = ','.join([VAF_NonCancer,VAF_Cancer])
			CCF = ','.join([CCF_NonCancer,CCF_Cancer])
			return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'KEEP'
		else:
			BCS = [int(i) for i in BC.split('|')]
			if min(BCS)/max(BCS)<0.02: # one alt 50x larger than the other)
				index = np.argmax(BCS)
			else:
				return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
				
			ALT = ALT.split('|')[index]
			BC = BC.split('|')[index]
			CC = CC.split('|')[index]
			VAF = VAF.split('|')[index]
			CCF = CCF.split('|')[index]
			return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'KEEP'
	else:	
		if 'Multi-allelic' in FILTER: #Dissonant alleles
			return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'DELETE'
		return ALT, FILTER, CTYPES, BC, CC, VAF, CCF,'KEEP'

def colnames(file):
	with open(file, 'r') as f:
		for line in f:
			if line.startswith('#'):
				# Extract output file column names
				if '#CHROM' in line:
					line = line.rstrip('\n')
					output_column_names = line.split('\t')
			else:
				break
	return output_column_names


def idk(input_df):

	input_df['INDEX'] = input_df['#CHROM'].astype(str) + ':' + input_df['Start'].astype(str) + ':' + input_df['ALT'].str.split(',', n=1, expand=True)[0]
	
	# Only keeping sites with alternative reads in cancer
	input_df = input_df[input_df['Cell_types']!='NonCancer']

	# Filtering multiallelic sites, keeping only one allele if it's 50x larger than the other
	# Otherwise deleting the line (it messes with filters later on)
	# This is important as SComatic will often detect low expr secondary alt allele in very high expr. (e.g. in chrM)
	input_df[['ALT', 'FILTER', 'Cell_types', 'Bc', 'Cc', 'VAF', 'CCF','MultiAllelic_filter']] =  input_df.apply(lambda x: MultiAllelic_filtering(x['ALT'], x['FILTER'], x['Cell_types'],x['Bc'], x['Cc'], x['VAF'], x['CCF']), axis=1, result_type="expand") 
	input_df = input_df[input_df['MultiAllelic_filter']=='KEEP']
	input_df = input_df[output_column_names + ['INDEX','ID']]

	TOT = len(input_df)
	print('TOT',TOT)
	
	# Filter 1: Non-cancer coverage filter
	input_df = input_df[~input_df['FILTER'].str.contains('Min_cell_types')]
	PASS_F1 = len(input_df)
	MIN_CTYPE = TOT - PASS_F1
	
	# Filter 2: Alt reads and cells in cancer filter
	input_df['Bc'] = [i.split(',')[1] if ',' in i else i for i in input_df['Bc']] #only keeping cancer info
	input_df['Cc'] = [i.split(',')[1] if ',' in i else i for i in input_df['Cc']] #only keeping cancer info
	input_df = input_df.astype({'Bc':'int','Cc':'int'})
	input_df = input_df[input_df['Bc']>=3]
	input_df = input_df[input_df['Cc']>=2]
	PASS_F2 = len(input_df)
	MINMUT = PASS_F1 - PASS_F2

	# Filter 3: Betabinomial filter in cancer cells
	input_df = input_df[~input_df['Cell_type_Filter'].str.contains(',Non-Significant|,Low-Significant', regex = True)]
	input_df = input_df[~input_df['Cell_type_Filter'].isin(['Non-Significant','Non-Significant'])]
	PASS_F3 = len(input_df)
	BETABIN_CANCER = PASS_F2 - PASS_F3

	# Filter 4: Noise filter:
	input_df = input_df[~input_df['FILTER'].str.contains('Noisy_site')] #multi allelic sites are filtered above
	PASS_F4 = len(input_df)
	NOISE =  PASS_F3 - PASS_F4

	# Filter 5: Homopolymer filter:
	input_df = input_df[~input_df['FILTER'].str.contains('LC_Upstream|LC_Downstream', regex = True)] #multi allelic sites are filtered above
	PASS_F5 = len(input_df)
	HOMO =  PASS_F4 - PASS_F5

	# Filter 6:RNA-Editing filter
	input_df = input_df[~input_df['FILTER'].str.contains('RNA_editing_db', regex = True)]
	PASS_F6 = len(input_df)
	EDIT =  PASS_F5 - PASS_F6

	# Filter 7: PoN filter
	input_df = input_df[~input_df['FILTER'].str.contains('PoN', regex = True)]
	PASS_F7 = len(input_df)
	PON = PASS_F6 - PASS_F7

	# Filter 8: Betabinomial filter in non-cancer cells
	input_df = input_df[~input_df['Cell_type_Filter'].str.contains('Low-Significant,|PASS,', regex = True)]
	input_df = input_df[~input_df['FILTER'].str.contains('Cell_type_noise|Multiple_cell_types|Noisy_site', regex = True)]
	PASS_F8 = len(input_df)
	BETABIN_NONCANCER =  PASS_F7 - PASS_F8

	# Filter 9: gnomAD filter
	input_df = input_df[~input_df['FILTER'].str.contains('gnomAD', regex = True)]
	PASS_F9 = len(input_df)
	GNOMAD =  PASS_F8 - PASS_F9
	
	input_df['FINAL_FILTER'] = tag_clustered_SNVs(input_df, 10000)
	DIST = len(input_df[input_df['FINAL_FILTER']=='Clust_dist10000'])
	PASS_F10 = len(input_df[input_df['FINAL_FILTER']=='PASS'])

	return [TOT,MIN_CTYPE,PASS_F1,MINMUT,PASS_F2,BETABIN_CANCER,PASS_F3,NOISE,PASS_F4,HOMO,PASS_F5,EDIT,PASS_F6,PON,PASS_F7,BETABIN_NONCANCER,PASS_F8,GNOMAD,PASS_F9,DIST,PASS_F10]

def tag_clustered_SNVs(input_df_tot, clust_dist):
	l = []
	for id in ['B486','B497','B500']:
		input_df = input_df_tot[input_df_tot['ID']==id]
		input_df2 = input_df[input_df['FILTER']=='PASS']
		idx = input_df2['INDEX']
		a=[]
		for i in idx:
			chr,pos,base=i.split(':')
			a.append((chr,pos,base))
		b = sorted(a, key=lambda x: (x[0],x[1]))

		trash=[]

		for (chr,pos,base), (chr2,pos2,base2) in zip(b, b[1:]):
			if chr==chr2:
				if chr == 'chrM':
					continue
				if abs(int(pos)-int(pos2))<clust_dist:
					trash.append(':'.join([chr,pos,base]))
					trash.append(':'.join([chr2,pos2,base2]))
		trash = set(trash)
		input_df['FINAL_FILTER'] = input_df.apply(lambda x : 
								modify_filter(x['INDEX'],x['FILTER'], clust_dist, trash),
								axis=1)
		l+=list(input_df['FINAL_FILTER'])
		
	return l

def modify_filter(INDEX, FILTER, clust_dist, trash):
	clustered = 'Clust_dist{}'.format(str(clust_dist))
	if INDEX in trash:
		if FILTER == 'PASS':
			return clustered
		else:
			return FILTER + ',' + clustered
	else:
		return FILTER

p = '/cluster/work/bewi/members/dondia/projects/long_reads_tree/LongSom_out/OvCa_LR/SNVCalling/BaseCellCalling/'

output_column_names = colnames(p+'B486.calling.step2.tsv')
input_dfs = []
for id in ['B486','B497','B500']:
	input_df =  pd.read_csv(p+id+'.calling.step2.tsv', sep='\t',comment='#',names=output_column_names)
	input_df['ID'] = id
	print(id,len(input_df))
	input_dfs.append(input_df)

input_df = pd.concat(input_dfs)
print('all input_df',len(input_df))
values = idk(input_df)

labels =  ["Candidates", "Low non-cancer<br>coverage","Pass<br>filter 1", "Low cancer<br>variant reads", "Pass<br>filter 2", "Failed cancer<br>beta-binom test",  "Pass<br>filter 3", 'Multi-allelic<br>noise',  "Pass<br>filter 4", "Homopolymer", "Pass<br>filter 5","In RNA-edit db", "Pass<br>filter 6","In PoN SR/LR","Pass<br>filter 7", "Failed non-cancer<br>beta-binom test", "Pass<br>filter 8", "In gnomAD db", "Pass<br>filter 9","Low distance<br>between loci", "Pass filter 10<br>Somatic calls"]
#               0                         1                 2                        3                            4                           5                         6                         7                      8                   9              10                   11                12             13               14                       15                              16                17              18                      19                            20

label_to_value = dict(zip(labels,values))
labels = [f"{k}<br>({v})" for k,v in label_to_value.items()]

fig = go.Figure(data=[go.Sankey(arrangement = "snap",
    node = dict(
      pad = 30,
      thickness = 10,
      line = dict(color = "black", width = 0.5),
      x=[0, 0.09, 0.09, 0.18, 0.18, 0.27, 0.27, 0.36, 0.36, 0.45, 0.45, 0.54, 0.54, 0.63, 0.63, 0.72, 0.72, 0.81, 0.81, 0.9, 1], 
      y=[0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.55],
      label = labels,
	  color=['green'] + ['red','green']*10
      ),
    link = dict(
      source = [0, 0, 2, 2, 4, 4, 6, 6, 8, 8,  10, 10, 12, 12, 14, 14, 16, 16, 18, 18], 
      target = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
      value = values[1:]   
  ))])

fig.update_layout(title_text="LongSom SNV Loci Filtering Workflow", font_size=18)
fig.write_image("/cluster/work/bewi/members/dondia/projects/long_reads_tree/LongSom_out/OvCa_LR/QC/Sankey/SankeyFilteringWorkflow.png",width=2000, height=400)
