from pywaffle import Waffle
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd

prefix = ''
func = []
for id in ['B486','B497','B500']:
    idx = []
    p = pd.read_csv(prefix+id+'.hg38_multianno.txt',sep='\t',na_values=['.'])
    func+=list(p['Func.refGene'])

func = [i.split(';')[0]  if ';' in i else i for i in func]
func = [i.split('ncRNA_')[1]  if 'ncRNA_' in i else i for i in func]
d = {'upstream':'intergenic','downstream':'intergenic'}
func = [d[i] if i in d else i for i in func]

data = Counter(func)
data = {k:data[k] for k in ['exonic','splicing','UTR3','UTR5', 'intronic','intergenic']}
tot = sum(data.values())
data_div10 = {k:round(v/10) for k,v in data.items()}

fig = plt.figure(
    FigureClass=Waffle,
    rows=25,
    values=data_div10,
    title={'label': '', 'loc': 'left'},
    labels=["{} ({:2.1%})".format(k,v/tot) for k, v in data.items()],
    legend={'loc': 'lower left', 'bbox_to_anchor': (1, 0.5), 'ncol': 1, 'framealpha': 1},
    starting_location='NW',
    vertical=True,
    cmap_name="Set3",
    block_arranging_style='snake',
    figsize = (3,3), dpi=600
    
)

plt.tight_layout()
plt.savefig('Functional_Annotation_legend.png',dpi=600)
plt.close()

fig = plt.figure(
    FigureClass=Waffle,
    rows=32,
    values=data,
    title={'label': 'Functional annotation of LongSom calls', 'loc': 'left','size':40},
    labels=["{} ({:2.1%})".format(k,v/tot) for k, v in data.items()],
    legend={'loc': 'lower left', 'bbox_to_anchor': (1, 0.5), 'ncol': 1, 'framealpha': 1},
    starting_location='NW',
    vertical=True,
    block_arranging_style='normal',
    cmap_name="Set3",
    figsize = (15,15), dpi=600
    
)

plt.tight_layout()
plt.savefig('Functional_Annotation.png',dpi=600)