import pandas as pd

# Run
PON=config['Run']['PoN']
REANNO=config['Run']['CellTypeReannotation']
SCOMATIC=config['Run']['SNVCalling']
CTATFUSION=config['Run']['FusionCalling']
BNPC=config['Run']['CellClustering']
INFERCNV=config['Run']['CNACalling']

INPUT = config['User']['input_dir']

# import sample map and retrieve sample names
samples_table = pd.read_table(config["User"]["sample_map"], header=0)
samples = samples_table.set_index("sample", drop=False)
IDS = samples_table["sample"].tolist()

def get_BetaBinEstimates(input, value):
    df = pd.read_csv(input, sep='\t')
    d = df.squeeze().to_dict()
    return d[value]