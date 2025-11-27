import pandas as pd
from numpy import mean, log2
import seaborn as sns


EXPRESSION_FILE = 'DataChapterTwo/Data/salmon.merged.gene_counts_gonad_length_scaled.tsv'
# GROUPS = ['SEX.FEMALE_STAGE.1_PRE', 'SEX.FEMALE_STAGE.2_VIT', 'SEX.FEMALE_STAGE.3_MAT', 'SEX.FEMALE_STAGE.4_REG']
# READABLE_GROUPS = ['Previtellogenic', 'Vitellogenic', 'Mature', 'Regression']
# GROUPS = ['SEX.UNKNOWN', 'SEX.MALE', 'SEX.FEMALE']
GROUPS = ['_F', '_M', '1_J', '2_J', '3_J', '4_J', '5_J', '6_J', '7_J', '8_J', '9_M', '10_J', '11_J', '12_J', '14_J', '15_J', '16_J', '17_J', '18_J', '20_J', '21_J', '35_F', '36_F']
READABLE_GROUPS =['Female', 'Male', '1_UN', '2_UN', '3_UN', '4_UN', '5_M', '6_UN', '7_M', '8_UN', '9_M', '10_UN', '11_F', '12_UN', '14_UN', '15_M', '16_UN', '17_M', '18_M', '20_F', '21_UN', '35_F', '36_F']
REFERENCE_GROUP = None#'SEX.FEMALE_STAGE.1_PRE'
OUTPUT = 'PCA.html'
TRANSFORMATION = lambda x: log2(x+1)
#TRANSCRIPTS = "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/Python/DifferentialExpression/PCA/all_juv_degs.csv"
TRANSCRIPTS = "DataChapterTwo/DESeq/Results/adult_degs.csv"#['gene-foxl2a', 'gene-amh', 'gene-dmrt1']
FIGSIZE = (10, 12)
FONTSIZE = 1.3
ANNOT = False
SAMPLES_TO_EXCLUDE = ['90_M3', '13_J0']


sns.set(font_scale=FONTSIZE)

df = pd.read_csv(EXPRESSION_FILE, sep='\t')
if type(TRANSCRIPTS) == str:
    transcripts = open(TRANSCRIPTS, 'r').read().replace('\n', '').split(',')
else:
    transcripts = TRANSCRIPTS
df = df.loc[df['gene_id'].isin(transcripts)]
df.index = df['gene_id']
df = df.reindex(transcripts)
df = df[pd.Series(df.columns).loc[~df.columns.isin(SAMPLES_TO_EXCLUDE)]]

group_dfs = {}
for group in GROUPS:
    group_dfs[group] = df[[x for x in df.columns if group in x]]

heatmap_data = {}
gene_names = []

for group in GROUPS:
    heatmap_data[group] = []

for idx in range(len(df)):
    try:
        gene = df.iloc[idx]['gene_id'].replace('gene-', '')
    except AttributeError:
        continue
    gene_names.append(gene)
    
    for group in GROUPS:
        group_df = group_dfs[group]
        line = group_df.iloc[idx]
        
        if TRANSFORMATION and not REFERENCE_GROUP:
            log2_mean = TRANSFORMATION(mean([x for x in line]))
        else:
            log2_mean = mean([x for x in line if x])
        heatmap_data[group].append(round(log2_mean, 2))
        
heatmap_data = pd.DataFrame(heatmap_data, index=gene_names)

if REFERENCE_GROUP:
    old_heatmap_data = heatmap_data.copy()
    for group in GROUPS:
        if TRANSFORMATION:
            heatmap_data[group] = TRANSFORMATION(old_heatmap_data[group] / old_heatmap_data[REFERENCE_GROUP])
        else:
            heatmap_data[group] = old_heatmap_data[group] / old_heatmap_data[REFERENCE_GROUP]

if READABLE_GROUPS:
    heatmap_data.columns = READABLE_GROUPS
    
#fig, ax = plt.subplots(figsize=FIGSIZE)
heatmap = sns.clustermap(heatmap_data, annot = ANNOT, figsize = FIGSIZE)
