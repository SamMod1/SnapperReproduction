import create_annotation
import sys

INPUT_ANNOTATION = sys.argv[1]
GENE_INFO = sys.argv[2]
ORTHOFINDER_RESULTS = sys.argv[3]
OUTPUT_DIR = sys.argv[4]
OUTPUT_FILENAME = sys.argv[5]
SPECIES_NAME = sys.argv[6]

GENE_PREFIXES_TO_IGNORE = set(['LOC', 'si:', 'zgc:', 'zmp:', '7955_'])

ORTHOGROUPS = ORTHOFINDER_RESULTS + 'Orthogroups/Orthogroups.tsv'
GENE_TREES_DIR = ORTHOFINDER_RESULTS + 'Gene_Trees'


gene_info = open(GENE_INFO, 'r').readlines()

dr_gene_data = {}
count = 0

for line in gene_info[1:]:
    line = line.split('\t')
    odb_gene = line[5].replace(':', '_')
    gene_symbol = line[2]
    gene_name = line[3]
    
    if not any(gene_symbol.startswith(prefix) for prefix in GENE_PREFIXES_TO_IGNORE):
        if odb_gene in dr_gene_data.keys():
            count += 1
        else:
            dr_gene_data[odb_gene] = (gene_symbol, gene_name)

orthogroups_data = create_annotation.parse_orthogroups(ORTHOGROUPS, dr_genes_to_keep=dr_gene_data)
#orthogroups_data = {og: orthogroups_data[og] for og in [x for x in orthogroups_data.keys()][2500:2520]}  # Uncomment for testing
og_gene_mappings = create_annotation.assign_closest_genes(orthogroups_data, GENE_TREES_DIR)
gene_mappings = create_annotation.rename_genes(og_gene_mappings)
new_gene_mappings = {}
gene_full_names = {}

for ca_gene, dr_gene in gene_mappings.items():
    real_dr_gene = dr_gene_data[dr_gene][0]
    if '%%' in ca_gene:
        ca_gene, gene_number = ca_gene.split('%%')
        real_dr_gene += '_' + SPECIES_NAME + str(int(gene_number) + 1)
    gene_full_names[ca_gene] = dr_gene_data[dr_gene][1]
    new_gene_mappings[ca_gene] = real_dr_gene
    
gene_mappings = new_gene_mappings
    
orthogroups = open(ORTHOGROUPS, 'r').read()
for mapping in dr_gene_data:
    orthogroups = orthogroups.replace(mapping, dr_gene_data[mapping][0])
    
open(OUTPUT_DIR + 'orthogroups_with_names.tsv', 'w').write(orthogroups)

annot = create_annotation.Gff(filename = INPUT_ANNOTATION)

mappings_len = len(gene_mappings)
for count, mapping in enumerate(gene_mappings):
    print('\r' + create_annotation.progress_bar(count, mappings_len), end='\r')
    gene = annot.get_feature(annot.get_feature(mapping)['Parent'])
    origional_gene = gene['ID'].replace('gene-', '')
    annot.add_attribute(gene['ID'], 'full_name', gene_full_names[mapping].replace(';', ','))
    if not gene_mappings[mapping] == origional_gene:
        try:
            annot.find_replace(gene['ID'], gene_mappings[mapping], origional_gene)
        except KeyError:
            gene_name = gene_mappings[mapping].split('_' + SPECIES_NAME)
            gene_num = 1
            if len(gene_name) > 1:
                gene_num = int(gene_name[1])
            annot.find_replace('gene-' + gene_mappings[mapping], 
                               gene_mappings[mapping] + '_'  + SPECIES_NAME + str(gene_num), gene_mappings[mapping])
            annot.find_replace(gene['ID'], gene_mappings[mapping] + '_'  + SPECIES_NAME + str(gene_num + 1), origional_gene)
                
    
open(OUTPUT_DIR + OUTPUT_FILENAME, 'w').write(annot.to_gff(True))
