import create_annotation
import SamBio
import sys


# Things you probably want to change
GENOME = sys.argv[1]
MAPPINGS_FILE = sys.argv[2]
LIFTOFF_GFF = sys.argv[3]
BRAKER_GTF = sys.argv[4]
OUTPUT = sys.argv[5]

# Thinks you might want to change
#OUTPUT = "merge_annotation/"  # Make sure this directory exists already
ANNOTATION_TO_PRIORITISE = "BRAKER"


print('Loading BRAKER annotation...')
braker_gtf = open(BRAKER_GTF, 'r').read()
mappings = open(MAPPINGS_FILE, 'r').read()
braker_gtf = create_annotation.restore_sequence_names(braker_gtf, mappings)
open(OUTPUT + 'braker_names_restored.gtf', 'w').write(braker_gtf)

print("\nConverting BRAKER annotation to GFF format:")
braker_gff = create_annotation.gtf_to_gff(braker_gtf, True)
open(OUTPUT + 'braker.gff', 'w').write(braker_gff)



print('\nLoading liftoff annotation...')
liftoff_gff_lines = open(LIFTOFF_GFF, 'r').readlines()
liftoff_gff = create_annotation.Gff(filelines = liftoff_gff_lines)


print('\nSplitting fused genes in the liftoff annotation:')
liftoff_gff = create_annotation.separate_fused_genes(liftoff_gff, True)
open(OUTPUT + 'liftoff_gff_defused.gff', 'w').write(liftoff_gff.to_gff())


liftoff_gff = create_annotation.Gff(OUTPUT + 'liftoff_gff_defused.gff')
braker_gff = open(OUTPUT + 'braker.gff', 'r').read()


print('\n\nCombining annotations...')
braker_gff_lines = braker_gff.split('\n')
braker_gff = create_annotation.Gff(filelines = braker_gff_lines)

braker_gff.strip_comments()
liftoff_gff.strip_comments()

combined_annotation = liftoff_gff.copy()
combined_annotation.extend(braker_gff)


print('\nMerging annotations:\n')
gene_clusters = create_annotation.cluster_overlapping_genes(combined_annotation, show_progress_bar = True)

new_annotation = create_annotation.Gff()
new_gene_clusters = {}
for sequence in gene_clusters:
    for idx, cluster in enumerate(gene_clusters[sequence]):
        if len(cluster) == 1:
            new_annotation.extend(combined_annotation[list(cluster)[0]])
        else:
            try:
                new_gene_clusters[sequence].append(cluster)
            except KeyError:
                new_gene_clusters[sequence] = [cluster,]
    
chosen_genes = create_annotation.merge_annotation(combined_annotation, new_gene_clusters, True)
new_annotation.extend(chosen_genes)

print('\n\nReordering annotation...')
mega_df = new_annotation.to_df()
reordered_annotation = create_annotation.Gff()
for sequence in gene_clusters:
    sequence_df = mega_df.loc[mega_df['sequence'] == sequence]
    sequence_df = sequence_df.loc[sequence_df['ID'].str.startswith('gene-')]
    sequence_df = sequence_df.sort_values(by='start')
    for gene in sequence_df['ID']:
        reordered_annotation.extend(new_annotation[gene])
        
print('\nSaving annotation...')
open(OUTPUT + 'merged.gff', 'w').write(reordered_annotation.to_gff())

annot = create_annotation.Gff(OUTPUT + 'merged.gff')
genome = SamBio.read_fasta(GENOME)
linked = create_annotation.FastaLink(annot, genome)

proteome = linked.create_proteome(alt_transcripts = False)
problematic_sequences = SamBio.MultiSequence()
ok_sequences = SamBio.MultiSequence()
all_sequences = SamBio.MultiSequence()

for protein in proteome:
    rna_name = protein.name.replace(' protein', '')
    sequence = linked[rna_name]
    if '*' in protein.sequence[:-1]:
        problematic_sequences[rna_name] = sequence
    else:
        ok_sequences[protein.name] = protein
    all_sequences[rna_name] = sequence

ok_sequences.write_fasta(OUTPUT + 'merged_ok.pep')
problematic_sequences.write_fasta(OUTPUT + 'merged_bad.fasta')
all_sequences.write_fasta(OUTPUT + 'merged_all.fasta')
proteome.write_fasta(OUTPUT + 'merged_all.pep')

ok_sequence_text = open(OUTPUT + 'merged_ok.pep', 'r').read()
ok_sequence_without_stop = ok_sequence_text.replace('*', '')
open(OUTPUT + 'merged_ok.pep', 'w').write(ok_sequence_without_stop)
