from .gff import Gff
from .general_funcs import progress_bar, feature_overlap


def extract_sequences(gff, show_progress_bar = False):
    sequences = {}  # The annotation will be split by sequence
    
    for idx, feature_name in enumerate(gff):
        feature = gff.get_feature(feature_name)
        
        if feature['feature_type'] == 'gene':
            if feature['gene_biotype'] == 'protein_coding':
                if show_progress_bar:
                    print('\rExtracting unique sequences: ' + progress_bar(idx, len(gff)), end='\r')
                    
                try:
                    sequences[feature['sequence']].extend(gff[feature_name])
                except KeyError:
                    sequences[feature['sequence']] = Gff()
                    sequences[feature['sequence']].extend(gff[feature_name]) 
    return sequences


def cluster_overlapping_genes(gff, show_progress_bar = False):
    sequences = extract_sequences(gff, show_progress_bar)
        
    if show_progress_bar:
        print('')
        
    clusters = {}
    
    for idx, sequence in enumerate(sequences):
        sequence_clusters = []
        if show_progress_bar:
            print('\rClustering overlapping genes: ' + progress_bar(idx, len(sequences)), end='\r')
            
        sequence_df = sequences[sequence].to_df()
        dfs = [sequence_df.loc[sequence_df['strand'] == '+'], sequence_df.loc[sequence_df['strand'] == '-']]
        for strand_df in dfs:
            strand_df = strand_df[strand_df['feature_type'] == 'gene'][['ID', 'start', 'end']]
            
            if len(strand_df) <= 1:
                if len(strand_df) == 1:
                    sequence_clusters.append(set([strand_df.iloc[0]['ID']]))
                continue
            
            strand_df['start'] = strand_df['start'].astype(int)
            strand_df['end'] = strand_df['end'].astype(int)
            
            starts_sorted = strand_df.sort_values(by='start')[['ID', 'start']].to_numpy()
            ends_sorted = strand_df.sort_values(by='end')[['ID', 'end']].to_numpy()
                    
            for gene_idx in range(len(strand_df)):
                gene = strand_df.iloc[gene_idx]
                # gene_name = gene['ID']
                
                # for idx, cluster in enumerate(sequence_clusters):
                #     if gene_name in cluster:
                #         found = True
                #         gene_cluster = idx
                #         break
                
                # if not found:
                #     sequence_clusters.append([gene_name])
                #     gene_cluster = len(sequence_clusters) - 1
                
                start = gene['start']
                end = gene['end']
                
                overlaps = feature_overlap(start, end, starts_sorted, ends_sorted)
                
                found = False
                for overlap in overlaps:
                    for idx, cluster in enumerate(sequence_clusters):
                        if overlap in cluster:
                            found = True
                            sequence_clusters[idx] = cluster.union(overlaps)
                            break
                    if found:
                        break
                    
                if not found:
                    sequence_clusters.append(overlaps)
                    
        if sequence_clusters:
            clusters[sequence] = sequence_clusters
        
    if show_progress_bar:
        print('\rClustering overlapping genes: ' + progress_bar(100, 100), end='\r')
        print('')
    
    return clusters
