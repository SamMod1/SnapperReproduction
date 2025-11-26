import pandas as pd
import csv
from ete3 import Tree


def progress_bar(current, total):
    percent = int(round((current / total) * 100))
    return '#' * int(percent / 10) + '-' * (10-int(percent/10)) + ' (' + str(percent) + ' %)'


def restore_sequence_names(target_file, mappings_text, reverse = False):
    mappings_text = mappings_text.split('\n')
    mappings_text = [x for x in mappings_text if x]
    mappings = {}
    
    for mapping in mappings_text:
        mapping = mapping.split('\t')
        if reverse:
            mappings[mapping[1]] = mapping[0]
        else:
            mappings[mapping[0]] = mapping[1]
        
    for key in mappings:
        target_file = target_file.replace(mappings[key] + '\t', key + '\t')
    
    return target_file


def find_insert_position(arr, target):
    """
    Perform a binary search to find the index where target falls within a sorted array arr.

    :param arr: Sorted list of integers.
    :param target: Target integer to find the insert position for.
    :return: Index at which target should be inserted to maintain sorted order.
    """
    left, right = 0, len(arr)

    while left < right:
        mid = (left + right) // 2
        if arr[mid] < target:
            left = mid + 1
        else:
            right = mid

    return left



def feature_overlap(target_start, target_end, starts_sorted, ends_sorted):
    """
    Given the start and end of a feature, find all overlapping features from two sorted lists of other features
    
    :param target_start: Integer - The start position of the target feature sequence
    :param target_end: Integer - The end position of the target feature sequence
    :param starts_sorted: Nested List - Each element contains an iterable with the name of features as the first
    element in each element and the start location as the second element upon which the list is sorted (ascending)
    :param ends_sorted: Nested List - Each element contains an iterable with the name of features as the first
    element in each element and the end location as the second element upon which the list is sorted (ascending)
    :return: List containing the names of all overlapping features
    """
    start_before_end = find_insert_position([x[1] for x in starts_sorted], target_end)
    start_before_end = set([x[0] for x in starts_sorted[:start_before_end]])
    end_after_start = find_insert_position([x[1] for x in ends_sorted], target_start)
    end_after_start = set([x[0] for x in ends_sorted[end_after_start:]])
    
    return start_before_end.intersection(end_after_start)


def cds_overlap_clusters(gff, features_to_cluster):
    cds_df = pd.DataFrame()
    clusters = []
    
    for feature in features_to_cluster:
        feature_df = gff[feature].to_df()[['feature_type', 'start', 'end', 'ID']]
        feature_df['start'] = feature_df['start'].astype(int)
        feature_df['end'] = feature_df['end'].astype(int)
        feature_df = feature_df.loc[feature_df['feature_type'] == 'CDS']
        if not feature_df.empty:
            feature_df['_count_'] = [str(x) for x in range(len(feature_df))]
            feature_df['transcript'] = feature
            feature_df['cds_id'] = feature_df['_count_'] + '_' + feature
            cds_df = pd.concat([cds_df, feature_df], ignore_index=True)
                    
    for feature in features_to_cluster:
        feature_df = cds_df.loc[cds_df['transcript'] == feature][['start', 'end']]
        
        if feature_df.empty:
            continue
        
        found = False
        for idx, x in enumerate(clusters):
            if feature in x:
                transcript_group = clusters[idx]
                found = True
                break
            
        if found:
            other_transcripts = cds_df.loc[~(cds_df['transcript'].isin(transcript_group))][['start', 'end', 'cds_id']]
        else:
            other_transcripts = cds_df.loc[~(cds_df['transcript'] == feature)][['start', 'end', 'cds_id']]
            
        starts_sorted = other_transcripts.sort_values(by='start')[['cds_id', 'start']].to_numpy()
        ends_sorted = other_transcripts.sort_values(by='end')[['cds_id', 'end']].to_numpy()
        
        overlaps = set()
        for idx in range(len(feature_df)):
            cds = feature_df.iloc[idx]
            start = cds['start']
            end = cds['end']
            cds_overlaps = feature_overlap(start, end, starts_sorted, ends_sorted)
            if cds_overlaps:
                overlaps = overlaps.union(set([x.split('_', 1)[1] for x in cds_overlaps]))
            
        if found:
            transcript_group.extend(list(overlaps))
        else:
            overlaps.add(feature)
            clusters.append(list(overlaps))
            
    return clusters


def parse_orthogroups(orthogroups_file, dr_genes_to_keep = None):
    orthogroup_data = {}
    with open(orthogroups_file, "r") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ca_genes = row["Chrysophrys_auratus"].split(", ")
            dr_genes = row["Danio_rerio"].split(", ")
            if dr_genes_to_keep:
                dr_genes = [x for x in dr_genes if x in dr_genes_to_keep]
                if not dr_genes:
                    dr_genes = ['']
            orthogroup_data[row["Orthogroup"]] = {
                "Chrysophrys_auratus": ca_genes,
                "Danio_rerio": dr_genes
            }
    return orthogroup_data


def find_closest_zebrafish_gene(tree_file, ca_genes, dr_genes):
    tree = Tree(tree_file)
    closest_genes = {}
    
    for ca_gene in ca_genes:
        node = tree.search_nodes(name='Chrysophrys_auratus_' + ca_gene)[0]
        distances = {}
        for dr_gene in dr_genes:
            dr_node = tree.search_nodes(name='Danio_rerio_' + dr_gene)[0]
            distance = node.get_distance(dr_node)
            distances.setdefault(distance, []).append(dr_gene)
            
        smallest_distance = min([x for x in distances.keys()])
        closest_gene = distances[smallest_distance][0]
        closest_genes.setdefault(closest_gene, []).append((ca_gene, smallest_distance))
        
    return closest_genes


def assign_closest_genes(orthogroups_data, gene_trees_dir):
    og_gene_mappings = {}
    
    for orthogroup, genes in orthogroups_data.items():
        og = {}
        ca_genes = genes["Chrysophrys_auratus"]
        dr_genes = genes["Danio_rerio"]
        if not ca_genes[0] or not dr_genes[0]:
            continue
        
        
        if len(ca_genes) == 1 and len(dr_genes) == 1:
            og[dr_genes[0]] = [(ca_genes[0],)]
                    
        else:
            tree_file = f"{gene_trees_dir}/{orthogroup}_tree.txt"
            print('\r' + orthogroup, end = '\r')
            og.update(find_closest_zebrafish_gene(tree_file, ca_genes, dr_genes))
        
        
        # if len(dr_genes) == 1:
        #     if len(ca_genes) == 1:
        #         og[dr_genes[0]] = [ca_genes[0]]
        #     else:
        #         og[dr_genes[0]] = []
        #         for ca_gene in ca_genes:
        #             og[dr_genes[0]].append(ca_gene)
                    
        # else:
        #     tree_file = f"{gene_trees_dir}/{orthogroup}_tree.txt"
        #     print('\r' + orthogroup, end = '\r')
        #     og.update(find_closest_zebrafish_gene(tree_file, ca_genes, dr_genes))
                
        og_gene_mappings[orthogroup] = og
    
    return og_gene_mappings


def rename_genes(og_mappings):
    new_mappings = {}
    
    for orthogroup, gene_mappings in og_mappings.items():
        for dr_gene, ca_genes in gene_mappings.items():
            if len(ca_genes) == 1:
                new_mappings[ca_genes[0][0]] = dr_gene
            elif len(ca_genes) > 1:
                ca_genes.sort(key = lambda x: x[1])
                for count, ca_gene in enumerate(ca_genes):
                    new_mappings[ca_gene[0] + '%%' + str(count)] = dr_gene
    
    return new_mappings
