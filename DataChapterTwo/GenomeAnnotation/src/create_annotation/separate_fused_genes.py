from .general_funcs import progress_bar, cds_overlap_clusters
from .gff import Gff


def separate_fused_genes(gff, verbose = False):
    if verbose:
        print(progress_bar(0, 100), end='\r')
    
    gff = gff.copy()
    num_genes = gff.feature_counts()['gene']
    for idx, gene_name in enumerate([x for x in gff if x[:5] == 'gene-']):
        if verbose:
            print('\r' + progress_bar(idx, num_genes), end='\r')
        gene = gff.get_feature(gene_name)
        if not gene['gene_biotype'] == 'protein_coding':
            continue
        
        gene_gff = gff[gene_name]
        clusters = cds_overlap_clusters(gene_gff, [x for x in gene_gff if x[:4] == 'rna-'])
        
        if len(clusters) == 1:
            continue
        
        if verbose:
            print('\rFused gene found: ' + gene['ID'])
            print('\r' + progress_bar(idx, num_genes), end='\r')
        
        gene_idx = gff.index(gene_name)
        gene_gff = gff.pop(gene_name)
        
        for idx, cluster in enumerate(clusters):
            new_gene_id = gene['ID'] + '_fusion_' + str(len(clusters) - idx)
            new_gene_name = new_gene_id.replace('gene-', '')
            new_gene = dict(gene)
            new_gene['ID'] = new_gene_id
            new_gene['Name'] = new_gene_name
            new_gene_gff = Gff()
            new_gene_gff.append(new_gene)
            
            for transcript in cluster:
                transcript_gff = gene_gff[transcript].copy()
                transcript_gff.edit_feature_data(transcript, 'Parent', new_gene_id)
                #transcript_gff.remove_single_feature(gene_name)
                new_gene_gff.extend(transcript_gff)
                
            new_gene_gff.edit_attributes('gene', new_gene_name)
            
            gff.insert(new_gene_gff, gene_idx)
            
    return gff
