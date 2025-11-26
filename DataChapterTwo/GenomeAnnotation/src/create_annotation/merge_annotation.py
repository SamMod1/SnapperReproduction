from .gff import Gff
from .general_funcs import progress_bar, cds_overlap_clusters


def process_gene_subcluster(gff, chosen_genes, cluster):
    if len(cluster) == 1: # i.e. gene is alone
        chosen_genes.extend(gff[cluster[0]])
    
    else:
        liftoff_genes = []
        braker_genes = []
        
        for gene_name in cluster:
            gene = gff.get_feature(gene_name)
            if gene['source'] == 'Liftoff':
                liftoff_genes.append(gene_name)
            else:
                braker_genes.append(gene_name)
        
        cluster_str = str(liftoff_genes + braker_genes).replace('[', '').replace(']', '').replace("'", "").replace(' ', '')
        
        for gene_name in braker_genes:
            gene = gff[gene_name]
            gene.edit_feature_data(gene_name, 'cluster', cluster_str)
            chosen_genes.extend(gene)  # Lets just keep the Braker genes for now...


def merge_annotation(gff, gene_clusters, show_progress_bar = False):
    chosen_genes = Gff()
    
    for idx, sequence in enumerate(gene_clusters):
        if show_progress_bar:
            print('\rDeciding which genes to keep: ' + progress_bar(idx, len(gene_clusters)), end='\r')
            
        for cluster in gene_clusters[sequence]:
            cluster_gff = gff[cluster]
            subclusters = cds_overlap_clusters(cluster_gff, [x for x in cluster_gff if x[:5] == 'gene-'])
            
            for cluster in subclusters:
                process_gene_subcluster(gff, chosen_genes, cluster)
                
    if show_progress_bar:
        print('\rDeciding which genes to keep: ' + progress_bar(100, 100), end='\r')

    return chosen_genes
                