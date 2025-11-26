import re
from .general_funcs import progress_bar


def buld_gff3_line(line, attribute_dict):
    """Line should already be split by \t"""
    new_line = '\t'.join(line[:-1])
    
    attributes = ''
    for attribute in attribute_dict:
        attributes += attribute + '=' + attribute_dict[attribute] + ';'
    
    new_line += '\t' + attributes[:-1]
    
    return new_line


def extract_features(features):
    codons = []
    exons = []
    introns = []
    cds = []
    for feature in features:
        feature_type = feature[2]
        if feature_type == 'exon':
            exons.append(feature)
        elif feature_type == 'intron':
            introns.append(feature)
        elif feature_type == 'CDS':
            cds.append(feature)
        elif feature_type == 'mRNA':
            pass
        elif 'codon' in feature_type:
            codons.append((feature, feature_type))
        else:
            raise ValueError("Unknown feature: " + feature_type)
            
    if not len(codons) == 2:
        codons = None
            
    return codons, exons, introns, cds


def process_feature(feature, feature_type, transcript_name, gene_name, feature_number, gbkey):
    feature_name = transcript_name + f'_{feature_type}_' + feature_number
    feature_text = buld_gff3_line(feature, {
        'ID': feature_type + '-' + feature_name,
        'Parent': 'rna-' + transcript_name,
        'gbkey': gbkey,
        'gene': gene_name,
        'transcript_id': transcript_name
    })
    
    return feature_text


def gtf_to_gff(gtf, show_progress_bar = False):
    """
    NOTE: Designed to be used on the output of the BRAKER pipeline
    Reformats it in the stle of the Liftoff pipeline
    """
    if gtf[-1] == '\n':
        gtf = gtf[:-1]
    
    gtf = re.split(r'\n(.*\tgene\t)', '\n' + gtf)[1:]
    gtf = [gtf[x] + gtf[x + 1] for x in range(0, len(gtf), 2)]
    num_genes = len(gtf)
    gff_list = ['##gff-version 3\n#Converted from gtf\n']
    
    for idx, gene in enumerate(gtf):
        gene_text = ''
        if show_progress_bar:
            print('\r' + progress_bar(idx, num_genes), end='\r')
            
        gene_line = gene.split('\n', 1)[0].split('\t')
        gene_name = gene_line[-1]
        transcripts = re.split(r'\n(.*\ttranscript\t)', gene)
        transcripts = [transcripts[x] + transcripts[x + 1] for x in range(1, len(transcripts), 2)]
        
        gene_text += buld_gff3_line(gene_line, {
            'ID': 'gene-' + gene_name,
            'Name': gene_name,
            'gbkey': 'Gene',
            'gene': gene_name,
            'gene_biotype': 'protein_coding'
        }) + '\n'
        
        for transcript in transcripts:
            transcript = [x.split('\t') for x in transcript.split('\n')]
            transcript_line = transcript[0]
            transcript_name = transcript_line[-1]
            
            gene_text += buld_gff3_line(transcript_line, {
                'ID': 'rna-' + transcript_name, 
                'Parent': 'gene-' + gene_name, 
                'Name': transcript_name,
                'gbkey': 'mRNA',
                'gene': gene_name,
                'transcript_id': transcript_name
            }) + '\n'
            
            codons, exons, introns, cdss = extract_features(transcript[1:])
                
            if codons:
                gene_text += process_feature(codons[0][0], codons[0][1], transcript_name, gene_name, '', 'codon') + '\n'
            
            for idx, intron in enumerate(introns):
                gene_text += process_feature(intron, 'intron', transcript_name, gene_name, str(idx), 'intron') + '\n'
            
            for idx, exon in enumerate(exons):
                gene_text += process_feature(exon, 'exon', transcript_name, gene_name, str(idx), 'mRNA') + '\n'
                
            for idx, cds in enumerate(cdss):
                gene_text += process_feature(cds, 'cds', transcript_name, gene_name, str(idx), 'CDS') + '\n'
                
            if codons:
                gene_text += process_feature(codons[1][0], codons[1][1], transcript_name, gene_name, '', 'codon') + '\n'

        gff_list.append(gene_text)
    
    gff = ''.join(gff_list)
         
    if show_progress_bar:
        print('\r' + progress_bar(100, 100))
    
    return gff
