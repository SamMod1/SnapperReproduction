import pandas as pd


class Genome:
    def __init__(self, sequences):
        self.genome = sequences
        self.annotation = None
        
    def add_annotation(self, annotation_path):
        self.annotation = pd.read_csv(annotation_path)
        self.annotation['gene_name'] = self.annotation['TranscriptName'] + self.annotation['TranscriptPath']
        self.annotation = self.annotation[['gene_name', 'ScaffoldName', 'FromPosition', 'ToPosition']]
    
    def extract_genes(self):
        if not self.annotation:
            raise ValueError('No genome annotation has been loaded')
            
        whole_sequence = ''
        for sequence in self.genome.sequences:
            whole_sequence += self.genome.sequences[sequence].sequence
        
        for idx in range(len(self.annotation)):
            gene = self.annotation.iloc[idx]
            name = gene['gene_name']
            from_pos = gene['FromPosition']
            to_pos = gene['ToPosition']
            sequence = whole_sequence[from_pos]
    