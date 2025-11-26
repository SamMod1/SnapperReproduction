from SamBio import MultiSequence


TRANSCRIPT_FEATURES = ['mRNA', 'transcript']


class FastaLink:
    def __init__(self, gff, fasta):
        self.gff = gff
        self.fasta = fasta
        
    def __getitem__(self, feature_name):
        feature = self.gff.get_feature(feature_name)
        
        sequence = self.fasta[feature['sequence']][feature['start'] - 1:feature['end']]
        
        return sequence
        
    def get_codingsequence(self, transcript_name):
        transcript = self.gff[transcript_name]
        transcript_data = self.gff.get_feature(transcript_name)
        
        transcript_df = transcript.to_df()
        transcript_df = transcript_df.loc[transcript_df['feature_type'] == 'CDS']
        
        if transcript_df.empty:
            raise ValueError("The given transcript/feature had no coding sequences")
        
        strand = transcript_df['strand'].unique()
        
        if len(strand) > 1:
            raise ValueError("The given transcript/feature had CDS subfeatures on multiple different strands")
        
        strand = strand[0]
        
        if strand == '-':
            transcript_df = transcript_df.sort_values(by='start', ascending=False)
        else:
            transcript_df = transcript_df.sort_values('start')
            
        row = transcript_df.iloc[0]
        if strand == '-':
            sequence = self.fasta[transcript_data['sequence']][row['start'] - 1:row['end']].compliment()
            for idx in range(1, len(transcript_df)):
                row = transcript_df.iloc[idx]
                sequence += self.fasta[transcript_data['sequence']][row['start'] - 1:row['end']].compliment()
        else:
            sequence = self.fasta[transcript_data['sequence']][row['start'] - 1:row['end']]
            for idx in range(1, len(transcript_df)):
                row = transcript_df.iloc[idx]
                sequence += self.fasta[transcript_data['sequence']][row['start'] - 1:row['end']]
        
        sequence.name = transcript_name
        
        return sequence
    
    def get_protein(self, transcript_name, offset = 0):
        codingseq = self.get_codingsequence(transcript_name)
        return codingseq.to_protein(offset = offset)
    
    def create_proteome(self, alt_transcripts = True, as_cds = False):
        proteome = MultiSequence()
        proteins_to_add = []
        protein = None
        
        for feature in self.gff:
            feature = self.gff.get_feature(feature)
            if not alt_transcripts and feature['feature_type'] == 'gene':
                if proteins_to_add:
                    protein_lengths = [len(x) for x in proteins_to_add]
                    longest_protein = max(protein_lengths)
                    protein = proteins_to_add[protein_lengths.index(longest_protein)]
                    proteome[protein.name] = protein
                    proteins_to_add = []
                
            elif feature['feature_type'] in TRANSCRIPT_FEATURES:
                try:
                    if not as_cds:
                        protein = self.get_protein(feature['ID'])
                    else:
                        protein = self.get_codingsequence(feature['ID'])
                except ValueError:
                    pass
            
                if protein:
                    if not alt_transcripts:
                        proteins_to_add.append(protein)
                    else:
                        proteome[protein.name] = protein
            protein = None
        
        if proteins_to_add:
            protein_lengths = [len(x) for x in proteins_to_add]
            longest_protein = max(protein_lengths)
            protein = proteins_to_add[protein_lengths.index(longest_protein)]
                
        return proteome
            