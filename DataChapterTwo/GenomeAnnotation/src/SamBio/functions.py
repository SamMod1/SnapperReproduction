from .Sequence import Sequence
from .MultiSequence import MultiSequence



def read_fasta(filename):
    text = open(filename, 'r').read()

    n_sequences = len(text.split('>'))
    
    if n_sequences > 2:
        return MultiSequence(text)
    else:
        return Sequence(text)


def _process_feature(feature, genome):
    if feature:
        if not feature[0] == '#':
            fields = feature.split('\t')
            if len(fields) > 8:
                if fields[2] == 'CDS':
                    phase = int(fields[7])
                    if fields[6] == '+':
                        return genome[fields[0]][int(fields[3]) - 1:int(fields[4])].to_protein(phase)
                    else:
                        return genome[fields[0]][int(fields[3]) - 1:int(fields[4])].compliment().to_protein(phase)
            return

def genome_to_proteome(genome, gff_file):
    proteome = MultiSequence()
    file = open(gff_file, 'r').read()
    genes = file.split('\tmRNA\t')[1:]
    for gene in genes:
        gene_name = gene.split('\t')[5].split(';')[0][3:]
        gene_sequence = Sequence(name = gene_name)
        features = gene.split('\n')
        for feature in features:
            feature_sequence = _process_feature(feature, genome)
            if not feature_sequence is None:
                gene_sequence += feature_sequence
        proteome[gene_name] = gene_sequence
    return proteome
