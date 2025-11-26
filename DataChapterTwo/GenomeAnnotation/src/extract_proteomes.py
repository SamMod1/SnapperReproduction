import sys


ORTHODB_PARTITION = sys.argv[1]
SPECIES = sys.argv[2]
OUTPUT_DIR = sys.argv[3]


def main(orthodb_partition, species, outdir):
    if not outdir[-1] == '/':
        outdir += '/'
    proteins = open(orthodb_partition, 'r').read()
    species_to_extract = open(species, 'r').read().split('\n')

    species_ids = {}

    proteomes = {}

    for line in species_to_extract:
        line = line.split('\t', 1)
        species_ids[line[1]] = line[0]
        proteomes[line[1]] = ''

    proteins = proteins.split('>')[1:]
    
    for protein in proteins:
        name = protein.split('\n', 1)[0]
        species = name.split(':', 1)[0]

        if species in species_ids:
            proteomes[species] += '>' + protein

    for species in species_ids:
        open(outdir + species_ids[species] + '.pep', 'w').write(proteomes[species])


if __name__ == '__main__':
    main(ORTHODB_PARTITION, SPECIES, OUTPUT_DIR)
