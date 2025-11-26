import sys

genome = sys.argv[1]
output_gnome = sys.argv[2]
output_mappings = sys.argv[3]

genome = open(genome, 'r').read().split('>')

mappings = {}
new_genome = ''

counts = {'scaffold': 1, 'contig': 1}

for x in genome:
    if x:
        header, sequence = x.split('\n', 1)

        if 'scaffold' in header:
            new_header = 'scaffold' + str(counts['scaffold'])
            mappings[header] = new_header
            counts['scaffold'] += 1

        elif 'contig' in header:
            new_header = 'contig' + str(counts['contig'])
            mappings[header] = new_header
            counts['contig'] += 1

        else:
            new_header = header
            mappings[header] = header

        new_genome += '>' + new_header + '\n' + sequence + '\n'
    
        
mappings_txt = ''

for key in mappings:
    mappings_txt += key + '\t' + mappings[key] + '\n'
    
open(output_gnome, 'w').write(new_genome)
open(output_mappings, 'w').write(mappings_txt)
        