import sys

gff_file = sys.argv[1]
count = int(sys.argv[2])
output_file = sys.argv[3]

genes = []
with open(gff_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.split('\t')
        if len(parts) > 8 and parts[2] == 'gene':
            # Extract ID or Name
            attributes = parts[8]
            gene_id = None
            if 'Name=' in attributes:
                gene_id = attributes.split('Name=')[1].split(';')[0]
            elif 'ID=' in attributes:
                gene_id = attributes.split('ID=')[1].split(';')[0]
            
            if gene_id and gene_id not in genes:
                genes.append(gene_id)
            
            if len(genes) >= count:
                break

with open(output_file, 'w') as f:
    for gene in genes:
        f.write(f"{gene}\n")

print(f"Extracted {len(genes)} genes to {output_file}")
