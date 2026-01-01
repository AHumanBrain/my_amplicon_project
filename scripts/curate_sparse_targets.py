import sys
import random

def curate_sparse_targets(gff_file, output_file, count=100, min_dist=2500):
    genes = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            cols = line.split('\t')
            if len(cols) < 9: continue
            if cols[2] == 'gene':
                start = int(cols[3])
                end = int(cols[4])
                contig = cols[0]
                info = cols[8]
                gene_id = None
                if 'ID=' in info: gene_id = info.split('ID=')[1].split(';')[0]
                elif 'gene=' in info: gene_id = info.split('gene=')[1].split(';')[0]
                
                # Filter out very short genes to avoid the failures we saw earlier
                if (end - start) < 150: continue
                
                if gene_id:
                    # Clean the ID (remove 'gene-' prefix if present)
                    clean_id = gene_id.replace('gene-', '')
                    genes.append({'id': clean_id, 'contig': contig, 'start': start, 'end': end})

    # Sort by position
    genes.sort(key=lambda x: (x['contig'], x['start']))
    
    selected = []
    last_end = -min_dist
    last_contig = None
    
    # Shuffle slightly to get various genes but keep it mostly ordered for distance checking
    random.shuffle(genes)
    genes.sort(key=lambda x: (x['contig'], x['start']))

    for g in genes:
        if len(selected) >= count: break
        
        if g['contig'] != last_contig or g['start'] > (last_end + min_dist):
            selected.append(g['id'])
            last_end = g['end']
            last_contig = g['contig']
            
    if len(selected) < count:
        print(f"Warning: Only found {len(selected)} sparse targets.")
    
    with open(output_file, 'w') as f:
        for gid in selected:
            f.write(f"{gid}\n")
    
    print(f"Successfully curated {len(selected)} sparse targets into {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python curate_sparse_targets.py <gff_file> <output_file> [count] [min_dist]")
    else:
        cnt = int(sys.argv[3]) if len(sys.argv) > 3 else 100
        dist = int(sys.argv[4]) if len(sys.argv) > 4 else 2500
        curate_sparse_targets(sys.argv[1], sys.argv[2], cnt, dist)
