import os
import sys
import glob
import pandas as pd

def get_low_gene_density_regions(karyotype_df, gff_file, window_size=100000):
    """
    Identifies regions with zero gene density.
    """
    print(f"  - Reading GFF file: {gff_file}")
    try:
        gff_cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        gff_df = pd.read_csv(gff_file, sep='\t', comment='#', header=None, names=gff_cols)
        gene_df = gff_df[gff_df['type'] == 'gene'].copy()
        gene_df['start'] = gene_df['start'].astype(int)
        gene_df['end'] = gene_df['end'].astype(int)
    except Exception as e:
        print(f"Error reading or parsing GFF file: {e}", file=sys.stderr)
        return {}

    lgd_regions = {}
    print("  - Calculating gene density and finding low-density regions...")
    for _, row in karyotype_df.iterrows():
        chrom = row['Chr']
        chrom_len = int(row['End'])
        
        chrom_genes = gene_df[gene_df['seqid'] == chrom]
        
        bins = range(1, chrom_len + window_size, window_size)
        gene_counts = [0] * (len(bins) - 1)

        for _, gene_row in chrom_genes.iterrows():
            bin_index = (gene_row['start'] - 1) // window_size
            if bin_index < len(gene_counts):
                gene_counts[bin_index] += 1
        
        # Find consecutive windows with zero genes
        in_lgd = False
        lgd_start = 0
        chrom_lgd_list = []
        for i, count in enumerate(gene_counts):
            if count == 0 and not in_lgd:
                in_lgd = True
                lgd_start = i * window_size + 1
            elif count > 0 and in_lgd:
                in_lgd = False
                lgd_end = i * window_size
                chrom_lgd_list.append((lgd_start, lgd_end))
        
        if in_lgd: # handles case where LGD extends to chromosome end
            lgd_end = chrom_len
            chrom_lgd_list.append((lgd_start, lgd_end))
            
        lgd_regions[chrom] = chrom_lgd_list

    return lgd_regions


def get_tr_regions(candidates_dir):
    """
    Parses all .candidate files to get high TR coverage regions.
    """
    tr_regions = {}
    search_path = os.path.join(candidates_dir, '*.candidate')
    candidate_files = glob.glob(search_path)
    
    print(f"  - Found {len(candidate_files)} candidate files in {candidates_dir}")

    for filepath in candidate_files:
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip() or line.startswith((' ', '\t')):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 6:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        if chrom not in tr_regions:
                            tr_regions[chrom] = []
                        tr_regions[chrom].append((start, end))
        except Exception as e:
            print(f"    - Warning: Could not process file {os.path.basename(filepath)}: {e}", file=sys.stderr)
            
    return tr_regions


def find_centromeres(karyotype_file, gff_file, candidates_dir, output_file):
    """
    Main logic to identify centromeres by intersecting LGD and TR regions.
    """
    try:
        karyotype_df = pd.read_csv(karyotype_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Karyotype file not found at {karyotype_file}", file=sys.stderr)
        return

    lgd_regions_map = get_low_gene_density_regions(karyotype_df, gff_file)
    tr_regions_map = get_tr_regions(candidates_dir)

    print("  - Intersecting LGD regions with TR regions...")
    final_centromeres = []

    for chrom in karyotype_df['Chr']:
        lgd_list = lgd_regions_map.get(chrom, [])
        tr_list = tr_regions_map.get(chrom, [])
        
        if not lgd_list or not tr_list:
            continue

        best_candidate = None
        max_overlap_span = 0

        # Find the best contiguous block of TRs inside the largest LGD region
        for lgd_start, lgd_end in lgd_list:
            
            # Find all TRs that overlap with this LGD region
            overlapping_trs = []
            for tr_start, tr_end in tr_list:
                # Check for overlap: (StartA < EndB) and (StartB < EndA)
                if tr_start < lgd_end and lgd_start < tr_end:
                    overlapping_trs.append((tr_start, tr_end))
            
            if not overlapping_trs:
                continue

            # Merge overlapping TRs into one continuous block
            # The final centromere is the min-start to max-end of these TRs,
            # clamped within the LGD boundaries.
            min_tr_start = min(s for s, e in overlapping_trs)
            max_tr_end = max(e for s, e in overlapping_trs)
            
            # Clamp the boundaries
            final_start = max(min_tr_start, lgd_start)
            final_end = min(max_tr_end, lgd_end)
            
            span = final_end - final_start
            if span > max_overlap_span:
                max_overlap_span = span
                best_candidate = {
                    "Chr": chrom,
                    "CE_start": final_start,
                    "CE_end": final_end
                }
    
        if best_candidate:
            final_centromeres.append(best_candidate)

    # Write output
    if final_centromeres:
        out_df = pd.DataFrame(final_centromeres)
        out_df.to_csv(output_file, sep='\t', index=False, header=["Chr", "CE_start", "CE_end"])
        print(f"Success! Centromere summary file saved to: {output_file}")
    else:
        # Create an empty file with header if no centromeres are found
        with open(output_file, 'w') as f_out:
            f_out.write("Chr\tCE_start\tCE_end\n")
        print("Warning: No centromeres could be identified. An empty summary file was created.")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python 02_find_centromeres.py <karyotype.txt> <annotation.gff> <candidates_dir> <output.txt>", file=sys.stderr)
        sys.exit(1)
    
    find_centromeres(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
