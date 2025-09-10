import sys
from itertools import groupby
import argparse

def find_telomeres_from_fasta(fasta_file, output_file, telotype='plant', search_len=10000, min_repeats=30):
    """
    Finds telomeric repeats at the ends of sequences in a FASTA file and
    creates a RIdeogram-compatible label file.
    """
    print(f"  - Reading FASTA: {fasta_file}")
    print(f"  - Searching for '{telotype}' type telomeres...")

    if telotype == 'plant':
        fwd_repeat = 'TTTAGGG'
        rev_repeat = 'CCCTAAA'
    elif telotype == 'animal':
        fwd_repeat = 'TTAGGG'
        rev_repeat = 'CCCTAA'
    else:
        print(f"Error: Invalid telotype '{telotype}'. Must be 'plant' or 'animal'.", file=sys.stderr)
        sys.exit(1)
    
    telomere_data = []

    try:
        with open(fasta_file, 'r') as f_in:
            faiter = (x[1] for x in groupby(f_in, lambda line: line.startswith(">")))
            
            for header in faiter:
                header_line = next(header).strip()
                chrom = header_line[1:].split()[0]
                seq = "".join(s.strip().upper() for s in next(faiter))
                seq_len = len(seq)
                
                # --- Check 5' end (head) ---
                head = seq[:search_len]
                fwd_count = head.count(fwd_repeat)
                rev_count = head.count(rev_repeat)
                
                if fwd_count + rev_count >= min_repeats:
                    # To define the boundary, find the last repeat instance in the first 2kb
                    boundary_region = seq[:2000]
                    last_pos_fwd = boundary_region.rfind(fwd_repeat)
                    last_pos_rev = boundary_region.rfind(rev_repeat)
                    
                    # The end of the telomere is the end of the last found repeat
                    end_pos = max(last_pos_fwd + len(fwd_repeat), last_pos_rev + len(rev_repeat))
                    if end_pos > 0:
                        telomere_data.append(['Telomere', 'box', chrom, 1, end_pos, '377eb8']) # Blue color

                # --- Check 3' end (tail) ---
                tail = seq[-search_len:]
                fwd_count = tail.count(fwd_repeat)
                rev_count = tail.count(rev_repeat)

                if fwd_count + rev_count >= min_repeats:
                    # Find the first repeat instance in the last 2kb
                    boundary_region = seq[-2000:]
                    first_pos_fwd = boundary_region.find(fwd_repeat)
                    first_pos_rev = boundary_region.find(rev_repeat)
                    
                    start_pos_in_region = -1
                    # Find the minimum start position if both are found
                    if first_pos_fwd != -1 and first_pos_rev != -1:
                        start_pos_in_region = min(first_pos_fwd, first_pos_rev)
                    elif first_pos_fwd != -1:
                        start_pos_in_region = first_pos_fwd
                    elif first_pos_rev != -1:
                        start_pos_in_region = first_pos_rev
                    
                    if start_pos_in_region != -1:
                        # Convert position relative to the boundary region to position relative to the whole chromosome
                        start_pos = (seq_len - 2000) + start_pos_in_region + 1
                        telomere_data.append(['Telomere', 'box', chrom, start_pos, seq_len, '377eb8']) # Blue color

        with open(output_file, 'w') as f_out:
            f_out.write("Type\tShape\tChr\tStart\tEnd\tcolor\n")
            for row in telomere_data:
                f_out.write("\t".join(map(str, row)) + "\n")
        
        print(f"  - Telomere label file created successfully: {output_file}")

    except FileNotFoundError:
        print(f"Error: FASTA file not found at {fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while processing telomeres: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find telomeric repeats from a FASTA file for RIdeogram.")
    parser.add_argument("fasta_file", help="Input genome FASTA file.")
    parser.add_argument("output_file", help="Output path for the telomere labels text file.")
    parser.add_argument("--telotype", type=str, default='plant', choices=['plant', 'animal'], help="Type of telomere repeat to search for ('plant' or 'animal'). Default: plant.")
    parser.add_argument("--search_len", type=int, default=10000, help="How many bases to search from each end of the chromosome. Default: 10000.")
    parser.add_argument("--min_repeats", type=int, default=30, help="Minimum number of repeat monomers required in the search window to call a telomere. Default: 30.")

    args = parser.parse_args()
    
    find_telomeres_from_fasta(args.fasta_file, args.output_file, args.telotype, args.search_len, args.min_repeats)

