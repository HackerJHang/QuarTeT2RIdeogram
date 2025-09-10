import sys
import re
from itertools import groupby
import argparse

def find_gaps_from_fasta(fasta_file, output_file, min_gap_len=10):
    """
    Finds assembly gaps (stretches of 'N's) from a FASTA file and
    creates a RIdeogram-compatible label file.
    """
    print(f"  - Reading FASTA: {fasta_file}")
    print(f"  - Looking for gaps of at least {min_gap_len} 'N's.")
    
    try:
        with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
            f_out.write("Type\tShape\tChr\tStart\tEnd\tcolor\n")
            
            # Efficiently parse FASTA file
            faiter = (x[1] for x in groupby(f_in, lambda line: line.startswith(">")))
            
            for header in faiter:
                header_line = next(header).strip()
                chrom = header_line[1:].split()[0]
                
                # Join sequence lines for the current chromosome
                seq = "".join(s.strip().upper() for s in next(faiter))
                
                # Use regex to find all long stretches of 'N's
                for match in re.finditer(f'N{{{min_gap_len},}}', seq):
                    start = match.start() + 1 # Convert to 1-based coordinate
                    end = match.end()
                    f_out.write(f"Gap\tbox\t{chrom}\t{start}\t{end}\te41a1c\n") # Red color for gaps
        
        print(f"  - Gap label file created successfully: {output_file}")

    except FileNotFoundError:
        print(f"Error: FASTA file not found at {fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while processing gaps: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find assembly gaps from a FASTA file for RIdeogram.")
    parser.add_argument("fasta_file", help="Input genome FASTA file.")
    parser.add_argument("output_file", help="Output path for the gap labels text file.")
    parser.add_argument("--min_gap_len", type=int, default=10, help="Minimum number of consecutive 'N's to be considered a gap. Default: 10.")
    
    args = parser.parse_args()
    
    find_gaps_from_fasta(args.fasta_file, args.output_file, args.min_gap_len)

