import sys

def generate_karyotype_from_fasta(fasta_file, output_file):
    """
    Reads a FASTA file and generates a RIdeogram-compatible karyotype file.
    Format: Chr, Start, End
    """
    karyotype_data = []
    current_chrom = None
    current_len = 0

    try:
        with open(fasta_file, 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if line.startswith('>'):
                    if current_chrom is not None:
                        karyotype_data.append((current_chrom, 0, current_len))
                    
                    current_chrom = line[1:].split()[0]
                    current_len = 0
                else:
                    current_len += len(line)
            
            if current_chrom is not None:
                karyotype_data.append((current_chrom, 0, current_len))

        with open(output_file, 'w') as f_out:
            f_out.write("Chr\tStart\tEnd\n")
            for chrom, start, end in karyotype_data:
                f_out.write(f"{chrom}\t{start}\t{end}\n")
        
        print(f"Success! Base karyotype file saved to: {output_file}")

    except FileNotFoundError:
        print(f"Error: File not found {fasta_file}. Please check the path.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 01_generate_karyotype.py <genome.fasta> <output_karyotype.txt>", file=sys.stderr)
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_txt = sys.argv[2]
    
    generate_karyotype_from_fasta(input_fasta, output_txt)
