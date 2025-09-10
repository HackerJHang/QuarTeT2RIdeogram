import sys
import pandas as pd

def merge_files(base_karyotype_path, centromere_summary_path, output_path):
    """
    Merges the base karyotype file with the centromere summary file using pandas.
    """
    try:
        base_df = pd.read_csv(base_karyotype_path, sep='\t')
        
        # Check if centromere file is empty or just has a header
        try:
            centro_df = pd.read_csv(centromere_summary_path, sep='\t')
            if centro_df.empty:
                 print("  - Centromere summary file is empty, proceeding without centromere data.")
                 centro_df = pd.DataFrame(columns=['Chr', 'CE_start', 'CE_end'])
        except pd.errors.EmptyDataError:
            print("  - Centromere summary file is empty, proceeding without centromere data.")
            centro_df = pd.DataFrame(columns=['Chr', 'CE_start', 'CE_end'])

        merged_df = pd.merge(base_df, centro_df, on='Chr', how='left')
        
        # Ensure column order is correct for RIdeogram
        final_cols = ['Chr', 'Start', 'End', 'CE_start', 'CE_end']
        merged_df = merged_df.reindex(columns=final_cols)

        merged_df.to_csv(output_path, sep='\t', index=False, na_rep='')
        
        print(f"Success! Final karyotype file saved to: {output_path}")

    except FileNotFoundError as e:
        print(f"Error: File not found {e.filename}. Please check the path.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during merging: {e}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 05_merge_karyotype.py <base_karyotype.txt> <centromere_summary.txt> <output_final_karyotype.txt>", file=sys.stderr)
        sys.exit(1)
        
    merge_files(sys.argv[1], sys.argv[2], sys.argv[3])
