This Python script identifies candidate genes by integrating genomic data. It reads a file containing regions of high Fst (a measure of genetic differentiation)
and a separate gene map file that links genomic windows to gene IDs. The script then filters these high Fst regions to include only those on specified chromosomes (10 and 11), 
retrieves the associated gene IDs from the map, and outputs a list of these potential candidate genes along with their corresponding Fst values and genomic locations.

import csv
import os

# --- Configuration ---
# Define file paths
TOP_FST_FILE = "/lustre/miifs01/project/m2_jgu-salmosex/nina/results_updated/fst_dxy_pi/fst_dxy_pi_org/top1pct_fst_windows_ORG.csv"
GENE_MAP_FILE = "/lustre/miifs01/project/m2_jgu-salmosex/nina/results_updated/window_overlap/Bter1.0_window_gene_ids.txt"
OUTPUT_FILE = "/lustre/miifs01/project/m2_jgu-salmosex/nina/results_updated/window_overlap/candidate_genes_chr10_11.txt"

# Define chromosome IDs for Bter1.0 (chr10 & chr11)
CHR10 = "NC_015771.1"
CHR11 = "NC_015772.1"

# --- Main Logic ---
def extract_candidate_genes():
    """
    Extracts candidate genes by overlapping high Fst windows with a gene map.
    """
    gene_map = {} # This will store our gene IDs, with window as key

    # 1. Load the gene map file
    try:
        with open(GENE_MAP_FILE, 'r') as f_gene:
            for line_num, line in enumerate(f_gene, start=1):
                line = line.strip() # Remove leading/trailing whitespace and newline
                if not line: # Skip empty lines
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    window_key = parts[0].strip() # Key is the window (e.g., "NC_015778.1:550000-560000")
                    gene_id = parts[1].strip()   # Value is the gene ID(s)
                    gene_map[window_key] = gene_id
                else:
                    print(f"Warning: Skipping malformed line {line_num} in gene map: {line}", file=os.sys.stderr)
        print(f"Loaded {len(gene_map)} entries from gene map.", file=os.sys.stderr)

        # DEBUG: Print a few sample keys from the gene_map, especially for target chromosomes
        print("\n--- Gene Map Sample Keys (relevant to target chromosomes) ---", file=os.sys.stderr)
        sample_count = 0
        for key in gene_map:
            if CHR10 in key or CHR11 in key:
                print(f"Gene Map Key: {key}", file=os.sys.stderr)
                sample_count += 1
                if sample_count >= 5: # Print up to 5 relevant keys
                    break
        if sample_count == 0:
            print(f"No gene map keys found for {CHR10} or {CHR11}.", file=os.sys.stderr)
        print("----------------------------------------------------------\n", file=os.sys.stderr)


    except FileNotFoundError:
        print(f"Error: Gene map file not found at {GENE_MAP_FILE}", file=os.sys.stderr)
        return
    except Exception as e:
        print(f"Error loading gene map file: {e}", file=os.sys.stderr)
        return

    # 2. Process the TOP_FST file and write to output
    try:
        with open(OUTPUT_FILE, 'w') as f_out:
            # Write header to output file
            f_out.write("Chromosome\tWindow\tFst\tGenes\n")

            with open(TOP_FST_FILE, 'r', newline='') as f_fst:
                csv_reader = csv.reader(f_fst)
                header = next(csv_reader) # Read and skip the header line

                for row_num, row in enumerate(csv_reader, start=2): # Start row_num at 2 for data lines
                    if len(row) < 6:
                        print(f"Warning: Skipping malformed row {row_num} in TOP_FST: {row}", file=os.sys.stderr)
                        continue

                    try:
                        chr_val = row[2].strip()
                        # Convert start_val to integer, subtract 1 for 0-based, then convert back to string
                        start_val = str(int(row[3].strip()) - 1)
                        end_val = row[4].strip()
                        fst_val = row[5].strip()
                    except (ValueError, IndexError) as e:
                        print(f"Error processing columns in row {row_num} of TOP_FST: {row} - {e}", file=os.sys.stderr)
                        continue

                    if chr_val == CHR10 or chr_val == CHR11:
                        # Construct the window key with the adjusted start coordinate
                        window_key = f"{chr_val}:{start_val}-{end_val}"
                        # DEBUG: Print the constructed window key for each matching row
                        print(f"DEBUG: Constructed key from TOP_FST: {window_key}", file=os.sys.stderr)

                        genes = gene_map.get(window_key, "NA")

                        f_out.write(f"{chr_val}\t{window_key}\t{fst_val}\t{genes}\n")

    except FileNotFoundError:
        print(f"Error: TOP_FST file not found at {TOP_FST_FILE}", file=os.sys.stderr)
    except Exception as e:
        print(f"Error processing TOP_FST file: {e}", file=os.sys.stderr)

# Run the extraction function
if __name__ == "__main__":
    extract_candidate_genes()

