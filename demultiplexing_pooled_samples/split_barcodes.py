import os
import pandas as pd
import sys

def split_barcodes(pool: int):
    """
    Splits barcodes from a TSV file into separate CSV files for each sample in a given pool.

    Parameters:
    pool (int): The pool number to process.

    Output:
    CSV files containing barcodes for each sample in the specified pool.
    """

    barcode_dir = f"/scratch/alpine/elathouwers@xsede.org/hgsoc/demultiplexing/results/vireo/pool{pool}/"
    df = pd.read_csv(barcode_dir + "donor_ids.tsv", delimiter="\t")
    df = df[["cell", "donor_id"]]

    samples = df["donor_id"].unique()
    
    for s in samples:
        df_sample = df[df["donor_id"]==s]
        barcode_lst = df_sample["cell"].reset_index(drop=True)
        print(f"Found {len(barcode_lst)} barcodes for sample {s}.")

        short_sample_name = s.split("_")[0]
        file_path = f"{os.getcwd()}/htseq-count/pool{pool}/{short_sample_name}"
        os.makedirs(file_path, exist_ok=True)
        barcode_lst.to_csv(file_path + "/barcodes_per_sample.csv", index=False, header=False)
    
    return 


if __name__ == "__main__":
    for pool in range(1, 11):
        print(f"Processing pool {pool}.")
        split_barcodes(pool)
        