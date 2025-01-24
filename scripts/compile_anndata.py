import sys
import scanpy as sc
import anndata as an
import pandas as pd
import scipy.sparse as sparse

def combine_and_annotate_h5ad(h5ad_paths, gene_table_path, outpath):
    """
    Combines multiple h5ad files, merges gene metadata, and saves the result.

    This function takes a list of h5ad files (presumably containing count data
    for different chromosomes or genomic regions), concatenates them into a 
    single AnnData object. It then loads a gene metadata table and merges it 
    with the AnnData object's `var` attribute. Finally, it removes unnamed 
    genes and saves the resulting AnnData object to a specified path.

    Args:
        h5ad_paths (list): A list of paths to the input h5ad files.
        gene_table_path (str): The path to the gene metadata table (TSV format).
        outpath (str): The path where the combined and annotated AnnData object
                       will be saved.

    Raises:
        FileNotFoundError: If any of the input h5ad files or the gene table 
                           file does not exist.
        ValueError: If the gene table is not in the expected tab-separated 
                    format or if there are issues merging the data.
    """

    try:
        # Combine counts from each chromosome/region
        print("Reading and combining h5ad files...")
        adata_list = [sc.read_h5ad(x) for x in h5ad_paths]
        adata = an.concat(adata_list, axis=1, join="outer", merge='only')

        # Load the gene metadata
        print("Loading gene metadata...")
        df = pd.read_csv(gene_table_path, sep='\t', index_col=0) # Assuming gene names or IDs are in the first column as the index

        # Merge the gene metadata
        print("Merging gene metadata...")
        adata.var = pd.merge(
            adata.var, df, how='left',
            left_index=True,
            right_index=True,
        )

        # Drop unnamed genes
        print("Removing unnamed genes...")
        adata = adata[:, adata.var['gene_name'].notna()].copy()

        # Ensure sparse matrix
        if not sparse.issparse(adata.X):
            print("Converting to sparse matrix (CSR format)...")
            adata.X = sparse.csr_matrix(adata.X)

        # Save the combined AnnData object
        print(f"Saving combined AnnData object to {outpath}...")
        adata.write_h5ad(outpath)
        print("Done!")

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: Data processing issue - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script_name.py <outpath> <gene_table_path> <h5ad_path1> <h5ad_path2> ...")
        sys.exit(1)

    outpath = sys.argv[1]
    gene_table_path = sys.argv[2]
    h5ad_paths = sys.argv[3:]

    combine_and_annotate_h5ad(h5ad_paths, gene_table_path, outpath)