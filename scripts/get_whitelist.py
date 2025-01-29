import sys
import scanpy as sc

if __name__ == "__main__":
    anndata_path = sys.argv[1]
    output_path = sys.argv[2]

    # Load the Scanpy object
    adata = sc.read_h5ad(anndata_path)

    # Extract obs_names
    obs_names = adata.obs_names

    # Write obs_names to a text file
    with open(output_path, 'w') as f:
        for obs_name in obs_names:
            f.write(obs_name + '\n')