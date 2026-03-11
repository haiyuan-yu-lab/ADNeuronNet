import os
import sys
import pandas as pd
import random
from tool.enrichment import perform_gsea_human_enrich_dict
from tool.var_name import map_symbol_to_HGNC_to_NCBI, map_NCBI_to_HGNC_to_symbol, init_geneName_mapping
import tool.network_proximity as ppi

def sample_from_ppi(list_symbol, ppi, repeat=1, seed=42):
    """
    Perform degree-preserving sampling in a PPI network.
    """

    # load gene name mapping
    MappingGene = init_geneName_mapping()

    list_ncbi = [map_symbol_to_HGNC_to_NCBI(s, MappingGene) for s in list_symbol]
    print(f"Initial symbols: {len(list_ncbi)}")
    list_ncbi = [x for x in list_ncbi if x is not None]  # Remove None values
    print(f"Symbols after removing None: {len(list_ncbi)}")
    list_index = ppi.Name2Index(list_ncbi)
    print(f"PPI indices: {len(list_index)}")
    random.seed(seed)
    sampled_index = ppi.DegreePreserveNSampling(list_index, repeat)
    sampled_name = [ppi.Index2Name(s) for s in sampled_index]
    sampled_symbol = [
        [map_NCBI_to_HGNC_to_symbol(s, MappingGene) for s in sample_list]
        for sample_list in sampled_name
    ]
    sampled_symbol = [[x for x in sample_list if x is not None] for sample_list in sampled_symbol]
    return sampled_symbol

def read_gene_list(file_path):
    """
    Read a TSV file and extract the 'symbol' column as a list.
    """
    df = pd.read_csv(file_path, sep="\t", dtype=str)
    if 'symbol' not in df.columns:
        raise ValueError(f"'symbol' column not found in file: {file_path}")
    return df['symbol'].tolist()

def process_permutation_analysis(pathway_names, df_enrichment, path_output_dir):
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    Permutation_Variable = 'Odds Ratio'
    
    for pathway_name in pathway_names:
        # Filter the dataframe for the current pathway
        df_enrichment_filtered = df_enrichment[df_enrichment['Term'] == pathway_name]

        # save the filtered dataframe to a file
        output_file = os.path.join(path_output_dir, f"{pathway_name}_permutation_results.csv")
        df_enrichment_filtered.to_csv(output_file, index=False)
        
        if df_enrichment_filtered.empty:
            print(f"No data found for pathway: {pathway_name}")
            continue
        
        # Map query to adjusted p-value
        list_query = df_enrichment_filtered['query'].tolist()
        list_adj_p = df_enrichment_filtered[Permutation_Variable].tolist()
        dict_query_value = dict(zip(list_query, list_adj_p))
        
        # Extract observed and permuted values
        if 'observed' in dict_query_value.keys():
            observed_adj_p = dict_query_value.pop('observed')
        else:
            print(f"Observed value not found in the enrichment results for {pathway_name}")
            continue
        
        perm_adj_p = list(dict_query_value.values())
        
        # Calculate p-value
        p_value = sum([1 for i in perm_adj_p if i >= observed_adj_p]) / len(perm_adj_p)
        p_value = round(p_value, 3)

        # write p-value to f"{pathway_name}_enrichment_results.csv" as the first line, # p = 0.001
        with open(output_file, 'r') as f:
            data = f.read()
        with open(output_file, 'w') as f:
            f.write(f"# p = {p_value}\n")
            f.write(data)
        
        # Plot histogram
        plt.figure(figsize=(4, 3))
        sns.set(style="whitegrid")
        
        plt.hist(perm_adj_p, bins=10, color='blue', alpha=0.7)
        plt.axvline(x=observed_adj_p, color='red', linestyle='--')
        
        plt.xlabel(Permutation_Variable)
        plt.ylabel('Frequency')
        plt.title(f'{pathway_name} (p = {p_value})')
        
        # Save plot to a file
        output_file = os.path.join(path_output_dir, f"{pathway_name}_histogram.png")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        
        print(f"Histogram saved to {output_file}")


if __name__ == "__main__":
    # Ensure the correct number of arguments
    if len(sys.argv) < 8:
        print("Usage: python scr/GSEA_permutation_analysis.py <bait_interactors_dir> <bait_name> <endophenotype_gmt> <pathway_name> <output_dir> <random_seed> <Permutation_Number>")
        sys.exit(1)

    # Parse input arguments
    bait_interactors_dir = sys.argv[1]
    bait_name = sys.argv[2]
    endophenotype_gmt = sys.argv[3]
    pathway_names_file = sys.argv[4]
    output_dir = sys.argv[5]
    random_seed = int(sys.argv[6])
    permutation_number = int(sys.argv[7])

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Initialize PPI object
    ppi = ppi.Interactome(pathG='scr/tool/HumanInteractome.tsv')

    # Initialize GSEA dictionary
    dict_gsea = {}

    # Read bait interactors
    df_bait = pd.read_csv(os.path.join(bait_interactors_dir, f"{bait_name}.tsv"), sep="\t", dtype=str)
    querry_symbol_list = df_bait['symbol'].tolist()
    
    Permutation_Variable = 'Odds Ratio'

    # read pathway names in lines of pathway_names_file
    with open(pathway_names_file, 'r') as f:
        pathway_names = f.readlines()
        pathway_names = [x.strip() for x in pathway_names]
    print(f"Pathway names: {pathway_names}")

    # Perform GSEA for the endophenotype
    print(f"Running Permutation for {bait_name}...")
    sampled_symbol_list = sample_from_ppi(querry_symbol_list, ppi, repeat=permutation_number, seed=random_seed)
    # add observed list to the dictionary
    dict_gsea['observed'] = querry_symbol_list
    # add permutation data
    for i, sample in enumerate(sampled_symbol_list):
        dict_gsea[f'perm_{i}'] = sample

    # Perform GSEA
    df_enrichment = perform_gsea_human_enrich_dict(dict_gsea, endophenotype_gmt)

    # summarize permutation results
    process_permutation_analysis(pathway_names, df_enrichment, output_dir)

    