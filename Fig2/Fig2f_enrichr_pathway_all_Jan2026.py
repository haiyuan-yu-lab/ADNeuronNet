import scipy.stats as stats
import numpy as np
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

def gene_network_enrichment_analysis(
    pathway_gene_list,
    interactor_list,
    pathway_name,
    background_genes=None  # Changed to None default
):
    """
    Perform enrichment analysis to check if a pathway is enriched in PPI interactors.
    Parameters:
    -----------
    pathway_gene_list : list
        Genes in the pathway
    interactor_list : list
        PPI interactor genes
    pathway_name : str
        Name of the pathway
    background_genes : int or None
        Number of background genes for hypergeometric test.
        If None, will be calculated from all pathways and interactors.
    Returns:
    --------
    pandas DataFrame with enrichment results
    """
    # Convert input lists to sets, removing duplicates
    pathway_genes = set(pathway_gene_list)
    interactor_genes = set(interactor_list)
    
    # Calculate overlap
    overlap_genes = pathway_genes.intersection(interactor_genes)
    overlap = len(overlap_genes)
    
    # Hypergeometric test
    p_value = stats.hypergeom.sf(
        overlap - 1,
        background_genes,  # total population
        len(interactor_genes),  # number of interactors
        len(pathway_genes)  # size of pathway
    )
    
    # Store results
    enrichment_result = {
        'Pathway': pathway_name,
        'Overlap_Genes_Count': overlap,
        'Total_Interactors': len(interactor_genes),
        'Pathway_Size': len(pathway_genes),
        'P_Value': p_value,
        'Enrichment_Ratio': overlap / len(pathway_genes) if len(pathway_genes) > 0 else 0,
        'Overlapping_Genes': list(overlap_genes),
        'Background_Size': background_genes  # Add this for transparency
    }
    
    # Convert to DataFrame
    results_df = pd.DataFrame([enrichment_result])
    
    # Multiple testing correction will be done on the combined results
    return results_df

def load_pwy(df_pwy, pwy_column):
    """
    Loads pathway from files.
    Returns:
    A dictionary of pathways
    """
    pwy_dic = {}
    pwy_list = df_pwy[pwy_column].unique()
    
    for pwy in pwy_list:
        # Only include genes from the specific pathway
        pwy_dic[pwy] = df_pwy[df_pwy[pwy_column] == pwy]["Gene_Symbol"].tolist()
    
    return pwy_dic

def main():
    # --- File Paths and Constants ---
    PATH = "/Users/houy2/Documents/projects/AD/AD_PPI/haiyuan_ppi_v6_Jan2026"
    enrichr_dir = "pathway_enrichr"
    UpDate="Jan2026"
    
    # --- 1. Load PPI data ---
    df_nt = pd.read_csv(f"{PATH}/20260115_ADNeuronNet.csv", sep=",").rename(columns={"Genes":"Interactors"})
    Interactors_list = df_nt["Interactors"].unique().tolist()
    baits_list = df_nt["Bait"].unique().tolist()
    
    print(f"Loaded {len(Interactors_list)} unique Interactors")
    print(f"Loaded {len(baits_list)} unique baits")
    print("=====================")
    
    # --- 2. Load Pathways ---
    df_pwy = pd.read_csv(f"{PATH}/{enrichr_dir}/hsa00001.cleaned_Mar182025_withoutDis.tsv", sep="\t")
    hs_pwy_dic = load_pwy(df_pwy, "Pathway")
    
    # --- Calculate background gene count ---
    # Use all unique genes from both pathways and interactors
    all_pathway_genes = set()
    for genes in hs_pwy_dic.values():
        all_pathway_genes.update(genes)
    
    all_genes = all_pathway_genes.copy()
    all_genes.update(Interactors_list)
    background_genes = len(all_genes)
    
    print(f"Using {background_genes} genes as background")
    print("=====================")
    
    # --- 3. Pathway Enrichment Analysis ---
    all_enrich_results = []
    
    for pwy_name, pwy_genes in hs_pwy_dic.items():
        # Perform enrichment analysis
        pathway_results = gene_network_enrichment_analysis(
            pwy_genes, Interactors_list, pwy_name, background_genes
        )
        all_enrich_results.append(pathway_results)
    
    # --- 4. Combine Results ---
    combined_df = pd.concat(all_enrich_results, ignore_index=True)
    
    # --- 5. Multiple Testing Correction ---
    if not combined_df.empty:
        _, adjusted_p_values, _, _ = multipletests(
            combined_df['P_Value'],
            method='fdr_bh'
        )
        combined_df['Adjusted_P_Value'] = adjusted_p_values
    else:
        combined_df['Adjusted_P_Value'] = []
    
    # Sort by adjusted p-value
    combined_df = combined_df.sort_values('Adjusted_P_Value')
    
    # --- 6. Save Results ---
    output_file = f"{PATH}/{enrichr_dir}/PPI_pathway_enrichment_all_interactors_{UpDate}.tsv"
    combined_df.to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to: {output_file}")

if __name__ == '__main__':
    main()