import pandas as pd
import numpy as np
import json
import gzip
import re
import os
# from pathlib import Path
import pathlib
from typing import Dict, Tuple
from Bio.PDB import MMCIFParser
import Bio
from scipy.spatial.distance import pdist

def _load_af3_files(
        res_dir: pathlib.Path,
        model_idx: int) -> Tuple[Bio.PDB.Structure.Structure,
                                 Dict,
                                 Dict]:
    """
    Loads an AF3 model into memory

    Parameters
    ----------
    res_dir : Path
        Path to the directory containing .cif and .json files, downloaded from
        AF3 server
    model_idx : int
        0-based index of the models

    Returns
    -------
    structure : Biopython structure
    full_data : dictionaryk
        parsed from the *_full_data_<model_idx>.json
    summary : dictionary
        parsed from the *_summary_confidences_<model_idx>.json
    """
    options = list(res_dir.glob(f"*_model_{model_idx}.cif"))
    if len(options) != 1:
        raise ValueError(
            f"Model {model_idx} has {len(options)} cif candidates")
    cif_file = options[0]
    options = list(res_dir.glob(f"*_full_data_{model_idx}.json"))
    if len(options) != 1:
        raise ValueError(
            f"Model {model_idx} has {len(options)} full_data candidates")
    full_data_file = options[0]
    options = list(res_dir.glob(f"*_summary_confidences_{model_idx}.json"))
    if len(options) != 1:
        raise ValueError(
            f"Model {model_idx} has {len(options)} summary candidates")
    summary_file = options[0]
    parser = MMCIFParser()
    structure = parser.get_structure("AlphaFold", cif_file)
    full_data = json.load(full_data_file.open())
    summary = json.load(summary_file.open())
    return structure, full_data, summary


def calculate_mutation_residue_ca_distances(structure, protein1_chain: str, protein2_chain: str, mutation_residue_number: int, plddt_threshold: float = 25):
    """
    Calculate the minimum distance between the alpha carbon (Cα) of a mutated residue in protein1 
    and all residues in protein2, considering only atoms with pLDDT above a threshold.
    
    Parameters:
    -----------
    structure : Bio.PDB.Structure
        The protein structure loaded from the PDB file
    protein1_chain : str
        Chain identifier for protein1
    protein2_chain : str
        Chain identifier for protein2
    mutation_residue_number : int
        Residue number of the mutated residue in protein1
    plddt_threshold : float, optional
        Minimum pLDDT score to consider an atom (default: 25)
    
    Returns:
    --------
    distances : dict
        Dictionary with keys as residue numbers in protein2 and values as distances to the mutated residue in protein1
    """
    # Find the Cα atom of the mutated residue with pLDDT threshold
    mutated_residue = None
    for residue in structure[0][protein1_chain]:
        if residue.id[1] == mutation_residue_number:
            mutated_residue = residue
            break
    
    if mutated_residue is None or 'CA' not in mutated_residue or mutated_residue['CA'].bfactor < plddt_threshold:
        print(f"Mutation residue {mutation_residue_number} in chain {protein1_chain} does not have a valid Cα atom above the pLDDT threshold.")
        return None
        # raise ValueError(f"Mutation residue {mutation_residue_number} in chain {protein1_chain} does not have a valid Cα atom above the pLDDT threshold.")
    
    mutated_ca = mutated_residue['CA'].coord
    
    # Find Cα atoms of all residues in protein2 and compute distances
    distances = {}
    for residue in structure[0][protein2_chain]:
        if 'CA' in residue and residue['CA'].bfactor >= plddt_threshold:
            res_ca = residue['CA'].coord
            distance = np.linalg.norm(mutated_ca - res_ca)
            distances[residue.id[1]] = distance
    
    if not distances:
        print(f"No valid Calpha atoms found in chain {protein2_chain} above the pLDDT threshold.")
        return None
        # raise ValueError(f"No valid Calpha atoms found in chain {protein2_chain} above the pLDDT threshold.")
    
    return distances

if __name__ == "__main__":
    root_dir = "/fs/cbsuhy02/storage/tx82/"

    project_name = "final_interactome_noAPOEvar_20251210" # given the name of the project, which is the same as the name of the input file for AF3 and the output folder for AF3 results
    # load interactome 
    PPIs = pd.read_table(f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/{project_name}.csv", sep=',')
    PPIs["af3_job_name"] = PPIs.apply(lambda row: (row["Bait ID"]+"_"+row["Interactor ID"]).replace("-", "_").lower(), axis=1)
    # load af3 input uniprot IDs and sequences
    PPIs_seq = pd.read_table(f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/{project_name}_seq.txt", sep='\t')
    # load AD missense mutations
    AD_mutations_all = pd.read_table("/fs/cbsuhy02/storage/tx82/AD_interactome/results/AD_mutations_results/AD_mutations_all_AF3_updated.csv", sep=',')

    # define isoforms names 
    isoforms = ["o00499_9", "p49768_2", "O75530-2", "Q05084-3", "Q9H7Z6-2", "Q06413-5", "Q13492-3", "P49768-2", "Q9H0F6-2", "Q8NG11-2"]
    isoforms = [isoform.replace("-", "_").lower() for isoform in isoforms]

    # delete model 1/2/3/4, only keep model 0, for all jobs 
    for date_worker in os.listdir(f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/af3_output/"):
        af3_results_dir = f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/af3_output/{date_worker}/"
        if os.path.exists(af3_results_dir) and os.path.isdir(af3_results_dir):
            for af3_job_name in os.listdir(af3_results_dir):
                af3_job_results_dir = pathlib.Path(os.path.join(af3_results_dir, af3_job_name))
                # af3_job_results_dir = os.path.join(af3_results_dir, af3_job_name)
                if os.path.isdir(af3_job_results_dir):
                    # rm msas/
                    if os.path.exists(os.path.join(af3_job_results_dir, "msas")):
                        os.system(f"rm -r {os.path.join(af3_job_results_dir, 'msas')}")
                    
                    # rm templates/
                    if os.path.exists(os.path.join(af3_job_results_dir, "templates")):
                        os.system(f"rm -r {os.path.join(af3_job_results_dir, 'templates')}")
                    
                    # rm job_request
                    if os.path.exists(os.path.join(af3_job_results_dir, f"fold_{af3_job_name}_job_request.json")):
                        os.system(f"rm -r {os.path.join(af3_job_results_dir, 'fold_'+af3_job_name+'_job_request.json')}")

                    # rm cif files for model 1/2/3/4
                    options = list(af3_job_results_dir.glob(f"*_model_*.cif"))
                    for option in options:
                        if str(option)[-6:] != "_0.cif":
                            print(str(option))
                            os.system(f"rm {str(option)}")
                    
                    # rm full data json files for model 1/2/3/4
                    options = list(af3_job_results_dir.glob(f"*_full_data_*.json"))
                    for option in options:
                        if str(option)[-7:] != "_0.json":
                            print(str(option))
                            os.system(f"rm {str(option)}")
                    
                    # rm summary_confidences json files for model 1/2/3/4
                    options = list(af3_job_results_dir.glob(f"*_summary_confidences_*.json"))
                    for option in options:
                        if str(option)[-7:] != "_0.json":
                            print(str(option))
                            os.system(f"rm {str(option)}")

    # check if there are any repeated jobs
    all_submitted_job_names = []
    all_repeated_job_names = [] 
   
    all_af3_jobs_scores = []  # extract AF3 summary scores
    all_af3_jobs_mapped_mutations = [] # extract mapped mutations on predicted interfaces
    for date_worker in os.listdir(os.path.join(root_dir, f"AD_interactome/results/Alphafold3/{project_name}/af3_output/")):
        print(date_worker)
        af3_results_dir = os.path.join(f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/af3_output/", date_worker)
        if os.path.exists(af3_results_dir) and os.path.isdir(af3_results_dir):
            # for each date_worker, generate a results_summary for the mapped mutations
            af3_results_summary_table = [] # store mapped mutations

            # summarize all the job names and check if there are repeated jobs
            job_names = list(os.listdir(af3_results_dir))
            job_names = list(set(job_names) - set(['terms_of_use.md', 'results_summary.txt', 'results_summary_table.tsv']))
            if len(set(all_submitted_job_names) & set(job_names)) > 0: # check if there are any repeated jobs
                print(date_worker, sorted(list(set(all_submitted_job_names) & set(job_names))))
                all_repeated_job_names.extend(sorted(list(set(all_submitted_job_names) & set(job_names))))
            all_submitted_job_names.extend(job_names)

            with open(os.path.join(af3_results_dir, "results_summary.txt"), 'w') as f:
                for af3_job_name in os.listdir(af3_results_dir):
                    if os.path.isdir(os.path.join(af3_results_dir, af3_job_name)):
                        # print(af3_job_name)
                        f.write("========="+af3_job_name+"========="+'\n')

                        # get bait & interactor name
                        bait = af3_job_name.split("_")[0].upper()
                        interactor = af3_job_name.split("_")[1].upper()
                        for isoform in isoforms: # if there are isoforms in the PPI (e.g. O00499-9)
                            if af3_job_name.find(isoform) == 0:
                                bait = isoform.upper().replace("_", "-")
                                interactor = af3_job_name.split("_")[2].upper()
                            elif af3_job_name.find(isoform) > 0:
                                bait = af3_job_name.split("_")[0].upper()
                                interactor = isoform.upper().replace("_", "-")
                        
                        bait_gene = list(PPIs[(PPIs["Bait ID"] == bait) & (PPIs["Interactor ID"] == interactor)]["Bait"])[0]
                        interactor_gene = list(PPIs[(PPIs["Bait ID"] == bait) & (PPIs["Interactor ID"] == interactor)]["Interactors"])[0]

                        # extract sequences of bait and interactor
                        bait_seq = list(PPIs_seq[(PPIs_seq["prot1"] == bait) & (PPIs_seq["prot2"] == interactor)]["seq1"])[0]
                        bait_seq_len = len(bait_seq)
                        interactor_seq = list(PPIs_seq[(PPIs_seq["prot1"] == bait) & (PPIs_seq["prot2"] == interactor)]["seq2"])[0]
                        interactor_seq_len = len(interactor_seq)

                        # load af3 results
                        structure, full_data, summary = _load_af3_files(res_dir=pathlib.Path(os.path.join(af3_results_dir, af3_job_name)), model_idx=0)

                        # check token chain ids
                        token_chain_ids = full_data["token_chain_ids"]
                        f.write(f"Check if chain A is consistent with bait sequence: {token_chain_ids.count('A') == bait_seq_len}\n")
                        f.write(f"Check if chain B is consistent with interactor sequence: {token_chain_ids.count('B') == interactor_seq_len}\n")

                        # extract AF3 scores
                        iptm = summary["iptm"]
                        ptm = summary["ptm"]
                        ranking_score = summary["ranking_score"]
                        f.write(f"Predicted iptm: {iptm}\n")

                        # save AF3 scores
                        all_af3_jobs_scores.append([af3_job_name, bait, interactor, bait_gene, interactor_gene, iptm, ptm, ranking_score])

                        # map mutations
                        # use "Uniprot ID exact"!!!
                        AD_mutations_bait = AD_mutations_all[AD_mutations_all["Uniprot ID exact"] == bait].reset_index(drop=True)
                        if len(AD_mutations_bait) > 0:
                            for i in range(len(AD_mutations_bait)):
                                af3_results_summary_table_row = [af3_job_name, bait, interactor, bait_gene, interactor_gene, iptm]

                                # extract information about the mutation in bait
                                # use "protein position exact"!!!
                                mutation_pos = list(AD_mutations_bait["protein position exact"])[i]
                                mutation_aa_ref = list(AD_mutations_bait["amino acids"])[i].split("/")[0]
                                mutation_aa_alt = list(AD_mutations_bait["amino acids"])[i].split("/")[1]
                                f.write(f"mutation: {mutation_aa_ref}{mutation_pos}{mutation_aa_alt}"+'\n')
                                af3_results_summary_table_row.extend([mutation_pos, mutation_aa_ref, mutation_aa_alt])

                                # check if mutation_aa_ref == aa in the bait sequence
                                f.write(f"Check reference aa of the mutation: {bait_seq[mutation_pos-1] == mutation_aa_ref}\n")

                                # calculatethe max contact_prob between the mutation residue and the opposite chain
                                contact_probs = np.array(full_data["contact_probs"])
                                contact_probs_mut_with_interactor = contact_probs[mutation_pos-1, bait_seq_len:]
                                contact_probs_mut_with_interactor_max = max(contact_probs_mut_with_interactor)
                                max_contact_prob_residue = np.argmax(contact_probs_mut_with_interactor) + 1
                                f.write(f"Max contact prob between mutation position in the bait with the interactor: {contact_probs_mut_with_interactor_max}\n")
                                af3_results_summary_table_row.append(contact_probs_mut_with_interactor_max)
                                f.write(f"Residue on the interactor with max contact prob {max_contact_prob_residue}\n")
                                af3_results_summary_table_row.append(max_contact_prob_residue)

                                # calculate the min distance between the mutation residue and all the alpha carbon atoms from the opposite chain
                                distances = calculate_mutation_residue_ca_distances(structure, "A", "B", mutation_pos, plddt_threshold=25)
                                if distances is not None:
                                    min_distance = sorted(distances.items(), key=lambda item: item[1])[0][1]
                                    min_distance_residue = sorted(distances.items(), key=lambda item: item[1])[0][0]
                                else:
                                    min_distance = None
                                    min_distance_residue = None
                                f.write(f"Min distance between mutation position in the bait with the interactor: {min_distance}\n")
                                af3_results_summary_table_row.append(min_distance)
                                f.write(f"Residue on the interactor with min distance: {min_distance_residue}\n")
                                af3_results_summary_table_row.append(min_distance_residue)

                                af3_results_summary_table.append(af3_results_summary_table_row)
            f.close()
            # save the overall mapped mutations results table (for each date_worker)
            af3_results_summary_table = pd.DataFrame(af3_results_summary_table, columns=["af3_job_name", "bait_id", "interactor_id", "bait_gene", "interactor_gene", "iptm", "mutation_pos", "mutation_aa_ref", "mutation_aa_alt", \
                                                                                        "max_contact_prob", "max_contact_prob_residue", "min_distance", "min_distance_residue"])
            af3_results_summary_table.to_csv(os.path.join(af3_results_dir, "results_summary_table.tsv"), sep='\t', index=None)
            all_af3_jobs_mapped_mutations.append(af3_results_summary_table)

    # save the overall AF3 scores table (for all date_workers)
    all_af3_jobs_scores = pd.DataFrame(all_af3_jobs_scores, columns=["af3_job_name", "bait_id", "interactor_id", "bait_gene", "interactor_gene", "iptm", "ptm", "ranking_score"])
    all_af3_jobs_scores["gene_pair"] = all_af3_jobs_scores.apply(lambda row: "_".join(sorted([row["bait_gene"],row["interactor_gene"]])), axis=1)
    all_af3_jobs_scores["ppi"] = all_af3_jobs_scores.apply(lambda row: "_".join(sorted([row["bait_id"],row["interactor_id"]])), axis=1)
    # remove the duplicated PPIs (keep the one with the highest iptm)
    all_af3_jobs_scores = all_af3_jobs_scores.sort_values(by=["iptm"], ascending=False).drop_duplicates(subset=["ppi"]).reset_index(drop=True)
    
    all_af3_jobs_scores = all_af3_jobs_scores[["af3_job_name", "bait_id", "interactor_id", "bait_gene", "interactor_gene", "gene_pair", "ppi", "iptm", "ptm", "ranking_score"]]
    all_af3_jobs_scores.to_csv(f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/af3_output_summary/all_af3_jobs_scores.tsv", sep='\t', index=None)

    # save the overall mapped mutations results table (for all date_workers)
    all_af3_jobs_mapped_mutations = pd.concat(all_af3_jobs_mapped_mutations, axis=0).reset_index(drop=True)
    all_af3_jobs_mapped_mutations["ppi"] = all_af3_jobs_mapped_mutations.apply(lambda row: "_".join(sorted([row["bait_id"],row["interactor_id"]])), axis=1)
    all_af3_jobs_mapped_mutations["gene_pair"] = all_af3_jobs_mapped_mutations.apply(lambda row: "_".join(sorted([row["bait_gene"],row["interactor_gene"]])), axis=1)
    # # remove the duplicated PPIs (keep the one with the highest iptm)
    # all_af3_jobs_mapped_mutations = all_af3_jobs_mapped_mutations.sort_values(by=["iptm"], ascending=False).drop_duplicates(subset=["ppi"]).reset_index(drop=True)
    # threshold iptm and min_distance to select mapped mutations on predicted interfaces
    iptm_cutoff = 0
    # min_distance_cutoff = 30
    min_distance_cutoff = all_af3_jobs_mapped_mutations["min_distance"].max()
    if min_distance_cutoff == all_af3_jobs_mapped_mutations["min_distance"].max():
        min_distance_cutoff_str = "all"
    else:
        min_distance_cutoff_str = min_distance_cutoff

    # select all pairs with iptm > iptm_cutoff & min_distance <= min_distance_cutoff
    all_af3_jobs_mapped_mutations_cutoff = all_af3_jobs_mapped_mutations[(all_af3_jobs_mapped_mutations["iptm"] > iptm_cutoff) & (all_af3_jobs_mapped_mutations["min_distance"] <= min_distance_cutoff)].reset_index(drop=True)
    # select all pairs with iptm > iptm_cutoff & (min_distance <= min_distance_cutoff | max_contact_prob > 0)
    # all_af3_jobs_mapped_mutations_cutoff = all_af3_jobs_mapped_mutations[(all_af3_jobs_mapped_mutations["iptm"] > iptm_cutoff) & ((all_af3_jobs_mapped_mutations["min_distance"] <= min_distance_cutoff) | (all_af3_jobs_mapped_mutations["max_contact_prob"] > 0))].reset_index(drop=True)
    print(all_af3_jobs_mapped_mutations_cutoff[["af3_job_name", "bait_id", "interactor_id", "bait_gene", "interactor_gene", "iptm", 
                                                "mutation_pos", "mutation_aa_ref", "mutation_aa_alt", 
                                                "max_contact_prob", "min_distance"]].sort_values(by=["bait_gene", "min_distance"]).reset_index(drop=True).head())
    all_af3_jobs_mapped_mutations_cutoff.to_csv(f"/fs/cbsuhy02/storage/tx82/AD_interactome/results/Alphafold3/{project_name}/af3_output_summary/all_af3_jobs_mapped_mutations_iptm{str(iptm_cutoff)}_mindist{min_distance_cutoff_str}.tsv", sep='\t', index=None)