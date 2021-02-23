import os
import pandas as pd
import argparse
import pickle

def mergeResults(families_dir, output_dir, free_energy, second_run):
    families = os.listdir(families_dir)
    dataFrames = []
    for family in families:
        full_family_path = os.path.join(families_dir, family)
        dicGeneIsoform_path = os.path.join(full_family_path, "geneIsoform.pkl")
        with open (dicGeneIsoform_path, 'rb') as handle_gene:
            dicGeneIsoform = pickle.load(handle_gene)
        if len(dicGeneIsoform) == 1:
            continue
        if second_run:
            path_of_csv_output = os.path.join(full_family_path, "AmiRNAs_final_" +free_energy+ ".csv")
            if not os.path.exists(path_of_csv_output):
                path_of_csv_output = os.path.join(full_family_path, "AmiRNAs" + free_energy + ".csv")
        else:
            path_of_csv_output = os.path.join(full_family_path, "AmiRNAs" + free_energy + ".csv")
        if not os.path.exists(path_of_csv_output):
            continue
        df = pd.read_csv(path_of_csv_output, index_col=0)
        dataFrames.append(df)
    final_table = pd.concat(dataFrames)
    if second_run:
        output_path = os.path.join(output_dir, "AmiRNAs_final_" + free_energy + ".csv")
    else:
        output_path = os.path.join(output_dir, "AmiRNAs" +free_energy+ ".csv")
    final_table.to_csv(output_path)
    return








if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_file_path", "-i", help = "path for the WMD3 output file")
    parser.add_argument("--output_dir", "-o", help = "directory for the final output")
    parser.add_argument("--free_energy_thr", "-f", help = "use threshold for energey")
    parser.add_argument("--second_run", "-s", type = int, default=0, help = "run second time")
    args = parser.parse_args()
    input_file_path = args.input_file_path
    output_dir = args.output_dir
    free_energy = args.free_energy_thr
    second_run = args.second_run
    mergeResults(input_file_path, output_dir, free_energy, second_run)