import pandas as pd
import re
import os
import argparse

# path_paper_results = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\final_amiRNA_L10.csv"
# path_my_results = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\final_results_8_10_0.75.csv"
# res_path = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\filtered_final_results_8_10_0.75.csv"
def removeAmiRNAsFromPaperResults(path_my_results, path_paper_results, res_path, res_path_filtered):
    Bsal = 'GGTCTC'
    my_df = pd.read_csv(path_my_results, index_col=0)
    paper_df = pd.read_csv(path_paper_results, index_col=0)
    amiRNA_in_paper = set(list(paper_df.index))
    my_amiRNAs = my_df.index
    my_df['AmiRNAs'] = my_amiRNAs
    my_amiRNAs_set = set(my_amiRNAs)
    amiRNA_to_remain = my_amiRNAs_set.difference(amiRNA_in_paper)
    amiRNAs_to_filter = list(amiRNA_in_paper)
    amiRNAs_effectively_to_filter = list(my_amiRNAs_set.intersection(amiRNA_in_paper))

    df_final = my_df[my_df.AmiRNAs.isin(list(amiRNA_to_remain))]
    filtered_df = my_df[my_df.AmiRNAs.isin(amiRNAs_to_filter)]
    #df_final.index = [i for i in range(1, len(amiRNA_to_remain)+ 1)]

    columns = list(df_final.columns)
    cols = [columns[-1]] +  [columns[0]] + columns[1:-1]
    df_final = df_final[cols]
    filtered_df = filtered_df[cols]
    lst_reason_for_filtering = ["in other list" for i in range(len(amiRNAs_effectively_to_filter))]

    list_of_miRs_to_remain = []
    list_to_filter = []
    all_AmiRNAs = df_final['AmiRNAs']
    for AmiRNA in all_AmiRNAs:
        if not (Bsal in AmiRNA):
            list_of_miRs_to_remain.append(AmiRNA)
        else:
            list_to_filter.append(AmiRNA)
    df_filtered_Bsal = df_final[df_final.AmiRNAs.isin(list_to_filter)]
    lst_reason_for_filtering += ["Bsal" for i in range(len(list_to_filter))]
    dataFrames_filtered = [filtered_df, df_filtered_Bsal]
    final_filtered_table = pd.concat(dataFrames_filtered)
    final_filtered_table["Reason"] = lst_reason_for_filtering
    df_final = df_final[df_final.AmiRNAs.isin(list_of_miRs_to_remain)]
    df_final.index = [i for i in range(1, len(list(df_final['AmiRNAs'])) + 1)]
    final_filtered_table.index = [i for i in range(1, len(list(final_filtered_table['AmiRNAs'])) + 1)]

    df_final.to_csv(res_path)
    final_filtered_table.to_csv(res_path_filtered)

def removeBsalSeqFromAntisense(path_remained, path_filtered):
    Bsal = 'GGTCTC'
    df_final = pd.read_csv(path_remained, index_col=0)
    columns = list(df_final.columns)
    df_filter_withoutAtisense_col = df_final[columns[:-1]]
    df_filter_withoutAtisense_col["Reason"] = ["Bsal in Antisense" for i in range(len(list(df_final["AmiRNAs"])))]
    df_filtered = pd.read_csv(path_filtered, index_col=0)
    numberOfMiRS = len(df_final['AmiRNAs'])
    to_filter = []
    to_remain = []
    for i in range(numberOfMiRS):
        antiSense = df_final.iloc[i]['miR*s']
        amiRNA = df_final.iloc[i]['AmiRNAs']
        if Bsal in antiSense:
            to_filter.append(amiRNA)
        else:
            to_remain.append(amiRNA)
    df_Bsal_antisense_filtered = df_filter_withoutAtisense_col[df_filter_withoutAtisense_col.AmiRNAs.isin(to_filter)]
    dfs_to_concat = [df_filtered, df_Bsal_antisense_filtered]
    final_filtered_table = pd.concat(dfs_to_concat)
    final_filtered_table.index = [i for i in range(1, len(list(final_filtered_table['AmiRNAs'])) + 1)]
    final_filtered_table.to_csv(path_filtered)
    df_final = df_final[df_final.AmiRNAs.isin(list(to_remain))]
    df_final.index = [i for i in range(1, len(list(df_final['AmiRNAs'])) + 1)]
    df_final.to_csv(path_remained)


if __name__ == "__main__":
    path_paper_results = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\final_amiRNA_L10.csv"
    path_my_results = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\final_results_8_10_0.75.csv"
    res_path = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\final_design\\remained_miRs_results_8_10_0.75.csv"
    res_path_filtered = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\final_design\\deleted_miRs_results_8_10_0.75.csv"
    #removeAmiRNAsFromPaperResults(path_my_results, path_paper_results, res_path, res_path_filtered)
    removeBsalSeqFromAntisense(res_path, res_path_filtered)







