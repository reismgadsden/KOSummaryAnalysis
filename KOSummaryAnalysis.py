"""
KOSummaryAnalysis.py

This program will analyze the annotated .tsv/.csv as well as the
Trinity output file in order to pull the following statistics for
each KO
- Total number of transcripts
- Total number of unique gene ids
- Average of log fold change
- Minimum of log fold change
- Maximum of log fold change
- P-value of the highest log fold change

author - Reis Mercer Gadsden
version - 2023/02/01
institution - Appalachian State University

github -
"""

"""""""""""
" IMPORTS "
"""""""""""
import pandas
import numpy
import os
import json


class KOSummaryAnalysis:

    _init_df = ""
    _ko_dict = dict()

    def __init__(self, input_file, output_file=None):
        if not os.path.isfile(input_file):
            print("Please provide a valid input file!")
            exit(-1)

        if output_file is None:
            output_file = input_file.replace(".csv", "_KO_Summary_Analaysis_Output.csv")

        self._init_df = pandas.read_csv(input_file)
        target_cols = ["TrinityID", "KO", "log2FoldChange", "pvalue"]

        self._init_df = self._init_df[target_cols]
        self.construct_ko_dict()

    def construct_ko_dict(self):

        for index, row in self._init_df.iterrows():
            # print("Reading row: " + str(index))

            KO = row["KO"]
            TrinityID = row["TrinityID"]
            l2fc = row["log2FoldChange"]
            pv = row["pvalue"]

            if KO in self._ko_dict:
                self._ko_dict[KO]["Total"] += 1
                self._ko_dict[KO]["IDs"].append(TrinityID)
                self._ko_dict[KO]["Stats"].append((TrinityID, l2fc, pv))
            else:
                self._ko_dict[KO] = dict()
                self._ko_dict[KO]["Total"] = 1
                self._ko_dict[KO]["IDs"] = [TrinityID]
                self._ko_dict[KO]["Stats"] = [(TrinityID, l2fc, pv)]

            # optionally dump intermediate dictionary to a json
            # with open("stats.json", "w+") as file:
            #     json.dump(self._ko_dict, file, indent=2)
            #     file.close()


if __name__ == "__main__":
    ko = KOSummaryAnalysis("phegg_unique_stat.csv")