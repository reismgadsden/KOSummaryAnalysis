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
import re

"""""""""""
" IMPORTS "
"""""""""""
import pandas
import numpy
import os
import json


class KOSummaryAnalysis:

    _init_df = ""
    _output_df = ""
    _ko_dict = dict()
    _output_file = ""

    def __init__(self, input_file, output_file=None):
        if not os.path.isfile(input_file):
            print("Please provide a valid input file!")
            exit(-1)

        if output_file is None:
            output_file = input_file.replace(".csv", "_KO_Summary_Analaysis_Output.csv")
        self._output_file = output_file

        self._init_df = pandas.read_csv(input_file)
        target_cols = ["TrinityID", "KO", "log2FoldChange", "pvalue"]

        self._init_df = self._init_df[target_cols]
        self.construct_ko_dict()

        self.build_stats_df()

        self._output_df.to_csv(self._output_file)

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

    def build_stats_df(self):

        kos = []

        trans_count = []
        unique_count = []

        l2fc_avg = []

        l2fc_max = []
        pv_max = []
        trans_max = []

        l2fc_min = []
        pv_min = []
        trans_min = []

        for ko in self._ko_dict:

            kos.append(ko)
            trans_count.append(self._ko_dict[ko]["Total"])

            unique_ids = set()
            for trans in self._ko_dict[ko]["IDs"]:
                t_id = re.search(r"DN[0-9]*", trans)
                unique_ids.add(t_id.group(0))
            unique_count.append(len(unique_ids))

            l2fc_cum = []
            ko_min = ("xxx", 100, 0)
            ko_max = ("xxx", -100, 0)

            for tup in self._ko_dict[ko]["Stats"]:
                l2fc_cum.append(tup[1])

                if tup[1] < ko_min[1]:
                    ko_min = tup

                if tup[1] > ko_max[1]:
                    ko_max = tup

            l2fc_avg.append(sum(l2fc_cum)/len(l2fc_cum))

            l2fc_max.append(ko_max[1])
            pv_max.append(ko_max[2])
            trans_max.append(ko_max[0])

            l2fc_min.append(ko_min[1])
            pv_min.append(ko_min[2])
            trans_min.append(ko_min[0])


        output_dict = {
            "KO": kos,
            "Total Unique Transcript": trans_count,
            "Total Unique Gene IDs": unique_count,
            "l2fc Average": l2fc_avg,
            "l2fc Max": l2fc_max,
            "P-value Max": pv_max,
            "Max Transcript ID": trans_max,
            "l2fc Min": l2fc_min,
            "P-value Min": pv_min,
            "Min Transcript ID": trans_min
        }

        self._output_df = pandas.DataFrame(data=output_dict)


if __name__ == "__main__":
    ko = KOSummaryAnalysis("phegg_unique_stat.csv")