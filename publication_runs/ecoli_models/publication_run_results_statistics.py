#!/usr/bin/env python3
#
# Copyright 2020-2021 PSB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
This script generates the tables (as console text) of the Results section of
ASTHERISC's publication.
This script is meant to be started from the ASTHERISC package's main folder.
"""
from os import listdir
from os.path import isfile, join
from typing import Dict, List
from copy import deepcopy


def sort_filenames(run_result_filename):
    if "ecolicore2double" in run_result_filename:
        value = 10
    elif "ecolicore2triple" in run_result_filename:
        value = 1000000
    elif "iML1515double" in run_result_filename:
        value = 100000000

    if "Inf" in run_result_filename:
        value *= 0.05

    value *= float(run_result_filename.split("maxYield_")
                   [1].split("__numMaxExchanges")[0])
    return value


base_path = "./publication_runs/ecoli_models/run_results/"
run_result_filenames = [x for x in listdir(base_path) if (
    isfile(join(base_path, x)) and "report__" in x)]

original_mdfs_per_file: Dict[str, Dict[float, int]] = {}
community_mdfs_per_file: Dict[str, Dict[float, int]] = {}

sorted_run_result_filenames = sorted(run_result_filenames, key=sort_filenames)

base_list = []
line_list = []
for _ in range(4):
    line_list.append(float("nan"))
for _ in range(6):
    base_list.append(deepcopy(line_list))

mdf_advantage_list = deepcopy(base_list)
exchanges_list = deepcopy(base_list)
advantage_pct_list = deepcopy(base_list)
for run_result_filename in sorted_run_result_filenames:
    original_mdfs: Dict[float, int] = {}
    community_mdfs: Dict[float, int] = {}

    with open(base_path+run_result_filename, "r") as f:
        run_result_lines = f.readlines()
    run_result_lines = [x.replace("\n", "") for x in run_result_lines]

    in_result = False
    mdf_advantages: List[float] = []
    i = 0
    for run_result_line in run_result_lines:
        if run_result_line.startswith("==TARGET WITH COMMUNITY ADVANTAGE:"):
            in_result = True

        if in_result:
            if run_result_line.startswith("Optimal MDF without community with minimal necessary yield as minimal yield constraint:"):
                mdf_without_community = float(run_result_line.replace(
                    "Optimal MDF without community with minimal necessary yield as minimal yield constraint:", ""))

            if run_result_line.startswith("Optimal MDF with community with minimal necessary yield as minimal yield constraint (if no warning given):"):
                mdf_with_community = float(run_result_line.replace(
                    "Optimal MDF with community with minimal necessary yield as minimal yield constraint (if no warning given):", ""))

                mdf_advantages.append(mdf_with_community-mdf_without_community)
                in_result = False

        if run_result_line.startswith("Maximal extra exchanges:"):
            maximal_exchanges = run_result_lines[i+1]
        if run_result_line.startswith("Minimal extra exchanges:"):
            minimal_exchanges = run_result_lines[i+1]
        if run_result_line.startswith("Mean extra exchanges:"):
            mean_exchanges = run_result_line.split(": ")[1]  # TODO: run_result_lines[i+1]
        if run_result_line.startswith("Number metabolites with found community MDF benefit:"):
            number_with_mdf_advantage = int(run_result_lines[i+1])
        if run_result_line.startswith("Number metabolites without found community MDF benefit:"):
            number_without_mdf_advantage = int(run_result_lines[i+1])
        if run_result_line.startswith("Number (fraction) producible target metabolites without dG0=NaN reactions: "):
            number_producible_metabolites = int(run_result_line.split(": ")[1].split("(")[0].rstrip())  # TODO: run_result_lines[i+1].split("(")[0]
        i += 1

    mdf_advantage_str = str(round(min(mdf_advantages), 2))+"/"+str(round(sum(
        mdf_advantages)/len(mdf_advantages), 2))+"/"+str(round(max(mdf_advantages), 2))
    print("File: ", run_result_filename)
    print("MDF advantage (min/mean/max): ", mdf_advantage_str)
    exchanges_str = str(minimal_exchanges)+"/" + \
        str(round(float(mean_exchanges), 2))+"/"+str(maximal_exchanges)
    print("Exchanges (min/mean/max): ", exchanges_str)
    ### pct_str = str(round(100*(number_with_mdf_advantage/(number_with_mdf_advantage +
    ###                                                     number_without_mdf_advantage)), 2))+"% ("+str(number_with_mdf_advantage)+")"
    pct_str = str(round(100*(number_with_mdf_advantage/number_producible_metabolites), 2))+"% ("+str(number_with_mdf_advantage)+")"
    print("Percentage of metabolites (and absolute number) with community MDF advantage: ", pct_str)
    print("")

    if "ecolicore2double" in run_result_filename:
        list_line = 1
    elif "ecolicore2triple" in run_result_filename:
        list_line = 3
    elif "iML1515double" in run_result_filename:
        list_line = 5

    if "Inf" in run_result_filename:
        list_line += 1

    list_line -= 1

    minyield = float(run_result_filename.split("maxYield_")
                     [1].split("__numMaxExchanges")[0])
    if minyield == 0.4:
        list_row = 1
    elif minyield == 0.6:
        list_row = 2
    elif minyield == 0.8:
        list_row = 3
    elif minyield == 0.98:
        list_row = 4

    list_row -= 1

    mdf_advantage_list[list_line][list_row] = mdf_advantage_str
    exchanges_list[list_line][list_row] = exchanges_str
    advantage_pct_list[list_line][list_row] = pct_str

    for run_result_line in run_result_lines:
        part = "Optimal MDF without community with minimal necessary yield as minimal yield constraint:"
        if run_result_line.startswith(part):
            mdf = round(float(run_result_line.replace(part, "")), 2)
            if mdf in original_mdfs.keys():
                original_mdfs[mdf] += 1
            else:
                original_mdfs[mdf] = 1

        part = "Optimal MDF with community with minimal necessary yield as minimal yield constraint (if no warning given):"
        if run_result_line.startswith(part):
            mdf = round(float(run_result_line.replace(part, "")), 2)
            if mdf in community_mdfs.keys():
                community_mdfs[mdf] += 1
            else:
                community_mdfs[mdf] = 1
    original_mdfs_per_file[run_result_filename] = original_mdfs
    community_mdfs_per_file[run_result_filename] = community_mdfs


def print_list(x):
    for line in x:
        print(";".join(line))


def print_mdfs_per_file(community_mdfs_per_file):
    filenames = [
        "report__model_ecolicore2double__maxYield_0.8__numMaxExchanges_Inf.txt",
        "report__model_ecolicore2triple__maxYield_0.8__numMaxExchanges_Inf.txt",
        "report__model_iML1515double__maxYield_0.8__numMaxExchanges_Inf.txt",
    ]
    for filename in filenames:
        print(filename)
        entry = community_mdfs_per_file[filename]
        while True:
            max_value = max(entry.values())
            keys_with_max_value = []
            for key in entry.keys():
                if entry[key] == max_value:
                    keys_with_max_value.append(key)
            keys_with_max_value = sorted(keys_with_max_value)[::-1]
            for key in keys_with_max_value:
                print(str(key)+": "+str(entry[key]))
                del(entry[key])
            if len(entry.keys()) == 0:
                break


print("==FINAL STATISTICS==")
print("MDF advantage statistics:")
print_list(mdf_advantage_list)
print("Exchange number statistics:")
print_list(exchanges_list)
print("Statistics about number of metabolites with MDF advantage:")
print_list(advantage_pct_list)
print("")
print_mdfs_per_file(community_mdfs_per_file)
print("")
