import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def check_regions_bed(regions_bed_path: str) -> Tuple[int, int]:
    """
    Checks if the bed file follows certain conditions.

    The conditions checked are:
        1. The bed file should not have more than 6 columns.
        2. For each row, the second column's value should be greater than the first column's value.
        3. All regions should be of the same size and even.
        4. Column 3 should have unique labels.

    Returns
    -------
    int start of region relative to center.
    int end of region relative to center.
    """

    df = pd.read_csv(regions_bed_path, sep="\t", header=None)
    if df.shape[1] > 6:
        return sys.exit("Your bed file has more than 6 columns. Exiting")
    elif any(df[2] - df[1] < 0):
        return sys.exit(
            "Coordinate in column 3 should be greater than coordinate in column 2 (column indexing starts at 1, not 0). Exiting"
        )
    elif len(set(df[2] - df[1])) > 1:
        return sys.exit("All regions should be of the same size. Exiting")
    elif list(set(df[2] - df[1]))[0] % 2 != 0:
        return sys.exit("All regions should be of even length. Exiting")
    elif len(set(df[3])) != len(df[3]):
        return sys.exit(
            "All regions should have a unique identifier in column 4 (column indexing starts at 1, not 0). Exiting"
        )

    center = int((df[1].iloc[0] + df[2].iloc[0]) / 2)
    start = int(df[1].iloc[0] - center)
    end = int(df[2].iloc[0] - center)

    return start, end


def run_bedtools(args: Tuple[str, str, str, str]) -> str:
    """
    Runs bedtools intersect to get the number of reads in each region.

    Parameters
    ----------
    args : tuple

    Returns
    -------
    str path to bed file with 4 columns: genes, counts, weighted average, stdev
    """
    reads, name, regions, strandedness = args
    # create temporary path for bedtools output
    temp_data_bedtools = Path(Path.cwd(), name + ".bed")

    # run bedtools
    if strandedness:
        cmd = f"bedtools intersect -a {regions} -b {reads} -wa -wb -s > {temp_data_bedtools}"
    else:
        cmd = (
            f"bedtools intersect -a {regions} -b {reads} -wa -wb > {temp_data_bedtools}"
        )
    completed_process = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    if completed_process.returncode != 0:
        error_message = completed_process.stderr.decode("utf-8")
        print(f"{error_message}")
        sys.exit(1)

    return temp_data_bedtools


def run_weighted_avg_stdev(
    path_result_bed: List[str],
    names_corr: Dict[str, str],
    overlap_type: str,
    region_left: str,
    region_right: str,
    output_dir: str,
) -> None:
    """
    Runs readWeightedAverageStdev.cpp to get the weighted average and stdev for each region.

    Parameters
    ----------
    path_result_bed : list
        List of paths to bed files with 4 columns: genes, counts, weighted average, stdev
    names_corr : dict
        Dictionary of names and spike-ins
    overlap_type : str
        Type of data to analyze: 5 or 3 prime ends (for transcriptional data) or centers (for ChIP-seq data). The input should be full fragments, the program will automatically calculate the 5 or 3 prime ends or centers
    region_left : str
        Left region relative to center
    region_right : str
        Right region relative to center
    output_dir : str
        Path to output

    Returns
    -------
    None
    """
    # create json file
    dict_data = {}
    for data_path in path_result_bed:
        base_name = os.path.basename(data_path).split(".")[0]
        dict_data[str(data_path)] = names_corr[base_name]
    json_data = json.dumps(dict_data)

    if os.name == "posix":  # Linux or macOS
        readWeightedAverageStdev = "./readWeightedAverageStdev"
    elif os.name == "nt":  # Windows
        readWeightedAverageStdev = "./readWeightedAverageStdev.exe"

    # run readWeightedAverageStdev
    subprocess.run(
        [
            readWeightedAverageStdev,
            json_data,
            overlap_type,
            output_dir,
            str(region_left),
            str(region_right),
        ]
    )


def parse_args():
    """
    Get arguments
    """

    parser = argparse.ArgumentParser(
        prog="weighted_average_stdev.py",
        description="Calculated the total (corrected) counts, the weighted average and stdev for a user chosen region relative to the TSS",
    )
    parser.add_argument(
        "regions", type=str, help="Bed file of genomic regions of chosen length"
    )
    parser.add_argument(
        "-f",
        dest="fragments",
        metavar="\b",
        type=str,
        required=True,
        nargs="*",
        help="BED file containing reads/fragments. Providing multiple input read files will consider them as replicates. Use -f followed by the filenames, e.g., -f WT1.bed WT2.bed.",
    )
    parser.add_argument(
        "-r",
        dest="range",
        metavar="\b",
        type=int,
        nargs=2,
        required=True,
        help="Range relative to the TSS, for exmaple -r -20 500",
    )
    parser.add_argument(
        "-t",
        dest="overlap_type",
        choices=["5", "3", "centers"],
        required=True,
        help="Type of data to analyze: 5 or 3 prime ends (for transcriptional data) or centers (for ChIP-seq data). The input should be full fragments, the program will automatically calculate the 5 or 3 prime ends or centers",
    )
    parser.add_argument(
        "-c",
        dest="spikein",
        metavar="\b",
        type=float,
        nargs="*",
        required=True,
        help="Spike-in or correction factors",
    )
    parser.add_argument(
        "-s",
        dest="strandedness",
        action="store_true",
        default=False,
        help="If argument is invoked, the program will consider the strandedness of the reads. The default is unstranded",
    )
    parser.add_argument(
        "-o",
        dest="output_dir",
        metavar="\b",
        type=str,
        required=True,
        nargs=1,
        help="Path to output",
    )

    args = parser.parse_args()

    read_file = args.fragments
    region_left, region_right = args.range
    spikeins = args.spikein
    output_dir = args.output_dir[0]
    regions = args.regions

    if len(read_file) != len(spikeins):
        sys.exit("The number of bed files and spike-ins do not match")

    if not os.path.exists(output_dir):
        sys.exit("Output directory does not exist")

    if not os.path.exists(regions):
        sys.exit("Regions file does not exist")

    if int(region_right) < int(region_left):
        sys.exit("Right region must be greater than left region")

    return args


def main(args):
    read_file = args.fragments
    region_left, region_right = args.range
    spikeins = args.spikein
    output_dir = args.output_dir[0]
    overlap_type = args.overlap_type
    regions = args.regions
    strandedness = args.strandedness

    # create temp names
    identifier = [str(uuid.uuid4()) for _ in range(len(spikeins))]

    # dictionary -> name:spikein
    names_corr = dict(zip(identifier, spikeins))

    # check bed file, and get start and end of regions relative to center
    region_start, region_end = check_regions_bed(regions)

    # create iterable for multiprocessing
    data = [
        (
            read,
            name,
            regions,
            strandedness,
        )
        for read, name in zip(read_file, identifier)
    ]

    # run bedtools
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # path to bed file: 4 columns  --- genes, counts, weighted average, stdev
        path_result_bed = list(executor.map(run_bedtools, data))

    # run readWeightedAverageStdev.cpp
    run_weighted_avg_stdev(
        path_result_bed, names_corr, overlap_type, region_left, region_right, output_dir
    )

    # delete temporary files
    for file in path_result_bed:
        os.remove(file)


if __name__ == "__main__":
    args = parse_args()
    main(args)
