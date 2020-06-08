#!/bin/env/python3

import dataset
import program
import experiment
import argparse
import util
import itertools
import os
import sys
import csv
import rich
import rich.progress


def make_datasets(taxa_count, length, ds_count):
    return [
        dataset.lagrange_dataset(prefix, taxa_count=taxa_count, length=length)
        for prefix in util.base58_generator(ds_count)
    ]


parser = argparse.ArgumentParser()
parser.add_argument("--prefix", type=str, required=True)
parser.add_argument("--regions", type=int, nargs="+", required=True)
parser.add_argument("--taxa", type=int, nargs="+", required=True)
parser.add_argument("--iters", type=int, required=True)
parser.add_argument("--procs", type=int)
parser.add_argument("--program", type=str, required=True)
parser.add_argument("--profile", action='store_true', default=False)
args = parser.parse_args()

if os.path.exists(args.prefix):
    print("Prefix already exists")
    sys.exit(1)

exp_program = [
    program.lagrange(binary_path=os.path.abspath(args.program),
                     profile=args.profile)
]

exp = []

exp_name_format = "{taxa}taxa_{regions}regions"

with rich.progress.Progress() as progress_bar:

    total_datasets = len(args.regions) * len(args.taxa)
    total_work = total_datasets * args.iters
    overall_task = progress_bar.add_task("Overall...",
                                         total=total_datasets + 1)
    make_task = progress_bar.add_task("Making Datasets...",
                                      total=total_datasets)

    for r, t in itertools.product(args.regions, args.taxa):
        exp_path = os.path.join(args.prefix,
                                exp_name_format.format(regions=r, taxa=t))
        exp.append(
            experiment.experiment(exp_path, make_datasets(t, r, args.iters),
                                  exp_program))
        progress_bar.update(make_task, advance=1.0)

    progress_bar.update(overall_task, advance=1.0)

    for e in exp:
        e.run(progress_bar, args.procs)
        progress_bar.update(overall_task, advance=1.0)

    results = []
    for e in exp:
        results.extend(e.collect_results())

    with open(os.path.join(args.prefix, 'results.csv'), 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=results[0].header())
        for result in results:
            writer.writerow(result.write_row())
