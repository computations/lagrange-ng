#!/usr/bin/env python3
import dataset
import program
import experiment
import util
import itertools
import os
import sys
import csv
import rich
import rich.progress
import yaml
import subprocess
import flamegraph
import pandas
import seaborn
import matplotlib


def make_datasets(taxa_count, length, ds_count):
    return [
        dataset.lagrange_dataset(prefix, taxa_count=taxa_count, length=length)
        for prefix in util.base58_generator(ds_count)
    ]


def run(prefix, regions, taxa, iters, procs, program_path, profile,
        flamegraph_cmd):
    if os.path.exists(prefix):
        print("Prefix already exists")
        sys.exit(1)
    os.makedirs(prefix)

    exp_program = [
        program.lagrange(binary_path=os.path.abspath(program_path),
                         profile=profile)
    ]

    exp = []

    exp_name_format = "{taxa}taxa_{regions}regions"

    with rich.progress.Progress() as progress_bar:

        total_datasets = len(regions) * len(taxa)
        total_work = total_datasets * iters
        extra_work = 1
        if profile:
            extra_work += 1
        overall_task = progress_bar.add_task("Overall...",
                                             total=total_datasets + extra_work)
        make_task = progress_bar.add_task("Making datasets...",
                                          total=total_datasets)

        with open(os.path.join(prefix, 'parameters.yaml'), 'w') as yamlfile:
            yamlfile.write(
                yaml.dump(
                    {
                        'prefix': prefix,
                        'regions': regions,
                        'taxa': taxa,
                        'iters': iters,
                        'procs': procs,
                        'program_path': program_path,
                        'profile': profile,
                    },
                    explicit_start=True,
                    explicit_end=True))

        for r, t in itertools.product(regions, taxa):
            exp_path = os.path.join(prefix,
                                    exp_name_format.format(regions=r, taxa=t))
            exp.append(
                experiment.experiment(exp_path, make_datasets(t, r, iters),
                                      exp_program))
            progress_bar.update(make_task, advance=1.0)

        progress_bar.update(overall_task, advance=1.0)

        for e in exp:
            e.run(progress_bar, procs)
            progress_bar.update(overall_task, advance=1.0)

        if not profile:
            results = []
            for e in exp:
                results.extend(e.collect_results())

            with open(os.path.join(prefix, 'results.csv'), 'w') as csv_file:
                writer = csv.DictWriter(csv_file,
                                        fieldnames=results[0].header())
                writer.writeheader()
                for result in results:
                    writer.writerow(result.write_row())

            dataframe = pandas.read_csv(os.path.join(prefix, 'results.csv'))
            seaborn.set_style("whitegrid")
            plot = seaborn.FacetGrid(dataframe,
                                     row="taxa",
                                     col="regions",
                                     height=7,
                                     sharex=False,
                                     margin_titles=True).map(seaborn.distplot,
                                                             "time",
                                                             hist=False,
                                                             rug=True)
            plot.savefig(os.path.join(prefix, 'times.png'))

        else:
            fg_work = len(exp) * len(exp[0].datasets)
            fg_task = progress_bar.add_task("Making Flamegraphs...",
                                            total=fg_work)

            for e in exp:
                for d in e.datasets:
                    flamegraph.build(d)
                    progress_bar.update(fg_task, advance=1.0)

            progress_bar.update(overall_task, advance=1.0)
