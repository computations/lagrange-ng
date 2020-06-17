#!/usr/bin/env python3
import os
import tarfile
import argparse
import pathlib
import util
import lagrange
import rich
import rich.console
import rich.progress
import yaml


def file_type(filename):
    basename, extension = os.path.splitext(filename)
    if extension == '.log':
        return "log"
    if extension == '.tre':
        subbasename, subextension = os.path.splitext(basename)
        if subextension == '.bgkey':
            return 'bgkey'
        if subextension == '.bgstates':
            return 'bgstates'
    if extension == '.nwk':
        return 'newick'
    if extension == '.phy':
        return 'phylip'
    if extension == '.conf':
        return 'config'
    return 'unknown'


def replace_prefix_list(path, prefix):
    tmp_path = pathlib.Path(prefix[0])
    for i in range(1, len(path.parts)):
        if i < len(prefix):
            tmp_path /= prefix[i]
        else:
            tmp_path /= path.parts[i]
    return tmp_path


def replace_prefix_string(path, prefix):
    tmp_path = pathlib.Path(prefix)
    for i in range(1, len(path.parts)):
        tmp_path /= path.parts[i]
    return tmp_path


def replace_prefix(path, prefix):
    if type(prefix) is list:
        return replace_prefix_list(path, prefix)
    if type(prefix) is str:
        return replace_prefix_string(path, prefix)


def replace_basename(path, basename):
    tmp_path = pathlib.Path(path.parts[0])
    for part in path.parts[1:-1]:
        tmp_path /= part
    tmp_path /= basename
    return tmp_path


def extract_name_without_nonce(basename):
    prefix_with_nonce = basename.split('.')[0]
    prefix = '_'.join(prefix_with_nonce.split('_')[0:-1])
    for ext in basename.split('.')[1:]:
        prefix += '.'
        prefix += ext
    return prefix


def extract(tar, member, member_path, prefix):
    os.makedirs(member_path.parent, exist_ok=True)
    with member_path.open('wb') as outfile:
        outfile.write(tar.extractfile(member).read())


def extract_with_path(tar, member, prefix):
    member_path = pathlib.Path(member.name)
    member_path = replace_prefix(member_path, prefix)
    extract(tar, member, member_path, prefix)


def extract_with_clean_name(tar, member, prefix):
    member_path = pathlib.Path(member.name)
    member_path = replace_prefix(member_path, prefix)
    basename = member_path.name
    basename = extract_name_without_nonce(basename)
    member_path = replace_basename(member_path, basename)
    extract(tar, member, member_path, prefix)


def extract_with_set_name(tar, member, name, prefix):
    member_path = pathlib.Path(member.name)
    member_path = replace_prefix(member_path, prefix)
    basename = '.'.join([name] + member_path.name.split('.')[1:])
    member_path = replace_basename(member_path, basename)
    extract(tar, member, member_path, prefix)


def find_config_file(files):
    for f in files:
        if file_type(f) == 'config':
            return f
    return None


def select_files_with_prefix(prefix, files):
    ret = []
    for f in files:
        if f[:len(prefix)] == prefix:
            ret.append(f)
    return ret


def get_lagrange_ext(filename):
    return os.path.splitext(os.path.splitext(filename)[0])[1]


def order_bgkey_bgstates(files):
    if len(files) != 2:
        raise Exception(
            "This function requires a list of two arguements, got {}".format(
                len(files)))
    if get_lagrange_ext(files[0]) == '.bgkey':
        return (files[0], files[1])
    return (files[1], files[0])


def select_bgkey_bgstates(files):
    ret = []
    for f in files:
        if get_lagrange_ext(f) == '.bgkey' or get_lagrange_ext(
                f) == '.bgstates':
            ret.append(f)
    return ret


def select_expected(path, files):
    bgkey, bgstates = order_bgkey_bgstates(
        select_bgkey_bgstates(select_files_with_prefix('expected', files)))
    return (os.path.join(path, bgkey), os.path.join(path, bgstates))


def select_experiment(path, files):
    try:
        bgkey, bgstates = order_bgkey_bgstates(
            select_bgkey_bgstates(
                select_files_with_prefix('lagrange_exp', files)))
    except:
        raise Exception("Failed to find the experiment files in: '{}'".format(
            ','.join(files)))
    return (os.path.join(path, bgkey), os.path.join(path, bgstates))


def compare_results_expected(path):
    files = [
        f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
    ]
    expected_results = lagrange.lagrange_results(
        *select_experiment(path, files))
    experiment_results = lagrange.lagrange_results(
        *select_expected(path, files))
    return experiment_results == expected_results


def get_taxa_count(trial):
    taxa, _ = trial.split('_')
    taxa_count = taxa.strip('taxa')
    return int(taxa_count)


def get_region_count(trial):
    _, region = trial.split('_')
    region_count = region.strip('regions')
    return int(region_count)


def print_current_trial(console, path):
    trial, iteration = os.path.split(path)
    trial = os.path.relpath(trial)

    console.print("Running trial:", get_taxa_count(trial), "taxa",
                  get_region_count(trial), "regions", "iteration", iteration)


def run(prefix, archive, program):
    runner = lagrange.lagrange(program)
    failed_paths = []

    console = rich.console.Console()
    with rich.progress.Progress() as progress:
        tar = tarfile.open(archive)

        member_count = len(tar.getmembers())
        extract_task = progress.add_task("[red]Extracting...",
                                         total=member_count)

        for member in tar.getmembers():
            member_ft = file_type(member.name)
            if member_ft == 'newick' or member_ft == 'phylip' or member_ft ==\
                'config':
                extract_with_path(tar, member, prefix)
            if member_ft == 'log' or member_ft == 'bgkey' or member_ft ==\
                'bgstates':
                extract_with_set_name(tar, member, 'expected', prefix)
            progress.update(extract_task, advance=1.0)

        with util.directory_guard(prefix):
            work_paths = []
            #work_task = progress.add_task("[blue]Building tests...")
            for root, dirs, files in os.walk('.'):
                config_file = find_config_file(files)
                if not config_file is None:
                    work_paths.append((root, config_file))
                #progress.update(work_task)

            test_task = progress.add_task("[green]Testing...",
                                          total=len(work_paths))

            for path, config_file in work_paths:
                runner.run(path, config_file)
                progress.update(test_task, advance=1.0)
                if not compare_results_expected(path):
                    failed_paths.append(path)
    with open(os.path.join(prefix, "failed_paths.yaml"), "w") as outfile:
        outfile.write(yaml.dump(failed_paths))
    if len(failed_paths) != 0:
        console.print("failed paths:", sorted(failed_paths))
        console.print("Total of {} paths failed".format(len(failed_paths)))
    else:
        console.print("[bold green]All Clear!")
