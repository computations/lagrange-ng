#!/bin/env/python3

import argparse
import benchmark
import git
import os
import datetime
import util
import rich

def git_describe(repo):
    description = repo.head.commit.hexsha[0:7]
    for t in repo.tags:
        if t.commit == repo.head.commit:
            description = t.name
            break
    return description

if __name__ == "__main__":

    SOURCE_DIR = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
    DEFAULT_PROGRAM = os.path.abspath(
        os.path.join(SOURCE_DIR, "../../../bin/lagrange"))
    DEFAULT_PROGRAM_PROFILE = os.path.abspath(
        os.path.join(SOURCE_DIR, "../../../bin/lagrange_profile"))

    flamegraph_cmd = [
        "bash",
        os.path.abspath(os.path.join(SOURCE_DIR, '../../build_flamegraphs.sh'))
    ]

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str)
    parser.add_argument("--regions", type=int, nargs="+", default=[4, 5])
    parser.add_argument("--taxa", type=int, nargs="+", default=[10, 50])
    parser.add_argument("--iters", type=int, default=100)
    parser.add_argument("--procs", type=int)
    parser.add_argument("--program", type=str)
    parser.add_argument("--profile", action='store_true', default=False)
    args = parser.parse_args()

    if args.profile and args.program is None:
        args.program = DEFAULT_PROGRAM_PROFILE
    elif args.program is None:
        args.program = DEFAULT_PROGRAM

    if args.prefix is None:
        repo = git.Repo(os.path.abspath(os.path.join(SOURCE_DIR, '../../../')))
        commit_string = datetime.datetime.now().strftime('%Y-%m-%d') + "_"\
                + git_describe(repo) + "_"\
                + util.make_random_nonce()

        if args.profile:
            args.prefix = os.path.abspath(
                os.path.join(SOURCE_DIR, '../profiles', commit_string))
        else:
            args.prefix = os.path.abspath(
                os.path.join(SOURCE_DIR, '../timings', commit_string))
        rich.print("Placing results in [red bold]{}[/red bold]".
                format(os.path.relpath(args.prefix)))

    benchmark.run(args.prefix, args.regions, args.taxa, args.iters, args.procs,
                  args.program, args.profile, flamegraph_cmd)
