#!/bin/env/python3

import argparse
import benchmark

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--regions", type=int, nargs="+", required=True)
    parser.add_argument("--taxa", type=int, nargs="+", required=True)
    parser.add_argument("--iters", type=int, required=True)
    parser.add_argument("--procs", type=int)
    parser.add_argument("--program", type=str, required=True)
    parser.add_argument("--profile", action='store_true', default=False)
    args = parser.parse_args()

    benchmark.run(args.prefix, args.regions, args.taxa, args.iters, args.procs,
                  args.program, args.profile)
