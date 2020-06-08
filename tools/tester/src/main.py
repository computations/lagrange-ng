#!/usr/bin/env python3

import argparse
import tester

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, required=True)
    parser.add_argument('--archive', type=str, required=True)
    parser.add_argument('--program', type=str, required=True)
    args = parser.parse_args()

    tester.run(args.prefix, args.archive, args.program)
