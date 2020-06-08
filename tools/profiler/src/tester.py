#!/usr/bin/env python3

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=string, required=True)
parser.add_argument("--datapath", type=string, required=True)
args = parser.parse_args()
