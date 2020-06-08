#!/usr/bin/env python3

import dataset
import util
import subprocess
import os

perf_script = 'perf script'
fold_cmd = 'stackcollapse-perf'
flamegraph_cmd = 'flamegraph'


def build(dataset, keep=False):
    with util.directory_guard(dataset.path):
        with open(dataset.profile_text_path, 'w') as outfile:
            subprocess.run(perf_script.split(' '), stdout=outfile)

        with open(dataset.folded_profile_path, 'w') as outfile:
            subprocess.run([fold_cmd, dataset.profile_text_path],
                           stdout=outfile)

        with open(dataset.flamegraph_path, 'w') as outfile:
            subprocess.run([flamegraph_cmd, dataset.folded_profile_path],
                           stdout=outfile)

        if not keep:
            os.remove(dataset.profile_text_path)
            os.remove(dataset.folded_profile_path)
