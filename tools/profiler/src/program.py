#!/usr/bin/env python3

import subprocess
import os
import util
import result


class program:
    def __init__(self, **kwargs):
        self._binary_path = kwargs['binary_path']
        self._profile = kwargs['profile']

    def run(self, *args, **kwargs):
        raise NotImplementedError("Run is not implemented in the base class")

    @property
    def binary(self):
        return os.path.abspath(self._binary_path)

    @property
    def profile_cmd(self):
        if self._profile is True:
            return "perf record -g".split()
        return []


class lagrange(program):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, dataset):
        with util.directory_guard(dataset.path):
            with open('lagrange.log', 'w') as logfile:
                cmd = []
                cmd.extend(self.profile_cmd)
                cmd.extend([self.binary, dataset.lagrange_config_path])
                subprocess.run(cmd, stdout=logfile, stderr=logfile)

    def get_result(self, dataset):
        return lagrange_result(dataset)


class lagrange_result(result.result):
    _logfile_filename = "lagrange.log"
    _program = "lagrange"

    def __init__(self, dataset, **kwargs):
        super().__init__(**kwargs)
        self._dataset = dataset
        with open(self.logfile_path) as logfile:
            time_line = list(logfile)[-1]

        self._time = lagrange_result._parse_timeline(time_line)

    @property
    def logfile_path(self):
        return os.path.join(self._dataset.path, self._logfile_filename)

    @staticmethod
    def _parse_timeline(line):
        line = line.strip()
        if "Analysis took: " not in line:
            print(line)
            raise RuntimeError("The time line of the log file was malformed")
        prefix_length = len("Analysis took: ")
        return float(line[prefix_length:-2])

    def write_row(self):
        return {
            'program': self.program,
            'taxa': self._dataset.taxa_count,
            'regions': self._dataset.region_count,
            'time': self._time
        }

    def header(self):
        return ['program', 'taxa', 'regions', 'time']
