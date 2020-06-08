#!/usr/bin/env python3

import subprocess
import util
import filecmp


class lagrange:
    def __init__(self, lagrange_path):
        self._lagrange_path = lagrange_path

    def run(self, path, config_file):
        with util.directory_guard(path):
            with open('lagrange.log', 'w') as logfile:
                subprocess.run([self._lagrange_path, config_file],
                               stdout=logfile,
                               stderr=logfile)


class lagrange_results:
    def __init__(self, bgkey, bgstates):
        self._bgkey = bgkey
        self._bgstates = bgstates

    @property
    def bgkey(self):
        return self._bgkey

    @property
    def bgstates(self):
        return self._bgstates

    def __eq__(self, other):
        return filecmp.cmp(self._bgkey, other._bgkey, shallow=False) and\
                filecmp.cmp(self._bgstates, other._bgstates, shallow=False)
