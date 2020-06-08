#!/usr/bin/env python3

import os
import numpy
import ete3
import util
import base58


class dataset:
    _file_prefix = "generic_prog"

    def __init__(self, path, **kwargs):
        self._path = path

    def add_prefix_dir(self, prefix):
        self._path = os.path.join(prefix, self._path)

    def _lock_paths(self):
        self._path = os.path.abspath(self._path)

    @property
    def path(self):
        return self._path

    def write(self):
        raise NotImplementedError("Please implement the write function")

    def _make_path(self):
        self._lock_paths()
        os.makedirs(self.path, exist_ok=True)

    @property
    def profile_path(self):
        return os.path.join(self.path, 'perf.data')

    @property
    def profile_text_filename(self):
        return '{}.perf'.format(self._file_prefix)

    @property
    def profile_text_path(self):
        return os.path.join(self.path, self.profile_text_filename)

    @property
    def folded_profile_filename(self):
        return '{}.perf.folded'.format(self._file_prefix)

    @property
    def folded_profile_path(self):
        return os.path.join(self.path, self.folded_profile_filename)

    @property
    def flamegraph_filename(self):
        return '{}_flamegraph.svg'.format(self._file_prefix)

    @property
    def flamegraph_path(self):
        return os.path.join(self.path, self.flamegraph_filename)


class lagrange_dataset(dataset):
    _file_prefix = "lagrange_exp"
    _lagrange_config =\
    """
treefile = {treefile}
datafile = {datafile}
areanames = {areanames}
ancstate = _all_
states
"""

    def __init__(self, path, **kwargs):
        """Create a lagrange dataset

        Keyword arguments:
        length      -- Number of regions to generate
        taxa_count  -- Number of taxa to generate
        """
        super().__init__(path, **kwargs)
        self._file_prefix = self._file_prefix + "_" + util.make_random_nonce()
        self._tree = ete3.Tree()
        self._length = kwargs['length']
        self._tree.populate(
            kwargs['taxa_count'],
            names_library=[
                s for s in util.base26_generator(kwargs['taxa_count'])
            ])
        self._alignment = lagrange_dataset.generate_alignment(
            self.taxa_set, kwargs['length'])
        self._area_names = [
            "R" + str.upper(s) for s in util.base26_generator(self.length)
        ]

    def make_lagrange_file(self):
        return self._lagrange_config.format(
            treefile=self.tree_filename,
            datafile=self.alignment_filename,
            areanames=self.area_names_string,
        )

    def write(self):
        self._make_path()
        self._write_lagrange_conf()
        self._write_treefile()
        self._write_alignmentfile()

    def _write_alignmentfile(self):
        with open(self.alignment_path, 'w') as outfile:
            outfile.write("{taxa_count} {length}\n".format(taxa_count=len(
                self.taxa_set),
                                                           length=self.length))
            for taxa_name, seq in self._alignment.items():
                outfile.write("{taxa_name} {seq}\n".format(
                    taxa_name=taxa_name, seq=''.join([str(s) for s in seq])))

    def _write_lagrange_conf(self):
        with open(self.lagrange_config_path, 'w') as outfile:
            outfile.write(self.make_lagrange_file())

    def _write_treefile(self):
        with open(self.tree_path, 'w') as outfile:
            outfile.write(self._tree.write(format=5))

    @property
    def lagrange_config_filename(self):
        return self._file_prefix + ".conf"

    @property
    def lagrange_config_path(self):
        return os.path.join(self.path, self.lagrange_config_filename)

    @property
    def area_names_string(self):
        return ' '.join(self._area_names)

    @property
    def length(self):
        return self._length

    @property
    def region_count(self):
        return self.length

    @property
    def taxa_count(self):
        return len(self.taxa_set)

    @property
    def taxa_set(self):
        return [n.name for n in self._tree.get_leaves()]

    @property
    def alignment(self):
        return self._alignment

    @property
    def tree_filename(self):
        return self._file_prefix + ".nwk"

    @property
    def tree_path(self):
        return os.path.join(self.path, self.tree_filename)

    @property
    def alignment_filename(self):
        return self._file_prefix + ".phy"

    @property
    def alignment_path(self):
        return os.path.join(self.path, self.alignment_filename)

    @staticmethod
    def generate_alignment(taxa_set, length):
        return {t: numpy.random.choice(2, length) for t in taxa_set}
