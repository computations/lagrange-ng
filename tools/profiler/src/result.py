#!/usr/bin/env python3

class result:
    _program = "None"
    def __init__(self, **kwargs):
        pass

    def write_row(self):
        raise NotImplementedError

    def header(self):
        raise NotImplementedError

    @property
    def program(self):
        return self._program
