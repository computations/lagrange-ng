#!/usr/bin/env python3

import math
import string
import os
import base58

CAP_ALPHABET = 'ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz'
BITCOIN_ALPHABET = '123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz'

class directory_guard:
    def __init__(self, path):
        self._path = path

    def __enter__(self):
        self._old_dir = os.getcwd()
        os.chdir(self._path)
        return self

    def __exit__(self, *args):
        os.chdir(self._old_dir)


def base26_generator(maximum):
    if maximum == 1:
        yield 'a'
    iters = math.ceil(math.log(maximum, 26))
    for index in range(maximum):
        bases = [
            string.ascii_lowercase[(index % (26**(e + 1))) // (26**e)]
            for e in range(iters)
        ]
        bases.reverse()
        yield ''.join(bases)


def base49_generator(maximum):
    if maximum == 1:
        yield CAP_ALPHABET[0]
    iters = math.ceil(math.log(maximum, len(BITCOIN_ALPHABET)))
    for index in range(maximum):
        bases = [
            CAP_ALPHABET[(index % (len(CAP_ALPHABET)**(e + 1))) //
                             (len(CAP_ALPHABET)**e)] for e in range(iters)
        ]
        bases.reverse()
        yield ''.join(bases)

def base58_generator(maximum):
    if maximum == 1:
        yield BITCOIN_ALPHABET[0]
    iters = math.ceil(math.log(maximum, len(BITCOIN_ALPHABET)))
    for index in range(maximum):
        bases = [
            BITCOIN_ALPHABET[(index % (len(BITCOIN_ALPHABET)**(e + 1))) //
                             (len(BITCOIN_ALPHABET)**e)] for e in range(iters)
        ]
        bases.reverse()
        yield ''.join(bases)

def make_random_nonce():
    return base58.b58encode(os.urandom(8)).decode('utf-8')

