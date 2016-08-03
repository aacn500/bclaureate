#!/usr/bin/python2.7

from __future__ import print_function
import random
import os
from datetime import datetime
from Bio import bgzf
import struct

LANES = 4
SURFACES = 2
SWATHS = 3
TILES = 12
SECTIONS = 3
CYCLES = 100
TILES = 216
CLUSTERS = 2000

class RunParameters(object):
    def __init__(self):
        self.lanes = 4
        self.surfaces = 2
        self.swaths = 3
        self.tiles = 12
        self.sections = 3
        self.cycles = 100
        self.max_tiles = 216
        self.clusters = 2000
        self.machinetype = 'nextseq'
        self.machinename = 'testdata'
        self.reads = [{"num_cycles": 151, "is_indexed": False},
                      {"num_cycles": 8, "is_indexed": True},
                      {"num_cycles": 151, "is_indexed": False}]

PARAMS = RunParameters()

class Run(object):
    def __init__(self):
        self.lanes = [Lane() for lane in xrange(PARAMS.lanes)]
        self.reads = [Read(read["num_cycles"], read["is_indexed"])
                      for read in PARAMS.reads]


class Read(object):
    def __init__(self):
        self.cycles = [Cycle() for cycle in xrange(PARAMS.cycles)]

class Lane(object):
    def __init__(self):
        self.sections = [Section() for section in xrange(PARAMS.sections)]

class Cycle(object):
    def __init__(self):
        pass

class Section(object):
    def __init__(self):
        self.swaths = [Swath() for swath in xrange(PARAMS.swaths)]

class Swath(object):
    def __init__(self):
        self.tiles = [Tile() for tile in xrange(PARAMS.tiles)]

class Tile(object):
    def __init__(self):
        self.clusters = [Cluster() for cluster in xrange(PARAMS.clusters)]

    def as_bcl(self):
        s = struct.pack('<I', len(self.clusters))
        for cluster in self.clusters:
            s += cluster.byte
        return s

class Cluster(object):
    def __init__(self):
        self.byte = chr(random.randint(0, 255))
        self.filtered = random.randint(0, 1)
        self.x_coord = random.uniform(-10.0, 9999999.99)
        self.y_coord = random.uniform(-10.0, 9999999.99)


def int_32_to_little_endian_list(n):
    l = []
    for i in xrange(0, 4):
        l.append((n >> (i * 8)) & 255)
    return l
