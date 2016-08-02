#!/usr/bin/python2.7

from __future__ import print_function
from random import randint
import os
from datetime import datetime
from Bio import bgzf

base_values = {
        'a': 0,
        'c': 1,
        'g': 2,
        't': 3
        }

machine_types = [
        'miseq',
        'hiseq',
        'nextseq'
        ]


class Run(object):
    def __init__(self, max_lanes, machinetype, machinename="testdata"):
        assert machinetype in machine_types
        assert type(machinename) is str and len(machinename) >= 1
        assert type(max_lanes) is int and max_lanes > 0 and max_lanes < 1000

        self.machinetype = machinetype
        self.machinename = machinename
        self.date = datetime.now().strftime("%y%m%d")
        self.dir_id = str(randint(0, 9999)).zfill(4)
        self.dir_name = '_'.join([self.date, self.machinename,
                                  self.dir_id, 'FC'])
        self.lanes = []
        for lane in xrange(1, max_lanes + 1):
            self.lanes.append(Lane(lane))

    def make_bcls(self):
        for lane in self.lanes:
            for tile in lane.tiles:
                if self.machinetype is 'nextseq':
                    f = bgzf.BgzfWriter(filename=os.path.join(lane.bcpath,
                                                              tile.name +
                                                              '.bcl.bgzf'),
                                        mode='wb')
                    f.write(tile.as_bcl())
                    f.close()
                else:
                    f = open(os.path.join(lane.bcpath, tile.name + '.bcl'),
                             'wb')
                    f.write(tile.as_bcl())
                    f.close()

    def make_bcis(self):
        for lane in self.lanes:
            f = open(os.path.join(lane.bcpath, 's_%d.bci' % lane.lane_idx),
                     'wb')
            f.write(lane.as_bci())
            f.close()


class Lane(object):
    def __init__(self, lane_idx, max_cycles=310, max_tiles=216):
        self.lane_idx = lane_idx
        self.name = "L" + str(lane_idx).zfill(3)
        self.cycles = []
        self.bcpath = ""
        for cycle in xrange(1, max_cycles + 1):
            self.cycles.append(Cycle(cycle))
        self.tiles = []
        for tile in xrange(1, max_tiles + 1):
            self.tiles.append(Tile(tile))

    def as_bci(self):
        l = []
        for tile in self.tiles:
            l += int_32_to_little_endian_list(tile.tile_idx)
            l += int_32_to_little_endian_list(tile.max_clusters)
        return bytearray(l)


class Cycle(object):
    def __init__(self, cycle_idx):
        self.cycle_idx = cycle_idx
        self.name = str(cycle_idx).zfill(4)


class Tile(object):
    def __init__(self, tile_idx, max_clusters=4000):
        self.tile_idx = tile_idx
        self.name = str(tile_idx + 11100)
        self.clusters = []
        self.max_clusters = max_clusters
        for cluster in xrange(max_clusters):
            self.clusters.append(Cluster(cluster))

    def as_bcl(self):
        l = int_32_to_little_endian_list(self.max_clusters)
        for c in self.clusters:
            l.append(c.as_byte())
        return bytearray(l)


class Cluster(object):
    def __init__(self, cluster_idx):
        self.cluster_idx = cluster_idx
        self.name = "C%d.1" % cluster_idx
        self.base = randint(0, 3)
        self.quality = randint(0, 63)

    def as_byte(self):
        return (self.quality << 2 | self.base)


def build_directory_structure(run):
    """Creates correct directory structure for run's machine type.
    """
    os.mkdir(os.path.join(os.getcwd(), run.dir_name))
    os.chdir(os.path.join(os.getcwd(), run.dir_name))
    for l in run.lanes:
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities', l.name))
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities',
                                 'BaseCalls', l.name))
        l.bcpath = os.path.join(os.getcwd(), 'Data', 'Intensities',
                                'BaseCalls', l.name)


def int_32_to_little_endian_list(n):
    l = []
    for i in xrange(0, 4):
        l.append((n >> (i * 8)) & 255)
    return l

if __name__ == "__main__":
    print("Running")
    run = Run(1, 'nextseq')
    build_directory_structure(run)
    run.make_bcls()
    run.make_bcis()
