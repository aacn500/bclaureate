#!/usr/bin/python2.7

from __future__ import print_function
import random
import os
from datetime import datetime
from Bio import bgzf
import struct


BASE_VALUES = {
        'a': 0,
        'c': 1,
        'g': 2,
        't': 3
        }

MACHINE_TYPES = [
        'miseq',
        'hiseq',
        'nextseq'
        ]

LANES = 4
SURFACES = 2
SWATHS = 3
SEGMENTS = 3
CAMERAS = SEGMENTS * (LANES/2)  # 6
TILES = 12
MAX_CYCLES = 3
MAX_TILES = 216
CLUSTERS = 100


class Run(object):
    def __init__(self, reads, machinetype, machinename="testdata"):
        assert machinetype in MACHINE_TYPES
        assert type(machinename) is str and len(machinename) >= 1

        self.machinetype = machinetype
        self.machinename = machinename
        self.date = datetime.now().strftime("%y%m%d")
        self.dir_id = str(random.randint(0, 9999)).zfill(4)
        self.dir_name = '_'.join([self.date, self.machinename,
                                  self.dir_id, 'FC'])
        self.lanes = []
        for lane in xrange(1, LANES + 1):
            self.lanes.append(Lane(lane))
        self.reads = []
        for read in reads:
            self.reads.append(read)

    def add_read(self, read):
        assert type(read) is Read
        self.reads.append(read)

    def make_bcls(self):
        """These files contain the machine's base calls and quality scores for
        a given tile.
        bcl files from a nextseq machine are zipped in the blocked
        GNU zip format bgzf. All others are zipped as gzips.
        """
        for lane in self.lanes:
            for tile in lane.tiles:
                if self.machinetype is 'nextseq':
                    f = bgzf.BgzfWriter(filename=os.path.join(lane.bcpath,
                                                              lane.prefix + 
                                                              tile.name +
                                                              '.bcl.bgzf'),
                                        mode='wb')
                    f.write(tile.as_bcl())
                    f.close()
                else:
                    f = open(os.path.join(lane.bcpath, lane.prefix + 
                                          tile.name + '.bcl'),
                             'wb')
                    f.write(tile.as_bcl())
                    f.close()

    def make_bcis(self):
        """This file contains a record of how many clusters there are on a
        given tile for a lane.
        """
        for lane in self.lanes:
            f = open(os.path.join(lane.bcpath, 's_%d.bci' % lane.lane_idx),
                     'wb')
            f.write(lane.as_bci())
            f.close()

    def make_filters(self):
        """This file records which clusters passed the filter for a given lane.
        """
        for lane in self.lanes:
            f = open(os.path.join(lane.bcpath, "s_%d.filter" % lane.lane_idx),
                     'wb')
            f.write(lane.as_filter())
            f.close()

    def make_locs(self):
        """This file records the locations of clusters for a given lane.
        """
        for lane in self.lanes:
            f = open(os.path.join(lane.locspath, "s_%d.locs" % lane.lane_idx),
                     'wb')
            f.write(lane.as_locs())
            f.close()


class Read(object):
    def __init__(self, num_cycles, is_indexed):
        self.num_cycles = num_cycles
        self.is_indexed = is_indexed


class Lane(object):
    def __init__(self, lane_idx):
        self.lane_idx = lane_idx
        self.name = "L" + str(lane_idx).zfill(3)
        self.cycles = []
        self.bcpath = ""
        self.locspath = ""
        self.clusters = []
        self.prefix = "{:d}_".format(lane_idx)
        for cycle in xrange(1, MAX_CYCLES + 1):
            self.cycles.append(Cycle(cycle))
        self.tiles = []
        for tile in xrange(1, MAX_TILES + 1):
            self.tiles.append(Tile(tile, self))

    def as_bci(self):
        """Binary file containing the number of a tile and the number of 
        clusters upon it. 
            Bytes 0-3: tile number
            Bytes 4-7: number of clusters on tile.
        """
        l = []
        for tile in self.tiles:
            l += int_32_to_little_endian_list(tile.tile_idx)
            l += int_32_to_little_endian_list(tile.CLUSTERS)
        return bytearray(l)

    def as_filter(self):
        """Binary file recording which clusters successfully passed the filter.
            Bytes 0-3: zero value
            Bytes 4-7: Version of filter file format
            Bytes 8-11: Number of clusters
            Bytes 12-X: 8 bit integer: 0 if failed filter, else 1.
        """
        l = [0 for n in xrange(4)]
        l += int_32_to_little_endian_list(3)
        l += int_32_to_little_endian_list(len(self.clusters))
        for cluster in self.clusters:
            l.append(cluster.filtered)
        return bytearray(l)

    def as_locs(self):
        """Binary file containing each cluster's location on its tile for a 
        given lane.
            Bytes 1-4: Int equalling 1
            Bytes 5-8: Float equalling 1
            Bytes 9-12: unsigned int: Number of clusters
            Bytes 13-16: 32-bit float: X-coordinate of first cluster
            Bytes 17-20: 32-bit float: Y-coordinate of first cluster
            Following bytes provide the X and Y coordinates of following
            clusters.
        """
        s = struct.pack('i', 1)
        s += struct.pack('f', 1.0)
        s += struct.pack('I', len(self.clusters))
        for cluster in self.clusters:
            s += struct.pack('f', cluster.x_coord)
            s += struct.pack('f', cluster.y_coord)
        return s


class Cycle(object):
    def __init__(self, cycle_idx):
        self.cycle_idx = cycle_idx
        self.name = str(cycle_idx).zfill(4)


class Tile(object):
    def __init__(self, tile_idx, lane, surface, swath, camera,
                 max_clusters=CLUSTERS):
        self.lane = lane
        self.tile_idx = tile_idx
        self.surface = surface
        self.swath = swath
        self.camera = camera
        self.name = "{:d}{:d}{:d}{:02d}".format(surface, swath, camera, tile)
        self.clusters = []
        self.max_clusters = max_clusters
        for cluster in xrange(max_clusters):
            cl = Cluster(cluster)
            self.clusters.append(cl)
            self.lane.clusters.append(cl)

    def as_bcl(self):
        l = int_32_to_little_endian_list(self.max_clusters)
        for c in self.clusters:
            l.append(c.as_byte())
        return bytearray(l)


class Cluster(object):
    def __init__(self, cluster_idx):
        self.cluster_idx = cluster_idx
        self.name = "C%d.1" % cluster_idx
        self.base = random.randint(0, 3)
        self.quality = random.randint(0, 63)
        self.filtered = random.randint(0, 1)
        self.x_coord = random.uniform(-10.0, 9999999.99)
        self.y_coord = random.uniform(-10.0, 9999999.99)

    def as_byte(self):
        return (self.quality << 2 | self.base)


def build_directory_structure(run):
    """Creates correct directory structure for run's machine type.
    """
    os.mkdir(os.path.join(os.getcwd(), run.dir_name))
    os.chdir(os.path.join(os.getcwd(), run.dir_name))
    for l in run.lanes:
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities', l.name))
        l.locspath = os.path.join(os.getcwd(), 'Data', 'Intensities', l.name)
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
    READS = [Read(151, False), Read(8, True), Read(151, False)]
    print("Running")
    run = Run(READS, 1, 'nextseq')
    run.add_read(Read(151, False))
    run.add_read(Read(8, True))
    run.add_read(Read(151, False))
    build_directory_structure(run)
    run.make_bcls()
    run.make_bcis()
    run.make_filters()
    run.make_locs()
