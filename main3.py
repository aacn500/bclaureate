#!/usr/bin/python2.7

from __future__ import print_function
import random
import os
from datetime import datetime
from xml.etree import ElementTree
from xml.dom import minidom
import struct
import subprocess


class RunParameters(object):
    def __init__(self):
        self.lanes = 4
        self.surfaces = 2
        self.swaths = 3
        self.tiles = 12
        self.sections = 3
        self.max_tiles = 216
        self.clusters = 2000
        self.machinetype = 'nextseq'
        self.machinename = 'NSTESTMACHINE'
        self.flowcellname = "TESTFLOWCELL"
        self.date = datetime.now().strftime("%y%m%d")
        self.reads = [{"num_cycles": 1, "is_indexed": False},
                      {"num_cycles": 1, "is_indexed": False}]

PARAMS = RunParameters()

class Run(object):
    def __init__(self):
        self.lanes = [Lane() for lane in xrange(PARAMS.lanes)]
        self.reads = [Read(read["num_cycles"], read["is_indexed"])
                      for read in PARAMS.reads]
        self.id = "{}_{}_{:04d}".format(PARAMS.date, PARAMS.machinename,
                                    random.randint(0, 9999))
        self.infopath = ""

    def make_runinfo(self):
        root = ElementTree.Element("RunInfo", {
                "xmlns:xsd": "http://www.w3.org/2001/XMLSchema",
                "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
                "Version": "2"
                })

        run = ElementTree.SubElement(root, "Run", {
                "Id": "{}".format(self.id),
                "Number": "2"
                })

        flowcell = ElementTree.SubElement(root, "Flowcell")
        flowcell.text = PARAMS.flowcellname

        instrument = ElementTree.SubElement(root, "Instrument")
        instrument.text = PARAMS.machinename

        date = ElementTree.SubElement(run, "Date")
        date.text = PARAMS.date

        reads = ElementTree.SubElement(run, "Reads")
        for i in xrange(len(self.reads)):
            read = ElementTree.SubElement(reads, "Read", {
                    "Number": str(i+1),
                    "NumCycles": "{:d}".format(self.reads[i].num_cycles),
                    "IsIndexedRead": "Y" if self.reads[i].is_indexed
                                     else "N"
                    })

        flowcell = ElementTree.SubElement(run, "FlowcellLayout", {
            "LaneCount": "{:d}".format(PARAMS.lanes),
            "SurfaceCount": "{:d}".format(PARAMS.surfaces),
            "SwathCount": "{:d}".format(PARAMS.swaths),
            "TileCount": "{:d}".format(PARAMS.tiles),
            "SectionPerLane": "{:d}".format(PARAMS.sections),
            "LanePerSection": "2"
            })

        tileset = ElementTree.SubElement(flowcell, "TileSet", {
            "TileNamingConvention": "FiveDigit"})

        tiles = ElementTree.SubElement(tileset, "Tiles")

        for lane_idx in xrange(PARAMS.lanes):
            for section_idx in xrange(PARAMS.sections):
                for swath_idx in xrange(PARAMS.swaths):
                    for surface_idx in xrange(PARAMS.surfaces):
                        for tile_idx in xrange(PARAMS.tiles):
                            section_offset = 1 if lane_idx < 2 else 4
                            tile = ElementTree.SubElement(tiles, "Tile")
                            tile.text = "{:d}_{:d}{:d}{:d}{:02d}".format(
                                    lane_idx + 1,
                                    surface_idx + 1,
                                    swath_idx + 1,
                                    section_idx + section_offset,
                                    tile_idx + 1)

        image_dims = ElementTree.SubElement(run, "ImageDimensions", {
            "Width": "2592",
            "Height": "1944"
            })

        image_chans = ElementTree.SubElement(run, "ImageChannels")

        red = ElementTree.SubElement(image_chans, "Name")
        red.text = "Red"

        green = ElementTree.SubElement(image_chans, "Name")
        green.text = "Green"

        f = open(os.path.join(self.infopath, 'RunInfo.xml'), 'w')
        f.write(pp(root))
        f.close()

    def make_bcls(self):
        for lane in self.lanes:
            cycle_idx = 0
            for read in self.reads:
                for cycle in xrange(read.num_cycles):
                    cycle_idx += 1
                    f = open(os.path.join(lane.bcpath,
                        "{:04d}.bcl".format(cycle_idx)), mode='wb')
                    f.write(lane.bcl())
                    f.close()
                    subprocess.call(["/software/solexa/pkg/tabix/current/bgzip",
                                     os.path.join(lane.bcpath,
                                         "{:04d}.bcl".format(cycle_idx))])
                    subprocess.call(["mv", os.path.join(lane.bcpath,
                                        "{:04d}.bcl.gz".format(cycle_idx)),
                                     os.path.join(lane.bcpath, 
                                        "{:04d}.bcl.bgzf".format(cycle_idx))])

    def make_bcis(self):
        cn = struct.pack("<I", PARAMS.clusters)
        for lane_idx in xrange(PARAMS.lanes):
            lane = self.lanes[lane_idx]
            f = open(os.path.join(lane.bcpath,
                                  's_{:d}.bci'.format(lane_idx + 1)), 'wb')
            s = ''
            for section_idx in xrange(PARAMS.sections):
                for swath_idx in xrange(PARAMS.swaths):
                    for surface_idx in xrange(PARAMS.surfaces):
                        for tile_idx in xrange(PARAMS.tiles):
                            section_offset = 1 if lane_idx < 2 else 4
                            s += struct.pack('<I',
                                  int("{:d}{:d}{:d}{:02d}".format(
                                    surface_idx + 1,
                                    swath_idx + 1,
                                    section_idx + section_offset,
                                    tile_idx + 1)))
                            s += cn
            f.write(s)
            f.close()


    def make_filters(self):
        for lane_idx in xrange(PARAMS.lanes):
            lane = self.lanes[lane_idx]
            s = struct.pack("<I", 0)
            s += struct.pack("<I", 3)
            cluster_passes = ""
            c = 0
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            for cluster in tile.clusters:
                                c += 1
                                cluster_passes += struct.pack("B", cluster.filtered)
            s += struct.pack("<I", c)
            s += cluster_passes
            f = open(os.path.join(lane.bcpath,
                                  "s_{:d}.filter".format(lane_idx + 1)), "wb")
            f.write(s)
            f.close()

    def make_locs(self):
        for lane_idx in xrange(PARAMS.lanes):
            lane = self.lanes[lane_idx]
            s = struct.pack('<I', 1)
            s += struct.pack('<f', 1.0)
            cluster_locs = ""
            c = 0
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            for cluster in tile.clusters:
                                c += 1
                                cluster_locs += struct.pack("<f", cluster.x_coord)
                                cluster_locs += struct.pack("<f", cluster.y_coord)
            s += struct.pack("<I", c)
            s += cluster_locs
            
            f = open(os.path.join(lane.locspath, "s_{:d}.locs".format(lane_idx + 1)),
                    "wb")
            f.write(s)
            f.close()
                            


class Read(object):
    def __init__(self, num_cycles, is_indexed):
        self.num_cycles = num_cycles
        self.is_indexed = is_indexed

class Lane(object):
    def __init__(self):
        self.sections = [Section() for section in xrange(PARAMS.sections)]
        self.bcpath = ""

    def bcl(self):
        s = ""
        for section in self.sections:
            for swath in section.swaths:
                for surface in swath.surfaces:
                    for tile in surface.tiles:
                        for cluster in tile.clusters:
                            s += cluster.byte

        l = len(s)
        return struct.pack("<I", l) + s


class Cycle(object):
    def __init__(self):
        pass

class Section(object):
    def __init__(self):
        self.swaths = [Swath() for swath in xrange(PARAMS.swaths)]

class Swath(object):
    def __init__(self):
        self.surfaces = [Surface() for surface in xrange(PARAMS.surfaces)]


class Surface(object):
    def __init__(self):
        self.tiles = [Tile() for tile in xrange(PARAMS.tiles)]


class Tile(object):
    def __init__(self):
        self.clusters = [Cluster() for cluster in xrange(PARAMS.clusters)]

    def as_bcl(self):
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

def build_directory_structure(run):
    os.mkdir(os.path.join(os.getcwd(), run.id + "_FC"))
    os.chdir(os.path.join(os.getcwd(), run.id + "_FC"))
    run.infopath = os.getcwd()
    for l in xrange(len(run.lanes)):
        lane = run.lanes[l]
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities',
                                 'L{:03d}'.format(l + 1)))
        lane.locspath = os.path.join(os.getcwd(), 'Data', 'Intensities',
                                  'L{:03d}'.format(l + 1))
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities',
                                  'BaseCalls', 'L{:03d}'.format(l + 1)))
        lane.bcpath = os.path.join(os.getcwd(), 'Data', 'Intensities',
                                  'BaseCalls', 'L{:03d}'.format(l + 1))

def pp(elem):
    """Return a pretty-printed XML string for given element.
    Taken from:
    https://pymotw.com/2/xml/etree/ElementTree/create.html
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


if __name__ == "__main__":
    run = Run()
    build_directory_structure(run)
    run.make_runinfo()
    run.make_bcls()
    run.make_bcis()
    run.make_filters()
    run.make_locs()
