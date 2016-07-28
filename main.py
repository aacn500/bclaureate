#!/usr/bin/python
from datetime import datetime
from xml.etree import ElementTree
from xml.dom import minidom
import random
import os

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
    def __init__(self, machinetype, machinename="testdata"):
        assert machinetype in machine_types
        self.machinetype = machinetype
        self.lanes = []
        self.machinename = machinename
        self.date = datetime.now().strftime("%y%m%d")
        self.reads = []
        self.cycles = []

    def gen_run_dir_name(self):
        return self.date +\
                "_" + self.machinename + "_" +\
                str(random.randint(0, 9999)).zfill(4) +\
                "_FC"

    def add_lane(self):
        self.lanes.append(Lane('L' + str(len(lanes)+1).zfill(3)))

    def add_read(self, is_indexed_read=False):
        self.reads.append(Read(self, is_indexed_read))

    def create_miseq_runinfo_xml(self):
        root = ElementTree.Element('RunInfo',
                {
                    'xmlns:xsd': "http://www.w3.org/2001/XMLSchema",
                    'xmlns:xsi': "http://www.w3.org/2001/XMLSchema-instance",
                    'Version': "2"
                })

        run = ElementTree.SubElement(root, 'Run', 
                {
                    'Id': "LongComplexIdWithLotsOfNumbers",
                    'Number': "1"
                })

        flowcell = ElementTree.SubElement(run, 'Flowcell')
        flowcell.text = "000000000-TEXTY"

        instrument = ElementTree.SubElement(run, 'Instrument')
        instrument.text = "MSINSTRU"

        date = ElementTree.SubElement(run, 'Date')
        date.text = self.date

        reads = ElementTree.SubElement(run, 'Reads')

        for i in xrange(len(self.reads)):
            read = ElementTree.SubElement(reads, 'Read', 
                    {
                        "NumCycles": str(self.reads[i].num_cycles),
                        "Number": str(i+1),
                        "IsIndexedRead": \
                            'Y' if self.reads[i].is_indexed_read else 'N'
                    })
        flowcell = ElementTree.SubElement(run, 'FlowcellLayout',
                {
                    "LaneCount": "1",
                    "SurfaceCount": "2",
                    "SwathCount": "1",
                    "TileCount": "14"
                })

        f = open("tmp/RunInfo.xml", 'w')
        f.write(pp(root)+'\n')
        f.close()


class Lane(object):
    def __init__(self, name):
        self.tiles = []
        self.name = name

    def add_tile(self):
        self.tiles.append(Tile())


class Tile(object):
    def __init__(self):
        self.clusters = []


class Read(object):
    def __init__(self, run, is_indexed_read=False):
        assert type(run) is Run
        assert type(is_indexed_read) is bool
        self.run = run
        self.cycles = []
        self.num_cycles = 0
        self.is_indexed_read = is_indexed_read

    def add_cycle(self):
        self.cycles.append(Cycle())
        self.run.cycles.append(Cycle())
        self.num_cycles += 1


class Cycle(object):
    def __init__(self):
        self.cluster_count = 5000

    def create_bcl(self):
        f = open('tmp/s_1_1101.bcl', 'wb+')
        l = int_32_to_little_endian(self.cluster_count)
        for n in xrange(self.cluster_count):
            l.append(Cluster().to_byte())
        b = bytearray(l)
        f.write(b)
        f.write('\n')
        f.close()


class Cluster(object):
    def __init__(self):
        self.base = random.randint(0, 3)
        self.quality = random.randint(0, 63)

    def to_byte(self):
        assert self.base >= 0 and self.base <= 3
        assert self.quality >= 0 and self.quality <= 63
        return (self.quality << 2 | self.base)


def cluster_to_byte(base, quality):
    assert base in base_values.keys()
    quality = clamp(quality, 0, 63)
    quality <<= 2
    bin_quality = bin(quality)
    return bin(quality | base_values[base])[2:]


def clamp(n, minn, maxn):
    return min(max(n, minn), maxn)


def int_to_32b_little_endian(n):
    c = bin(n)
    c = c[2:].zfill(32)
    l = [c[:8], c[8:16], c[16:24], c[24:32]]
    l.reverse()
    s = ""
    for c in l:
        s += chr(int(c, 2))
    return s

def int_32_to_little_endian(n):
    l = []
    for i in xrange(0, 4):
        l.append((n >> (i * 8)) & 255)
    return l


def pp(elem):
    """Return a pretty-printed XML string for the Element.
    Taken from:
    https://pymotw.com/2/xml/etree/ElementTree/create.html
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


if __name__ == "__main__":
#    run = Run()
#    print run.gen_run_dir_name()
#    os.mkdir(run.gen_run_dir_name())
#    print int_to_bin_chr(1094861636)
#    print cluster_to_byte('c', 31)
    cycle = Cycle()
    cycle.create_bcl()
#    run = Run('miseq')
#    run.add_read()
#    run.add_read(True)
#    run.add_read()
#    run.add_read()
#    run.create_miseq_runinfo_xml()
