import struct
import sys

def int_to_bytes(n):
    s = bin(n)[2:]
    l = [s[:8], s[8:16], s[16:24], s[24:]]
    return ''.join([chr(int(c, 2)) for c in l])

def intrep_to_float(n):
    return struct.unpack('f', int_to_bytes(n))[0]


if __name__ == "__main__":
#    file_name = sys.argv[1]
#    f = open('example.locs', 'r')
#    li = f.readlines()
#    f.close()
    li = sys.stdin.readlines()
    for line in xrange(len(li)):
        li[line] = li[line].strip().split('      ')
        for coord in xrange(len(li[line])):
            li[line][coord] = intrep_to_float(int(li[line][coord]))
        print li[line]

