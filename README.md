# bclaureate
Create (random) data sets matching Illumina sequencing file formats.

WIP.

Requires bgzip on the path.
https://github.com/samtools/htslib


Run with `$ ./bclaureate.py -m <machinetype>`

machinetype is one of:
    nextseq
    hiseqx
    hiseq4000
    hiseq2500
    miseq

Change number of lanes used, clusters per tile etc. by editing values in `PARAMS` dictionary near the top of `bclaureate.py` script.
