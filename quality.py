import argparse
import sys
import os
import time
from libs.call_subprocess import call_samtools, call_bwa
from libs.reads import read_run
from base import run
from libs.preprocessing import load_fragments

def parser(args):
    parser = argparse.ArgumentParser(
        description=" ",
        usage="quality -a [GENOME ASSEMBLY] -ref [REFERENCE] -o [OUTPUT_NAME] [OTHER OPTIONS]",
        add_help=False,
    )

    parser.add_argument(
        "-a",
        type=str,
        metavar="GENOME ASSEMBLY",
        help="genome assembly"
    )
    parser.add_argument(
        "-ref",
        type=str,
        metavar="REFERENCE",
        help="reference genome"
    )
    parser.add_argument(
        "-reads",
        type=str,
        metavar="READS",
        help="corresponding reads of given genome assembly"
    )
    parser.add_argument(
        "-l",
        type=str,
        metavar="LINEAGE",
        help="Specify the name of the BUSCO lineage to be used."
    )
    parser.add_argument(
        "-o",
        type=str,
        metavar="OUTPUT_NAME",
        help="the output path"
    )

    parser.add_argument(
        "-b",
        type=str,
        metavar="BAM FILE",
        help="the given bam file helps calculate the coverage"
    )

    return (parser.parse_args())


# requirements:
#   install samtools, bamCoverage, entrez direct, prefetch, fasterq-dump
def main(args):
    global auto_l
    args = parser(args)
    assembly = args.a
    reads = args.reads
    reference = args.ref
    lineage = args.l
    bam = args.b
    output_path = args.o
    if not output_path:
        output_path = './results/' + time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        if not os.path.exists(output_path):
            os.makedirs(output_path)
    # fragment = call_minimap2(assembly, reference, output_path, "fragments")
    #output_name = 'fragment'

    accession = ('_').join(str(assembly).split('/')[-1].split('_')[:2])
    #inputdir = os.path.abspath(os.path.dirname(os.path.abspath(assembly)) + os.path.sep + '.')
    if reads == None:
        readsacc = read_run(accession, output_path)
        reads = output_path + '/' + readsacc[0]
    if bam == None:
        call_bwa(assembly, reads, output_path, accession)
        bam = call_samtools(output_path + '/' + accession + '.sam', output_path, accession)

    run(output_path, assembly, lineage, bam)

if __name__ == '__main__':
    results = main(sys.argv[1:])
    exit(results)
