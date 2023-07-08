import argparse
import sys
import os
import time
from libs.call_subprocess import call_samtools, call_bwa, call_subprocess
from libs.reads import read_run
from base import run
from libs.call_subprocess import call_bamtocov, call_calculate_cov
from libs.contamination import Contamination
from libs.alignment import get_trans_score

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
        "-b",
        type=str,
        metavar="BAM FILE",
        help="the given bam file helps calculate the coverage"
    )

    parser.add_argument(
        "-p",
        type=str,
        metavar="PROTEIN FILE",
        help="this file helps aligment"
    )

    parser.add_argument(
        "-t",
        type=str,
        metavar="TRANSCRIPTOME FILE",
        help="this file helps aligment"
    )

    parser.add_argument(
        "-o",
        type=str,
        metavar="OUTPUT_NAME",
        help="the output path"
    )

    parser.add_argument(
        "-anno",
        type=str,
        metavar="ANNOTATION_FILE",
        help="This file can help find annotation mistakes"
    )


    return (parser.parse_args())


# requirements:
#   install samtools, bamtocov, entrez direct, prefetch, fasterq-dump, makeblastdb, kalign,
def main(args):
    global auto_l
    args = parser(args)
    assembly = args.a
    reads = args.reads
    reference = args.ref
    lineage = args.l
    bam = args.b
    protein = args.p
    trans = args.t
    output_path = args.o
    if not output_path:
        output_path = './results/' + time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    assembly_name = assembly.split('/')[-1]

    accession = ('_').join(str(assembly).split('/')[-1].split('_')[:2])
    #inputdir = os.path.abspath(os.path.dirname(os.path.abspath(assembly)) + os.path.sep + '.')


    if reads == None:
        readsacc = read_run(accession, output_path)
        if readsacc == None:
            reads = []
        else:
            reads = output_path + '/reads/' + readsacc[0] + '.fastq'
    if bam == None:
        if len(reads) != 0:
            call_bwa(assembly, reads, output_path, accession)
            bam = call_samtools(output_path + '/' + accession + '.sam', output_path, accession)
        else:
            bam = None
    #print(output_path)
    fragments = call_subprocess(["grep '>'", assembly], out = True)
    fragments = [x[1:].split('.')[0] + '.' + x[1:].split('.')[1][0] for x in fragments]
    coverage_dir = call_bamtocov(bam, output_path)
    coverage = [call_calculate_cov(x, coverage_dir) for x in fragments]
    #print(coverage)
    run(assembly_name, output_path, assembly, lineage, bam, protein, trans, None)


if __name__ == '__main__':
    results = main(sys.argv[1:])
    exit(results)
    #time_start = time.time()

    #run_aligment("../Downloads/GCF_000001735.4_TAIR10.1_protein.faa", './seq_X', "./kalign")
    #time_end = time.time()
    #print('time:', time_end - time_start)