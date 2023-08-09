import argparse
import sys
import os
import time
from libs.call_subprocess import call_samtools, call_bwa
from libs.reads import read_run
from base import run

def parser(args):
    parser = argparse.ArgumentParser(
        description=" ",
        usage="quality -a [GENOME ASSEMBLY] -l [LINEAGE] -o [OUTPUT_NAME] [OTHER OPTIONS]",
        add_help=False,
    )

    parser.add_argument(
        "-a",
        type=str,
        metavar="GENOME ASSEMBLY",
        help="Genome assembly"
    )

    parser.add_argument(
        "-reads",
        type=str,
        metavar="READS",
        help="Corresponding reads of given genome assembly"
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
        help="The given bam file helps calculate the coverage and accelerate the evaluation"
    )

    parser.add_argument(
        "-p",
        type=str,
        metavar="PROTEIN FILE",
        help="This file helps aligment"
    )

    parser.add_argument(
        "-t",
        type=str,
        metavar="TRANSCRIPTOME FILE",
        help="This file helps aligment"
    )

    parser.add_argument(
        "-o",
        type=str,
        metavar="OUTPUT_NAME",
        help="The output path"
    )

    parser.add_argument(
        "-anno",
        type=str,
        metavar="ANNOTATION_FILE",
        help="This file can help find annotation mistakes"
    )


    return (parser.parse_args())


# requirements:
#   install samtools, bamtocov, entrez direct, prefetch, fasterq-dump, makeblastdb, diamond, Rscript
def main(args):
    global auto_l
    args = parser(args)
    assembly = args.a
    reads = args.reads
    lineage = args.l
    bam = args.b
    protein = args.p
    trans = args.t
    anno = args.anno
    output_path = args.o
    if not output_path:
        output_path = './results/' + time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    assembly_name = assembly.split('/')[-1]

    accession = ('_').join(str(assembly).split('/')[-1].split('_')[:2])

    read = []
    if reads == 'auto':
        readsacc = read_run(accession, output_path)
        if readsacc != None:
            read = output_path + '/reads/' + readsacc[0] + '.fastq'

    if bam == None:
        if len(read) != 0:
            call_bwa(assembly, read, output_path, accession)
            bam = call_samtools(output_path + '/' + accession + '.sam', output_path, accession)
        else:
            bam = None

    run(assembly_name, output_path, assembly, lineage, bam, protein, trans, anno)


if __name__ == '__main__':
    results = main(sys.argv[1:])
    exit(results)
