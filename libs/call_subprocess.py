import subprocess
import os
import pandas as pd

def call_subprocess(command, out=False, string=None, contam=False, align = False):
    if out:
        pip = subprocess.run(' '.join(command), shell=True, stdout=subprocess.PIPE)
        std = pip.stdout

        if (string != None) & (contam == False):
            std = str(std)[1:].replace(' ', '').split("'")[1].split('\\n')[:-1]
            std = [x.split(string)[0] for x in std]
        if align == True:

            std = str(std)[1:].replace(' ', '').split("'")[1].split('\\n')
            align = pd.DataFrame([x.split('\\t') for x in std if '#' not in x and len(x) != 0])
            align.columns = ['queryacc.ver', 'subjectacc.ver', 'identity', 'alignment length', 'mismatches', 'gapopens',
                         'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bitscore']
            return align
        elif contam == True:
            std = str(std)[1:].split('# BLAST processed')[0].replace(' ', '').split("'")[1].split('#Query:')[1:]
            name = [x.split('\\n')[0].split('.')[0] + '.' + x.split('\\n')[0].split('.')[1][0] for x in std if "0hitsfound" not in x]
            hits = [x.split('BLASTN')[0].split('\\n')[4:-1] for x in std if "0hitsfound" not in x]
            strong_hits = [[list(map(int, y.split('\\t')[8:10])) for y in x] for x in hits]
            for x in strong_hits:
                x.sort(key = lambda x: (x.sort(), x[0], x[1]))

            return(name, strong_hits)
        else:
            if len(str(std).split('\\n')) > 1:
                std = str(std)[1:].replace(' ', '').split("'")[1].split('\\n')[:-1]
            else:
                std = str(std)[1:].replace(' ', '').split("'")[1]
                std = str(std)[1:].replace(' ', '').split("'")[1]

        return (std)
    else:
        subprocess.run(' '.join(command), shell=True, check=True)


def call_bwa(assembly_path, reads_path, output_path, output_name):
    out = output_path + '/' + output_name + '.sam'
    command_index = [os.getcwd() + '/libs/bwa/bwa', 'index', assembly_path]
    command_mem = [os.getcwd() + '/libs/bwa/bwa', 'mem -t 8', assembly_path, reads_path, '>', out]
    call_subprocess(command_index)
    call_subprocess(command_mem)
    return (out)


def call_samtools(sam_file, output_path, output_name):
    out = output_path + '/' + output_name + '.bam'
    command_sort = ['samtools sort', sam_file, '>', out]
    command_index = ['samtools index', out]
    call_subprocess(command_sort)
    call_subprocess(command_index)
    return (out)


def call_bamtocov(bam_file, output_path):
    out = output_path + '/' + 'coverage.bed'
    command = ['bamtocov -T 8', bam_file, '>', out]
    call_subprocess(command)

    return out


def call_calculate_cov(fragment, coverage_file_dir):
    command = ['grep', fragment, coverage_file_dir, "| awk '{x+=$4} END {printf (\"%.2f\",x/NR)}' "]
    coverage = call_subprocess(command, out=True)
    return coverage

def call_makeblastdb(seq, type, output_dir):
    command = ['makeblastdb', '-in', seq, '-dbtype', type, '-out', output_dir + '/' + str(type) + '/' + str(type), "-parse_seqids"]
    call_subprocess(command)

    return(output_dir + '/' + str(type) + '/' + str(type))


def call_diamond(seq, protein, output_dir):
    command = ['diamond makedb', '--in', protein, '-d', output_dir + '/nr', '>', output_dir + '/log.txt']
    call_subprocess(command)

    command = ['diamond blastx -d', output_dir + '/nr', '-q', seq, '-f 6 -o', output_dir + '/matches.m8', '>', output_dir + '/log.txt']
    call_subprocess(command)

