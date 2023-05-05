import subprocess
import os
import numpy as np

def call_subprocess(command, out=False, string=None, contam=False, mode = None):
    if out:
        pip = subprocess.run(' '.join(command), shell=True, stdout=subprocess.PIPE)
        std = pip.stdout

        if (string != None) & (contam == False):
            std = str(std).replace('b','').replace(' ', '').split("'")[1].split('\\n')[:-1]
            std = [x.split(string)[0] for x in std]
        elif contam == True:
            std = str(std)[1:].split('# BLAST processed')[0].replace(' ', '').split("'")[1].split('#Query:')[1:]
            name = [x.split('\\n')[0] for x in std if "0hitsfound" not in x]
            #if mode == "vector":
            #    hits = [x.split('BLASTN')[0].split('\\n')[3:-1] for x in std if "0hitsfound" not in x]
            #    strong_hits = [[list(map(int, y.split('\\t')[8:10])) for y in x if np.floor(float(y.split('\\t')[-1]) / 2) >= 19] for x in hits]
            #    weak_hits = [[list(map(int, y.split('\\t')[8:10])) for y in x if np.floor(float(y.split('\\t')[-1]) / 2) < 19] for x in hits]
            #    for x in strong_hits:
            #        x.sort(key = lambda x: (x.sort(), x[0], x[1]))
            #
            #    return (name, strong_hits, weak_hits)
            hits = [x.split('BLASTN')[0].split('\\n')[4:-1] for x in std if "0hitsfound" not in x]
            strong_hits = [[list(map(int, y.split('\\t')[8:10])) for y in x] for x in hits]
            for x in strong_hits:
                x.sort(key = lambda x: (x.sort(), x[0], x[1]))

            return(name, strong_hits)
        else:
            if len(str(std).split('\\n')) > 1:
                std = str(std).replace('b', '').replace(' ', '').split("'")[1].split('\\n')[:-1]
            else:
                std = str(std).replace('b', '').replace(' ', '').split("'")[1]

        return (std)
    else:
        subprocess.run(' '.join(command), shell=True, check=True)


def call_bwa(assembly_path, reads_path, output_path, output_name):
    out = output_path + '/' + output_name + '.sam'
    command_index = [os.getcwd() + '/libs/bwa/bwa', 'index', assembly_path]
    command_mem = [os.getcwd() + '/libs/bwa/bwa', 'mem', assembly_path, reads_path, '>', out]
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
    command = ['bamtocov', bam_file, '>', out]
    call_subprocess(command)

    return out


def call_calculate_cov(fragment, coverage_file_dir):
    command = ['grep', fragment, coverage_file_dir, "| awk '{x+=$4} END {printf (\"%.2f\",x/NR)}' "]
    coverage = call_subprocess(command, out=True)
    return coverage
