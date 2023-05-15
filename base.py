import os
import numpy as np
from libs.BUSCO.src.busco.run_BUSCO import run_busco
from libs.contamination import Contamination
from libs.basic_statistics import run_statistic
from libs.call_subprocess import call_bamtocov, call_calculate_cov, call_subprocess


def make_table(fragments, output_path, completeness, contamination, coverage, transcriptom = None, proteome = None, multi = None, score = None):
    table = output_path + '/score_table.txt'
    score_table = open(table, 'w')
    column_name = ('\t').join(['fragments', 'completeness score', 'contamination score', 'coverage score', 'transcriptome score', 'proteome score', 'multi-assembly score', 'overall score'])
    score_table.write(column_name + '\n')

    for i in range(len(fragments)):
        line = ('\t').join([str(fragments[i]), str(round(np.float64(completeness[fragments[i]]), 2)), str(contamination[i]), str(coverage[i])])
        score_table.write(line + '\n')
    score_table.close()
    return



def run(output_path, assembly, lineage, bam, reads=None, reference=None):
    auto_l = True
    if lineage:
        auto_l = False
    out = run_busco(assembly, output_path + '/busco', lineage, auto_l)

    fragments = call_subprocess(["grep '>'", assembly], out = True)
    fragments = [x[1:] for x in fragments]


    Cont = Contamination(os.getcwd() + '/database', assembly)
    trim_name, trim_area = Cont.run(fragments)
    statistics, sequence, new_scaffolds_l, cleaned_scaffolds_l = run_statistic(assembly, trim_name, trim_area, output_path)
    Cont.contamination_score(new_scaffolds_l, cleaned_scaffolds_l)

    coverage_dir = call_bamtocov(bam, output_path)
    coverage = [call_calculate_cov(x, coverage_dir) for x in fragments]

    ini = [0] * len(fragments)
    complete = dict(zip(fragments, ini))
    freq = 1/len(out)
    for frag in out:
        if len(frag) > 1:
           if frag[0] == 'Complete':
               complete[frag[1]] = complete[frag[1]] + freq
           if frag[0] == 'Duplicated':
               complete[frag[1]] = complete[frag[1]] + freq

    make_table(fragments, output_path, complete, Cont.score, coverage)


    return