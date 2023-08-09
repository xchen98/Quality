import os
import pandas as pd

from libs.BUSCO.src.busco.run_BUSCO import run_busco
from libs.contamination import Contamination
from libs.basic_statistics import run_statistic
from libs.call_subprocess import call_bamtocov, call_calculate_cov, call_subprocess
from libs.alignment import get_protein_score, get_trans_score, check_alignment




def make_table(assembly_name, statistics, output_path, fragments, completeness, contamination, coverage, transcriptom, proteome):
    s_table = output_path + '/score_table.txt'
    score_table = open(s_table, 'w')
    score_table.write("Fundamental statistics of genome: " + assembly_name + '\n')
    score_table.write("---------------------------" + '\n')
    for i in statistics:
        score_table.write(str(i) + ": " + str(statistics[i]) + '\n')
    score_table.write("---------------------------" + '\n' + '\n')

    table = pd.DataFrame()
    table['fragments'] = fragments
    table['completeness'] = completeness['score'].tolist()
    table['contamination'] = list(contamination)
    table['transcriptome_score'] = transcriptom['score'].tolist()
    table['protein_score'] = proteome['score'].tolist()

    if coverage != None:
        table['coverage'] = coverage
        table['coverage'] = table['coverage'].astype(float)

        max_sum = max(table['completeness']) + max(table['contamination']) + max(table['transcriptome_score']) \
                  + max(table['protein_score']) + max(table['coverage'])
        table['score'] = [(table.loc[x]['completeness'] + table.loc[x]['contamination']
                           + table.loc[x]['protein_score'] + table.loc[x]['transcriptome_score']
                           + table.loc[x]['coverage']) / max_sum for x in
                          table.index.tolist()]

    else:
        max_sum = max(table['completeness']) + max(table['contamination']) + max(table['transcriptome_score'])\
                  + max(table['protein_score'])
        table['score'] = [(table.loc[x]['completeness'] + table.loc[x]['contamination']
                           + table.loc[x]['protein_score'] + table.loc[x]['transcriptome_score']) / max_sum for x in table.index.tolist()]

    table_string = table.to_string()
    score_table.write(table_string)

    return

def completeness(complete, fragments):
    complete = pd.DataFrame(complete, columns=['hmm_query', 'status', 'fragments', 'start', 'end', 'score'])
    complete = complete.sort_values(by=['fragments'])
    complete['score'] = complete['score'].astype(float)
    complete.start = complete.start.astype(int)
    complete.end = complete.end.astype(int)

    complete['length'] = abs((complete['end']) - (complete['start']))

    for i in set(complete.fragments):
        complete.loc[complete.fragments == i, 'final_score'] = (complete.loc[complete.fragments == i, 'score'] / complete.loc[
            complete.fragments == i, 'length']) * complete.loc[complete.fragments == i, 'length'].describe().mean()
    out = pd.DataFrame(fragments, columns=['fragments'])
    out['score'] = [0] * len(fragments)
    completeness = pd.DataFrame(complete['fragments'].unique(), columns=['fragments'])
    completeness['hmm_score'] = [complete.loc[complete['fragments'] == frag]['final_score'].describe().mean() for frag in completeness['fragments']]
    completeness['%complete_duplicated'] = [len(complete.loc[complete['fragments'] == frag]['hmm_query'].unique())/len(complete['hmm_query'].unique()) for frag in complete['fragments'].unique()]

    completeness['score'] = completeness['hmm_score'] * completeness['%complete_duplicated']
    completeness.round(2)

    for frag in completeness.fragments:
        if frag in list(completeness.fragments):
            out.loc[out.fragments == frag, 'score'] = completeness[completeness.fragments == frag].score.tolist()[0]

    return(out)

def run(assembly_name, output_path, assembly, lineage, bam, protein, transcriptome, annotation = None):

    auto_l = True
    if lineage:
        auto_l = False
    out = run_busco(assembly, output_path + '/busco', lineage, auto_l)

    fragments = call_subprocess(["grep '>'", assembly], out = True)
    fragments = [x[1:].split('.')[0] + '.' + x[1:].split('.')[1][0] for x in fragments]

    complete = completeness(out, fragments)

    Cont = Contamination(os.getcwd() + '/database', assembly)
    trim_name, trim_area = Cont.run(fragments)
    statistics, sequence, new_scaffolds_l, cleaned_scaffolds_l = run_statistic(assembly, trim_name, trim_area, output_path)

    Cont.contamination_score(new_scaffolds_l, cleaned_scaffolds_l)

    protein_score = get_protein_score(assembly, protein, fragments, output_path, annotation)
    transcriptome_score = get_trans_score(assembly, transcriptome, fragments, output_path, annotation)



    if bam == None:
        coverage = None
    else:
        coverage_dir = call_bamtocov(bam, output_path)
        coverage = [call_calculate_cov(x, coverage_dir) for x in fragments]


    make_table(assembly_name, statistics, output_path, fragments, complete, Cont.score, coverage, transcriptome_score, protein_score)



