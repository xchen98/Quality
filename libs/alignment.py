import os
import numpy as np
import pandas as pd

from Bio.Seq import Seq
from libs.call_subprocess import call_makeblastdb, call_subprocess, call_kalign, call_hmmbuild, call_hmmsearch, call_diamond


def alignment_transcriptome(sequences, transcriptome, output_dir):
    transcriptome_dir = call_makeblastdb(transcriptome, 'nucl', output_dir)
    command = ['blastn', '-query', sequences, '-db', transcriptome_dir,
               '-evalue 10e-6 -num_threads 8 -outfmt 7']

    out = call_subprocess(command, out=True, align=True)
    return(out)


def get_protein_score(seq, protein, fragments, output_dir):
    call_diamond(seq, protein, output_dir)
    matching = pd.read_csv(output_dir + '/matches.m8', sep="\t")
    matching.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                        'send', 'evalue', 'bitscore']

    #score = []
    #for fragment in fragments:
    #    sorted_matching = matching[matching['qseqid'] == fragment].sort_values(by=['bitscore', 'length', 'pident'],
    #                                                                           ascending=False)
    #    score_f = pd.DataFrame(sorted_matching['bitscore'] / sorted_matching['length'] * sorted_matching['pident'],
    #                           columns=['score'])
    #    score.append([sorted_matching.loc[
    #                      score_f.sort_values(by=['score'], ascending=False).iloc[0:1].index.tolist()[0]]['sseqid'],
    #                  score_f.sort_values(by=['score'], ascending=False).max().max()])
    #score = pd.DataFrame(score, columns=['protein_sseqid', 'score'])
    #score.index = fragments

    score = pd.DataFrame()
    for fragment in fragments:
        f_matching = matching[matching['qseqid'] == fragment]
        score_f = pd.DataFrame([fragment, np.mean(
            f_matching['bitscore'].astype(float) / f_matching['length'].astype(float) * f_matching[
                'pident'].astype(float))]).T
        score = pd.concat([score, score_f])
    score.columns = ['fragments', 'score']
    score.index = fragments
    return score

def get_trans_score(seq, transcriptome, output_dir):
    matching = alignment_transcriptome(seq, transcriptome, output_dir)

    #score = []
    fragments = list(matching['queryacc.ver'].unique())

    #for fragment in fragments:
    #    sorted_matching = matching[matching['queryacc.ver'] == fragment].sort_values(
    #        by=['identity', 'bitscore', 'alignment length'], ascending=False).sort_values(by=['evalue', 'mismatches'])
    #    score_f = pd.DataFrame(
    #        sorted_matching['bitscore'].astype(float) / sorted_matching['alignment length'].astype(float) *
    #        sorted_matching['identity'].astype(float), columns=['score'])
    #    s_index = score_f.sort_values(by=['score'], ascending=False).iloc[0:1].index.tolist()[0]
    #    score.append([sorted_matching.loc[s_index]['subjectacc.ver'],
    #                  score_f.sort_values(by=['score'], ascending=False).max().max()])
    #score = pd.DataFrame(score, columns=['trans_sseqid', 'score'])
    #score.index = fragments

    score = pd.DataFrame()
    for fragment in fragments:
        f_matching = matching[matching['queryacc.ver'] == fragment]
        score_f = pd.DataFrame([fragment, np.mean(
            f_matching['bitscore'].astype(float) / f_matching['alignment length'].astype(float) * f_matching[
                'identity'].astype(float))]).T
        score = pd.concat([score, score_f])
    score.columns = ['fragments', 'score']
    score.index = fragments

    return (score)

def split_msf(msf_file, output_dir):
    msf_file = open(msf_file, 'r')
    global msf
    out_name = ''
    names = []
    for line in msf_file.readlines():
        line = line.strip()
        if '>' in line:
            name = line[1:].split(':')
            if out_name != '' and out_name not in name:
                msf.close()
            if not os.path.exists(output_dir + '/msf'):
                os.makedirs(output_dir + '/msf')
            out_name = name[0] + '.msf'
            names.append(name[0])
            msf = open(output_dir + '/msf/' + out_name, 'w')
            msf.write(line + '\n')
        else:
            msf.write(line + '\n')

    msf.close()

    return(names)


def run_aligment(protein, predicted_protein, output_dir):
    call_kalign(predicted_protein, output_dir)
    names = split_msf('./kalign/out.msf',output_dir)


    if not os.path.exists(output_dir + '/hmm'):
        os.makedirs(output_dir + '/hmm')
    if not os.path.exists(output_dir + '/out'):
        os.makedirs(output_dir + '/out')
    call_hmmbuild(names, output_dir)
    call_hmmsearch(protein, names, output_dir)





