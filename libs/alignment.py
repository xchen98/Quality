import os
import numpy as np
import pandas as pd

from libs.call_subprocess import call_makeblastdb, call_subprocess, call_diamond


def check_alignment(output_dir, anno_dir, mode):
    call_subprocess(['Rscript ./Granges.R', '-f', anno_dir, '-s', output_dir, '-t', mode, ' > ', output_dir + '/log.txt'])

def alignment_transcriptome(sequences, transcriptome, output_dir):
    transcriptome_dir = call_makeblastdb(transcriptome, 'nucl', output_dir)
    command = ['blastn', '-query', sequences, '-db', transcriptome_dir,
               '-evalue 10e-6 -num_threads 8 -outfmt 7']

    out = call_subprocess(command, out=True, align=True)
    return(out)


def get_protein_score(seq, protein, fragments, output_dir, annotation = None):
    call_diamond(seq, protein, output_dir)
    matching = pd.read_csv(output_dir + '/matches.m8', sep="\t")
    matching.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                        'send', 'evalue', 'bitscore']
    matching.to_csv(output_dir + '/prot.tsv', sep = '\t')

    if annotation:
        check_alignment(output_dir, annotation, 'prot')
        prot_out = pd.read_table(output_dir + '/prot_out.tsv', sep = '\t').x.tolist()
        matching['bitscore'] = matching['bitscore'] * prot_out

    score = pd.DataFrame()
    for fragment in fragments:
        f_matching = matching[matching['qseqid'] == fragment]
        score_f = pd.DataFrame([fragment, np.mean(
            f_matching['bitscore'].astype(float) / f_matching['length'].astype(float) * f_matching[
                'pident'].astype(float))]).T
        score = pd.concat([score, score_f])
    score.columns = ['fragments', 'score']
    score.index = fragments

    os.remove(output_dir + '/prot.tsv')
    os.remove(output_dir + '/prot_out.tsv')

    return score

def get_trans_score(seq, transcriptome, fragments, output_dir, annotation = None):
    alignment_transcriptome(seq, transcriptome, output_dir).to_csv(output_dir + '/trans.tsv', sep = '\t')

    matching = pd.read_table(output_dir + '/trans.tsv', sep = '\t')

    if annotation:
        check_alignment(output_dir, annotation, 'trans')
        trans_out = pd.read_table(output_dir + '/trans_out.tsv', sep = '\t').x.tolist()
        matching['bitscore'] = matching['bitscore'] * trans_out


    score = pd.DataFrame()
    for fragment in fragments:
        f_matching = matching[matching['queryacc.ver'] == fragment]
        score_f = pd.DataFrame([fragment, np.mean(
            f_matching['bitscore'].astype(float) / f_matching['alignment length'].astype(float) * f_matching[
                'identity'].astype(float))]).T
        score = pd.concat([score, score_f])
    score.columns = ['fragments', 'score']
    score.index = fragments

    os.remove(output_dir + '/trans.tsv')
    os.remove(output_dir + '/trans_out.tsv')

    return score

