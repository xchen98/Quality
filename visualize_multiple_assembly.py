import os
import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import warnings

def read_score_table(path):
    score = pd.read_table(path)
    assembly_name = score.columns[0].split(' ')[4].split('_')[0] + '_' + score.columns[0].split(' ')[4].split('_')[1]
    table = pd.DataFrame(score[12:])
    score_table = [table.loc[x,].values[0].split()[1:-1] for x in table.index]
    score_table = pd.DataFrame(score_table)
    col = score.loc[11,].values[0].split()

    col[col.index('transcriptome_score')] = "transcriptome"
    col[col.index('protein_score')] = "proteome"
    score_table.columns = col[:-1]
    score_table['assembly'] = [assembly_name] * len(score_table['fragments'])
    return(score_table)

def plot_assembly(score_table, output_path):
    palette = sns.color_palette("flare", len(score_tables))

    sns.relplot(
        data=score_table, x="fragment", y="score",
        col="variable", hue="assembly",
        kind="line", palette=palette,
        height=3, aspect=1, facet_kws=dict(sharex=False)
    ).set_titles("{col_name}")
    sns.despine(offset=10, trim=True)

    plt.savefig(output_path + '/multiple_assemblies.png', dpi=600, bbox_inches='tight')

def plot_multi_assemblies(score_table, output_path):

    sns.set(font_scale=0.3)
    sns.set_theme(style="white")
    sns.relplot(
        data=score_table, x="fragment", y="score",
        col="variable", row="assembly",
        kind="line", hue="assembly",
        height=3, aspect=1, facet_kws=dict(sharex=False),
    ).set_titles("{col_name}")


    plt.savefig(output_path + '/assemblies.png', dpi=600, bbox_inches='tight')

def run(output_path, score_tables):
    warnings.simplefilter("ignore", category=DeprecationWarning)

    mix_table = pd.DataFrame()

    for assembly_i in (score_tables):

        mix_table = pd.concat([mix_table, read_score_table(assembly_i)])

    mix_table = mix_table.fillna(0)
    mix_table.iloc[:, 1:-1] = mix_table.iloc[:, 1:-1].astype(float)
    if 'coverage' in mix_table.columns:
        max_score = max(mix_table.completeness) + max(mix_table.contamination) + max(mix_table.coverage) + max(mix_table.transcriptome) + max(
            mix_table.proteome)
    else:
        max_score = max(mix_table.completeness) + max(mix_table.contamination) + max(mix_table.transcriptome) + max(mix_table.proteome)

    mix_table['score'] = [mix_table.iloc[i, 1:-1].sum() / max_score for i in range(len(mix_table.fragments))]
    mix_table['idx'] = mix_table.index.astype(str)

    mix_table['completeness'] = mix_table['completeness'] / max(mix_table['completeness'])
    mix_table['contamination'] = mix_table['contamination'] / max(mix_table['contamination'])
    mix_table['transcriptome'] = mix_table['transcriptome'] / max(mix_table['transcriptome'])
    mix_table['proteome'] = mix_table['proteome'] / max(mix_table['proteome'])
    if 'coverage' in mix_table.columns:
        mix_table['coverage'] = mix_table['coverage'] / max(mix_table['coverage'])

    mix_melt = pd.melt(mix_table, id_vars=['fragments', 'idx', 'assembly'])
    mix_melt.idx = mix_melt.idx.astype(int) + 1
    mix_melt.columns = ['fragments_accession', 'fragment', 'assembly', 'variable', 'score']

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    plot_assembly(mix_melt, output_path)
    plot_multi_assemblies(mix_melt, output_path)


if __name__ == '__main__':
    path = sys.argv[1]
    score_tables = sys.argv[2:]
    run(path, score_tables)
