import copy
from libs.call_subprocess import call_subprocess




def scaffold(seq, trim_name, trim_area, output_path):
    scaffolds_l = []
    cleaned_scaffolds_l = []
    contigs_l = []
    trim = False
    index = None
    area = None
    n = 0
    seq_name = None
    s = ""
    sqs = []
    cleaned_file = output_path + '/cleaned_sequence.fna'
    cleaned_sequence = open(cleaned_file, "w")
    with open(seq) as file:
        scaffold_length = 0
        cleaned_scaffold_length = 0
        contig = 0
        for line in file:
            line = line.strip()
            if len(line) != 0:
                if '>' in line:
                    print(line)
                    scaffolds_l.append(scaffold_length)
                    sqs.append(s)
                    cleaned_scaffolds_l.append(cleaned_scaffold_length)
                    contigs_l.append(contig)
                    scaffold_length = 0
                    s = ""
                    cleaned_scaffold_length = 0
                    contig = 0
                    name = line.replace(' ', '').split('>')[1]
                    if name in trim_name :
                        trim = True
                        index = trim_name.index(name)
                        if len(trim_area[index]) != 0:
                            area = trim_area[index][0]
                            n = 0
                            seq_name = line
                        else:
                            trim = False
                        cleaned_sequence.write('_'.join([line, str(n)]) + '\n')

                    else:
                        trim = False
                        cleaned_sequence.write(line + '\n')

                else:
                    if trim == True:
                        o_scaffold_length = scaffold_length
                        s = s + line
                        scaffold_length = scaffold_length + len(line)
                        if area[0] == 1: area[0] = 0
                        if area[0] <= o_scaffold_length and area[1] <= scaffold_length and area[1] >= o_scaffold_length:
                            if len(trim_area[index]) > 1:
                                last_end = area[1]
                                trim_area[index] = trim_area[index][1:]
                                area = trim_area[index][0]
                                n = n + 1
                                new_name = '_'.join([seq_name, str(n)])
                                cleaned_sequence.write(new_name + '\n')
                            if len(trim_area[index]) == 1 and area[0] < len(line) and area[1] <= scaffold_length:
                                seq = line[0:area[0]] + line[area[1] - scaffold_length:]
                                cleaned_scaffold_length = cleaned_scaffold_length + len(seq)
                                cleaned_sequence.write(seq)
                            elif area[0] < scaffold_length and area[1] <= scaffold_length:
                                seq = line[o_scaffold_length - last_end - 1:o_scaffold_length - area[0] - 1] + line[area[1] - scaffold_length:]
                                cleaned_scaffold_length = cleaned_scaffold_length + len(seq)
                                if (area[0] - last_end - 1) + (scaffold_length - area[1]) < len(line):
                                    cleaned_sequence.write(seq)
                                else:
                                    cleaned_sequence.write(seq + '\n')

                            elif area[0] <= scaffold_length:
                                seq = line[last_end - o_scaffold_length - 1:area[0] - o_scaffold_length - 1]
                                cleaned_scaffold_length = cleaned_scaffold_length + len(seq)
                                if (area[0] - last_end - 1) < len(line):
                                    cleaned_sequence.write(seq)
                                else:
                                    cleaned_sequence.write(seq + '\n')

                            else:
                                seq = line[area[1] - scaffold_length:]
                                cleaned_scaffold_length = cleaned_scaffold_length + len(seq)
                                if (scaffold_length - area[1]) < len(line):
                                    cleaned_sequence.write(seq)
                                else:
                                    cleaned_sequence.write(seq + '\n')

                        elif o_scaffold_length < area[0] and area[1] <= scaffold_length:
                            seq = line[0:area[0] - o_scaffold_length - 1] + line[area[1] - scaffold_length:]
                            cleaned_sequence.write(seq + '\n')
                            cleaned_scaffold_length = cleaned_scaffold_length + len(seq)

                        elif o_scaffold_length < area[0] and area[1] > scaffold_length:
                            seq = line[0:area[0] - o_scaffold_length - 1]
                            cleaned_sequence.write(seq + '\n')
                            cleaned_scaffold_length = cleaned_scaffold_length + len(seq)
                        elif scaffold_length > area[1]:
                            seq = line
                            cleaned_sequence.write(seq + '\n')
                            cleaned_scaffold_length = cleaned_scaffold_length + len(seq)

                    else:
                        scaffold_length = scaffold_length + len(line)
                        s = s + line
                        cleaned_sequence.write(line + '\n')
                        cleaned_scaffold_length = scaffold_length

                    if 'N' in line:
                        splits = line.split('N')
                        for i in range(len(splits)):
                            if i == 0 and splits[i] == '' and contig != 0:
                                contigs_l.append(contig)
                                contig = 0
                            if splits[i] != '':
                                contig = contig + len(splits[i])
                                if i < len(splits) - 1:
                                    contigs_l.append(contig)
                                    contig = 0

                    else:
                        contig = contig + len(line)

        scaffolds_l.append(scaffold_length)
        sqs.append(s)
        cleaned_scaffolds_l.append(cleaned_scaffold_length)
        contigs_l.append(contig)

        scaffolds_l = scaffolds_l[1:]
        sqs = sqs[1:]
        cleaned_scaffolds_l = cleaned_scaffolds_l[1:]
        contigs_l = contigs_l[1:]

    file.close()
    cleaned_sequence.close()

    return (scaffolds_l, cleaned_scaffolds_l, contigs_l, sqs)

def find_X(list, X, Genome_size = None):
    if Genome_size == None:
        threshold = sum(list)*X/100
        add = 0
        for i in range(len(list)):
            add = add + list[i]
            if add >= threshold:
                return(list[i], i)


def NGX_X(scaffolds_l, contigs_l, X, Genome_size = None):
    stat = dict()

    scaffold_NX, scaffold_LX = find_X(scaffolds_l, X)
    contig_NX, contig_LX = find_X(contigs_l, X)

    stat['scaffold N' + str(X)] = scaffold_NX
    stat['scaffold L' + str(X)] = scaffold_LX + 1
    stat['contig N' + str(X)] = contig_NX
    stat['contig L' + str(X)] = contig_LX + 1

    if Genome_size:
        scaffold_NGX, scaffold_LGX = find_X(scaffolds_l, X, Genome_size)
        contig_NGX, contig_LGX = find_X(contigs_l, X, Genome_size)
        stat['scaffold NG' + str(X)] = scaffold_NGX
        stat['scaffold LG' + str(X)] = scaffold_LGX + 1
        stat['contig NG' + str(X)] = contig_NGX
        stat['contig LG' + str(X)] = contig_LGX + 1


    return(stat)

def run_statistic(seq, trim_name, trim_area, output_path):
    scaffolds_l, cleaned_scaffolds_l, contigs_l, sqs = scaffold(seq, trim_name, trim_area, output_path)
    new_scaffolds_l = copy.deepcopy(scaffolds_l)
    scaffolds_l.sort(reverse=True)
    contigs_l.sort(reverse=True)



    sum_scaffold = sum(scaffolds_l)
    sum_contig = sum(contigs_l)
    statistics = {}
    statistics['genome size'] = sum_scaffold
    statistics['total ungapped length'] = sum_contig
    statistics['gaps'] = (sum_scaffold - sum_contig)/sum_scaffold
    statistics.update(NGX_X(scaffolds_l, contigs_l, 50))
    statistics['scaffold total'] = len(scaffolds_l)
    statistics['contig total'] = len(contigs_l)


    return (statistics, sqs, new_scaffolds_l, cleaned_scaffolds_l)

