from libs.call_subprocess import call_subprocess

class Contamination:
    def __init__(self, database, query):
        self.database = database
        self.query = query
        self.score = []


    def contamination_score(self, scaffold_l, cleaned_scaffold_l):
        sub = list(map(lambda x: x[0] - x[1], zip(scaffold_l, cleaned_scaffold_l)))
        sum_l = sum(scaffold_l)
        self.score = [round((1-x/sum_l)*100, 2) for x in sub]

        return

    def Pipline(self, type, query, database):
        global blastn, command

        if type == "mito":
            command = ['blastn', '-query', query, '-db', database, '-task megablast',
                       "-word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1",
                       "-dust yes -evalue 1e-5 -perc_identity 98.6",
                       "-gapextend 2 -gapopen 4 -penalty -2 -reward 1 -num_threads 8",
                       '-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                       "| awk '$4>=150'"]

        if type == "yeast":
            command = ['blastn', '-query', query, '-db', database, '-task megablast',
                       "-word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1",
                       "-dust yes -evalue 1e-5 -perc_identity 95",
                       "-gapextend 2 -gapopen 4 -penalty -2 -reward 1 -num_threads 8",
                       '-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                       "| awk '$4>=100'"]

        if type == "Ecoli":
            command = ['blastn', '-query', query, '-db', database, '-task megablast',
                       "-word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1",
                       "-dust yes -evalue 1e-5 -perc_identity 95",
                       "-gapextend 2 -gapopen 4 -penalty -2 -reward 1 -num_threads 8",
                       '-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                       "| awk '$4>=100'"]

        blastn = call_subprocess(command, out=True, contam=True)
        return blastn

    def merge_area(self, area_1, area_2):
        merged_area = []
        if area_2[0] <= area_1[1]:
            if area_2[1] < area_1[1]:
                merged_area.append(area_1)
            else:
                new_area = [area_1[0], area_2[1]]
                merged_area.append(new_area)
        else:
            merged_area.append(area_1)
            merged_area.append(area_2)

        return merged_area

    def remove_duplicated(self, areas):
        cont_area = []
        for area in areas:
            if area not in cont_area:
                if len(cont_area) != 0:

                    last_area = cont_area[-1]
                    if area[0] <= last_area[1]:
                        cont_area = cont_area[:-1]
                        for area_i in self.merge_area(last_area, area):
                            cont_area.append(area_i)
                    else:
                        cont_area.append(area)
                else:
                    cont_area.append(area)
        cont_area.sort(key=lambda x: (x.sort(), x[0], x[1]))

        return cont_area

    def merge_contamination(self, cont_1, cont_2, fragments_name):
        all_cont = []
        all_cont_name = []
        all_cont_area = []
        l1_l = len(cont_1[0])
        l2_l = len(cont_2[0])
        if l1_l == 0:
            all_cont = cont_2
        elif l2_l == 0:
            all_cont = cont_1
        else:
            l = max(l1_l, l2_l)
            if l1_l >= l2_l:
                low_l = l2_l
                large = cont_1
                low = cont_2
            else:
                low_l = l1_l
                large = cont_2
                low = cont_1
            index_l = [large[0].index(x) for x in low[0] if x in large[0]]
            low_done = False
            low_j = 0
            for i in range(l):
                if low_done == False:
                    if i in index_l:
                        rc_areas = (low[1][index_l.index(i)] + large[1][i])
                        if len(self.remove_duplicated(rc_areas)) != 0:
                            all_cont_area.append(self.remove_duplicated(rc_areas))
                            all_cont_name.append(large[0][i])
                            low_j = low_j + 1


                    else:
                        if i < low_l and low[0][i] not in large[0][i:]:
                            if len(self.remove_duplicated(low[1][i])) != 0:
                                ind = fragments_name.index(low[0][i])
                                all_cont_name.insert(ind, low[0][i])
                                all_cont_area.insert(ind, (self.remove_duplicated(low[1][i])))
                                low_j = low_j + 1
                        elif i >= low_l and low_j <= low_l:
                            print( i, low_l,  low_j, low)
                            if len(self.remove_duplicated(low[1][low_j])) != 0:
                                ind = fragments_name.index(low[0][low_j])
                                all_cont_name.insert(ind, (low[0][low_j]))
                                all_cont_area.insert(ind, (self.remove_duplicated(low[1][low_j])))
                                low_j = low_j + 1

                        if len(self.remove_duplicated(large[1][i])) != 0:
                            ind = fragments_name.index(large[0][i])
                            all_cont_name.insert(ind, (large[0][i]))
                            all_cont_area.insert(ind, (self.remove_duplicated(large[1][i])))

                    if set(low[0]) < set(all_cont_name): low_done = True

                else:
                    if len(self.remove_duplicated(large[1][i])) != 0:
                        ind = fragments_name.index(large[0][i])
                        all_cont_name.insert(ind, (large[0][i]))
                        all_cont_area.insert(ind, (self.remove_duplicated(large[1][i])))

            all_cont.append(all_cont_name)
            all_cont.append(all_cont_area)

        return (all_cont)

    def run(self, fragments_name):
        mito = self.Pipline("mito", self.query, self.database + '/mito/mito')
        Ecoli = self.Pipline("Ecoli", self.query, self.database + '/Ecoli/Ecoli')
        yeast = self.Pipline("yeast", self.query, self.database + '/yeast/yeast')

        all_contamination = [[], []]

        if len(mito[0]) != 0:
            all_contamination = self.merge_contamination(all_contamination, mito, fragments_name)
        if len(Ecoli[0]) != 0:
            all_contamination = self.merge_contamination(all_contamination, Ecoli, fragments_name)
        if len(yeast[0]) != 0:
            all_contamination = self.merge_contamination(all_contamination, yeast, fragments_name)

        return all_contamination
