from libs.call_subprocess import call_subprocess

class Contamination:
    def __init__(self, database, query):
        self.database = database
        self.query = query
        self.score = []


    def contamination_score(self, scaffold_l, cleaned_scaffold_l):
        sub = list(map(lambda x: x[0] - x[1], zip(scaffold_l, cleaned_scaffold_l)))
        #因为busco的completeness是越大越好
        sum_l = sum(scaffold_l)
        self.score = [round(1-x/sum_l, 2) for x in sub]

        return

    def Pipline(self, type, query, database):
        global blastn, command
        if type == "vector" or "adaptor":
            command = ['blastn', '-query', query, '-db', database, '-reward 1 -penalty -5 -gapopen 3',
                       '-gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt',
                       '"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"',
                       "| awk '$3>=94.0'"]
        if type == "common_cont":
            command = ['blastn', '-query', query, '-db', database, '-task megablast', '-word_size 28',
                       '-evalue 1e-3', '-perc_identity 94', '-outfmt', '"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"']

        blastn = call_subprocess(command, out=True, contam=True, mode=type)
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

    def merge_contamination(self, cont_1, cont_2):
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
            for i in range(l):
                if i < low_l:
                    if low[0][i] == large[0][i]:
                        rc_areas = (low[1][i] + large[1][i])
                        if len(self.remove_duplicated(rc_areas)) != 0:
                            all_cont_area.append(self.remove_duplicated(rc_areas))
                            all_cont_name.append(large[0][i])

                    else:
                        if len(self.remove_duplicated(low[1][i])) != 0:
                            all_cont_name.append(low[0][i])
                            all_cont_area.append(self.remove_duplicated(low[1][i]))
                        if len(self.remove_duplicated(large[1][i])) != 0:
                            all_cont_name.append(large[0][i])
                            all_cont_area.append(self.remove_duplicated(large[1][i]))

                else:
                    if len(self.remove_duplicated(large[1][i])) != 0:
                        all_cont_name.append(large[0][i])
                        all_cont_area.append(self.remove_duplicated(large[1][i]))

            all_cont.append(all_cont_name)
            all_cont.append(all_cont_area)

        return (all_cont)

    def run(self):
        #vec = self.Pipline("vector", self.query, self.database + '/UniVecDB/UniVecDB')
        adaptor = self.Pipline("adaptor", self.query, self.database + '/Adaptors/adaptors_for_screening_euks')
        common_cont = self.Pipline("common_cont", self.query, self.database + '/contam_in_euks/contam_in_euks')
        rrna = self.Pipline("rrna", self.query, self.database + '/rrna/rrna')

        all_contamination = [[], []]
        if len(adaptor[0]) != 0:
            all_contamination = self.merge_contamination(all_contamination, adaptor)
        if len(common_cont[0]) != 0:
            all_contamination = self.merge_contamination(all_contamination, common_cont)
        if len(rrna[0]) != 0:
            all_contamination = self.merge_contamination(all_contamination, rrna)

        return all_contamination
