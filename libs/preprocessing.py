
def load_fragments(path):
    fragments = path + '/fragments.txt'
    fragment_names = []
    start_positions = []
    end_positions = []
    with open(fragments, 'r') as fragments:
        for line in fragments.readlines():
            line = (line.strip('\n')).split('\t')
            fragment_names.append(line[0])
            start_positions.append(line[7])
            end_positions.append(line[8])


