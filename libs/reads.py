import os
from libs.call_subprocess import call_subprocess

def search_Biosample(accession):
    query = '"' + str(accession) + '"'
    command = ['esearch -db assembly -query', query,
               '| esummary | xtract -pattern DocumentSummary -element BioSampleAccn']
    Biosample = call_subprocess(command, out=True)
    return(Biosample)

def search_readsacc(Biosample):
    query = Biosample[0]
    command = ['esearch -db biosample -query', query,
               '| elink -target sra | esummary | xtract -pattern Runs -first Run@acc']
    readsacc = call_subprocess(command, out=True)
    return(readsacc)

def read_run(accession, outputdir):
    Biosample = search_Biosample(accession)

    readsacc = search_readsacc(Biosample)
    if len(readsacc) == 0:
        print('None reads are found!')
        return None
    else:
        print('Done! We find the corresponding reads.', readsacc)
        print('Start to download reads.')
        if not os.path.exists(outputdir + '/reads'):
            os.makedirs(outputdir + '/reads')
        call_subprocess(['prefetch', "--max-size 10G", readsacc[0]])
        call_subprocess(['fasterq-dump --concatenate-reads', readsacc[0], '-O', outputdir + '/reads'])



        return(readsacc)

