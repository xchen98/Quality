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
    print('Done! We find the corresponding reads.', readsacc[1:])

    print('Start to download reads.')
    for acc in readsacc[1:]:
        call_subprocess(['prefetch', "--max-size 100G", acc])
        call_subprocess(['fasterq-dump --concatenate-reads', acc, '-O', outputdir])
        read = open(outputdir + '/' + acc + '.fastq', 'r')
        call_subprocess(["echo", read.readlines(), '>>', outputdir + '/reads.fastq'])
        read.close()

    return(readsacc)

