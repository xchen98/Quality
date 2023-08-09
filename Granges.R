library(GenomicFeatures)
library(rtracklayer)
library(getopt)

spec <- matrix(
  c("first",  "f", 1, "character", "This is the path for gtf file!",
    "second", "s", 1, "character",  "This is the path for the trans/pro file!",
    "third", "t", 1, "character", "Mode: trans/prot"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)


makeGrange = function(query, q_start, q_end){

  if (all(q_end > q_start) == FALSE){
    new_end = c()
    new_start = c()
    com = (q_end>q_start)
    for(i in 1:length(q_start)){
      if (com[i] == FALSE){
        new_end = append(new_end, q_start[i])
        new_start = append(new_start, q_end[i])
      }else{
        new_end = append(new_end, q_end[i])
        new_start = append(new_start, q_start[i])
      }
    }
  }else{
    new_end = q_end
    new_start = q_start
  }
  grange = GRanges(seqnames = query, ranges = IRanges(start = new_start, end = new_end))

  return(grange)
}

gtf <- read.table(opt$first, sep = '\t')
granges <- GRanges(seqnames = gtf$V1, ranges = IRanges(start = gtf$V4, end = gtf$V5), type = gtf$V3)

query <- read.table(paste0(opt$second, '/', opt$third, '.tsv'), sep = '\t', header = 1, row.names = 1)

if (opt$third == 'trans'){
    query_grange = makeGrange(query$queryacc.ver, query$q.start, query$q.end)
}else{
    query_grange = makeGrange(query$qseqid, query$qstart, query$qend)
}

overlap <- findOverlaps(query_grange, granges)


mis = (1:length(query_grange))[!(1:length(query_grange) %in% unique(overlap@from))]
factor = c()
last = 0
for(i in mis){
    factor = append(factor, rep(1, i-last-1))
    factor = append(factor, 0)
    last = i
}
factor = append(factor, rep(1, length(query_grange)-last))

write.table(factor, file = paste0(opt$second, '/', opt$third, '_out.tsv'), sep = '\t')
write.table(query_grange[mis], file = paste0(opt$second, '/', opt$third, '_anno_mistake.txt'))
