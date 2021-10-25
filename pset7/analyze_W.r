# Wiggins' RNA-seq analysis script
#
# To run:
#    % Rscript analyze_W.r
#
# You need to have your data in a file called "mydata.tbl"
# This is a tab-delimited file with a header line, followed
# by one line per gene: 
#   <genename>  <counts1> <counts2> ... <counts6>
# Script assumes that there are six samples, three from one
# condition (e.g. wild type), three from another (e.g. mutant).
#
# The script generates an output file "myresult.out".
# 

library(edgeR)
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

group  <- factor(c(1,1,1,2,2,2))     

x     <- read.table(args[1], sep='\t', row.names=1)
y     <- DGEList(counts=x,group=group)
keep  <- filterByExpr(y)        ## The filterByExpr function keeps rows that have worthwhile counts in 
                                ## a minumum number of samples. This one's probably unnecessary but 
                                ## let's add it just in case.
y     <- y[keep, , keep.lib.sizes=FALSE]
y     <- calcNormFactors(y)     ## normalization!
y     <- estimateDisp(y)
et    <- exactTest(y)
tab   <- topTags(et, nrow(x))

write.table(tab, file=args[2])