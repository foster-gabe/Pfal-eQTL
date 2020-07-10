###############################################################################
#
# Filter nominals by q-value thresholds
#
###############################################################################


library(argparser)

p <- arg_parser('Filter QTLtools nominal output by q-value')
p <- add_argument(p, 'nominal_file', help='QTLtools nominal output')
p <- add_argument(p, 'threshold_file', help = 'Q-value threshold file')
p <- add_argument(p, 'output_file', help = 'Output file name')
args <- parse_args(p)


nom = read.delim(args$nominal_file, sep = ' ', header = F)

cutoffs = read.delim(args$threshold_file, sep = ' ', header = F)

cutoffslist = cutoffs$V2
names(cutoffslist) = cutoffs$V1

test = nom[nom$V12 < cutoffslist[nom$V1],]

write.table(test, quote = F, row.names = F, col.names = F, sep = '\t')