###############################################################################
#
# Title: CRC and Covariate Generation
# Author: Gabe Foster
# Date: 6/11/2020
# Purpose: This script takes an expression file (BED format) and a covariate
# file, generates cyclical regression covariates, and adds them to the
# covariate file. This is preparation for running PEER.
#
###############################################################################

library(PFExpTools)
library(argparser)

setwd('/data')

p <- arg_parser('Prepare Covariates for PEER')
p <- add_argument(p, 'exp_file', help='Expression file to be staged, GxN named matrix')
p <- add_argument(p, 'cov_file', help = 'Covariate file, csv')
p <- add_argument(p, 'refcourse', help = 'reference time course, space delim')
args <- parse_args(p)
readGct <- function(exp_file){
  read.delim(exp_file,
             sep = '\t',
             #header = 3,
             #skip = 2,
			       header = 1,
             row.names = 'Geneid',
             stringsAsFactors = F)
}

prepareCovs = function(expdata, covfile, refcourse){
  
  stages = stagingByTranscription(expdata, refcourse, heatmap = F)
  crcs = CRCgeneration(stages)[,-1]
  cov = read.csv(covfile, row.names = 1, stringsAsFactors = F)
  final = merge(cov, crcs, by.x = 'row.names', by.y = 'row.names')
  rownames(final) = final$SampleID
  final$SampleID = NULL
  return(final)
}

prepareRefcourse = function(refcourse, times, newtimes){
  
  newcourse = splineCourse(refcourse, times, newtimes)
  return(newcourse)
  
}

###############################################################################
#
# Functional block- this piece does the work using the above functions
#
###############################################################################

# Read exp data file, drop unnecessary cols
expdata <- readGct(args$exp_file)


#################################################
#
# This is specific for a single ref time course;
# it needs to be edited for other reference
# courses. Sorry.
#
#################################################


refcourse = read.delim(args$refcourse,
                       row.names = 1,
                       sep = ' ',
                       stringsAsFactors = F)

refcourse = refcourse[,-c(1:3)]
fullcourse = prepareRefcourse(refcourse, seq(4,56,4), 1:48)

# Generate CRCs, and write it out in a PEER compatible format.

allcovs = prepareCovs(expdata, args$cov_file, fullcourse)
rownames(allcovs) = allcovs$Row.names
allcovs$Row.names = NULL

write.table(allcovs,
            paste(strsplit(args$cov_file, '\\.')[[1]][1],'_crc.txt', sep = ''),
            row.names = TRUE,
            col.names = TRUE,
            sep = '\t',
            quote = FALSE)





