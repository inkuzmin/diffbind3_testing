## Setup R error handling to go to stderr
options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr()); q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
Sys.setlocale("LC_TIME", "en_US.UTF-8")
Sys.setlocale("LC_MONETARY", "en_US.UTF-8")
Sys.setlocale("LC_PAPER", "en_US.UTF-8")
Sys.setlocale("LC_MEASUREMENT", "en_US.UTF-8")
Sys.setlocale("LC_NUMERIC", "C")
Sys.setlocale("LC_COLLATE", "en_US.UTF-8")
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages(library("DiffBind"));
suppressPackageStartupMessages(library("getopt"));
suppressPackageStartupMessages(library("jsonlite"));

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)


is.formula <- function(x) {
    inherits(x, "formula")
}

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    "infile", "i", 1, "character",
    "scorecol", "n", 1, "integer",
    "lowerbetter", "l", 1, "logical",
    "summits", "s", 1, "integer",
    "th", "t", 1, "double",
    "plots", "p", 2, "character",
    "bmatrix", "b", 0, "logical",
    "rdaOpt", "r", 0, "logical",
    "infoOpt", "a", 0, "logical",

    "filter", NA, 1, "integer",
    "mincount", NA, 1, "integer",
    "design", NA, 1, "logical",
    "formula", NA, 1, "character",

    "custompeaks", NA, 1, "character",

    "verbose", "v", 2, "integer",
    "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4);

opt <- getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE));
    q(status = 1);
}

sample <- dba(sampleSheet = opt$infile,
              scoreCol = opt$scorecol,
              bLowerScoreBetter = opt$lowerbetter)

sample$config$RunParallel <- FALSE

if ( !is.null(opt$summits) ) {
    if (is.null(opt$custompeaks)) {
       sample_count <- dba.count(sample, summits = opt$summits, filter = opt$filter, minCount = opt$mincount)
    } else {
       sample_count <- dba.count(sample, peaks=opt$custompeaks, summits = TRUE, filter = opt$filter, minCount = opt$mincount)
    }
} else {
    if (is.null(opt$custompeaks)) {
        sample_count <- dba.count(sample, filter = opt$filter, minCount = opt$mincount)
    } else {
        sample_count <- dba.count(sample, peaks=opt$custompeaks, summits = TRUE, filter = opt$filter, minCount = opt$mincount)
    }
}

if ( opt$design ) {
    sample$config$design <- TRUE

    if ( !is.formula(opt$formula) ) {
        sample_contrast <- dba.contrast(sample_count, design=TRUE)
    } else {
        sample_contrast <- dba.contrast(sample_count, design=opt$design)
    }
} else {
    sample$config$design <- FALSE
    sample_contrast <- dba.contrast(sample_count, categories = DBA_CONDITION, minMembers = 2)
}

sample_analyze <- dba.analyze(sample_contrast)

# Generate plots
if (!is.null(opt$plots)) {
    pdf(opt$plots)

    show <- dba.show(sample_analyze, bContrasts = TRUE, th = opt$th)
    n_rows <- nrow(show)

    for (i in 1:n_rows) {
        sites <- show[i, "DB.DESeq2"]
        factor <- show[i, "Factor"]
        group <- show[i, "Group"]
        group2 <- show[i, "Group2"]


        if (sites > 0) {
            plot.new()
            title(paste(factor, ":", group, group2))


            tryCatch(
                expr = {
                    dba.plotHeatmap(sample_analyze, contrast = i, correlations = FALSE, cexCol = 0.8, th = opt$th)
                },
                error=function(cond){
                    message(cond)
                    return(NA)
                }
            )

            tryCatch(
                expr = {
                    dba.plotPCA(sample_analyze, contrast = i, th = opt$th, label = DBA_ID, labelSize = 0.3)
                },
                error=function(cond){
                    message(cond)
                    return(NA)
                }
            )
            
            tryCatch(
                expr = {
                    dba.plotMA(sample_analyze, contrast = i, th = opt$th)
                },
                error=function(cond){
                    message(cond)
                    return(NA)
                }
            )
            
            tryCatch(
                expr = {
                    dba.plotVolcano(sample_analyze,  contrast = i, th = opt$th)
                },
                error=function(cond){
                    message(cond)
                    return(NA)
                }
            )
            
            tryCatch(
                expr = {
                    dba.plotBox(sample_analyze, contrast = i, th = opt$th)
                },
                error=function(cond){
                    message(cond)
                    return(NA)
                }
            )

            # dba.plotVenn(sample_analyze, contrast = i, th = opt$th)

        }
    }

    dev.off()
}

# Output reports
    dir.create("reports")
    original_dir <- getwd()
    setwd("reports")

    show <- dba.show(sample_analyze, bContrasts = TRUE, th = opt$th)
    n_rows <- nrow(show)

    for (i in 1:n_rows) {
        sites <- show[i, "DB.DESeq2"]
        factor <- show[i, "Factor"]
        group <- show[i, "Group"]
        group2 <- show[i, "Group2"]

        if (sites > 0) {
            dba.report(sample_analyze, th = opt$th, contrast = i, file = paste(factor, group, group2, sep="__"))
        }
    }

    setwd(original_dir)



# Output binding affinity scores
if (!is.null(opt$bmatrix)) {
    bmat <- dba.peakset(sample_count, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
    # Output as 0-based tabular
    bmat <- data.frame(Chrom = bmat[, 1],
        Start = bmat[, 2] - 1,
        End = bmat[, 3],
        bmat[, 4:ncol(bmat)])
    write.table(bmat, file = "bmatrix.tab", sep = "\t", quote = FALSE, row.names = FALSE)
}

# Output RData file
if (!is.null(opt$rdaOpt)) {
    save.image(file = "DiffBind_analysis.RData")
}

# Output analysis info
if (!is.null(opt$infoOpt)) {
    info <- "DiffBind_analysis_info.txt"
    cat("dba.count Info\n\n", file = info, append = TRUE)
    capture.output(sample, file = info, append = TRUE)
    cat("\ndba.analyze Info\n\n", file = info, append = TRUE)
    capture.output(sample_analyze, file = info, append = TRUE)
    cat("\nSessionInfo\n\n", file = info, append = TRUE)
    capture.output(sessionInfo(), file = info, append = TRUE)
}
