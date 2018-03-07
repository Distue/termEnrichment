#' termEnrichment
#'
#' This function calculates enrichment for a set of features against a set
#' of background features.
#'
#' @param termTable data.frame or tibble ID where first column is an ID column and the second column are associated terms to be tested against
#' @param foregroundIDs vector of IDs defining the foreground (must be same IDs as is first column of termTable)
#' @param backgroundIDs vector of IDs defining the background (must be same IDs as is first column of termTable)
#' @param annotation data.frame or tibble annotating terms (first column must be term identifier)
#' @param terms specifies what terms will be tested for. tests all existing terms if NULL
#' @param permutations number of permutations used for caluclating mean background ranks for random samples data
#' @param correction logical. TRUE if p-value correction is applied
#' @param padj.cutoff cutoff for p.adjusted value
#' @param padj.method p adjustment method, default is "ihw" see ?ihw, other options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" from ?p.adjust
#' @param alternative alternative hypothesis for fisher exact test, see ?fisher.test
#' @param parallel logical the testing be parallised with BiocParallel
#' @param quiet non verbose output
#' @param removeDuplicatedTerms removes duplicated terms in termTable if necessary
#' @param ... parameter list forwarded to fisher.test functions
#' @author Thomas Schwarzl <schwarzl@embl.de>
#' @import tibble
#' @import dplyr
#' @import graphics
#' @import utils
#' @importFrom tibble as.tibble
#' @importFrom graphics barplot
#' @importFrom graphics par
#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
#' @importFrom utils head
#' @return tibble annotated analysis results
#' @export
#' @usage termEnrichment(termTable, foregroundIDs, backgroundIDs, annotation, terms, padj.cutoff, padj.method, alternative, parallel, removeDuplicatedTerms)
#' @examples
#'
#' data(yeastGO)
#' data(backgroundIDs)
#' data(foregroundIDs)
#' data(yeastGOdesc)
#'
#' termEnrichment(yeastGO, foregroundIDs, backgroundIDs,
#'                annotation = yeastGOdesc)

termEnrichment <- function(termTable, foregroundIDs, backgroundIDs = NULL,
                           annotation = NULL, terms = NULL, permutations = 1e2,
                           correction = F, padj.cutoff = 0.05, padj.method = "IHW",
                           alternative = "greater", parallel = F,
                           dropIDs = F,  quiet = F, removeDuplicatedTerms = T,
                           ...) {
    # ----------- INPUT CHECK -------------
    if(! (nrow(termTable) > 0 && ncol(termTable) == 2))
        stop("'termTable' has to have 2 columns, first one being the ID column, second one the term column.")

    if(! (is.null(backgroundIDs) || length(backgroundIDs) > 0))
        stop(paste0("'backgroundIDs' cannot be empty, if you have no background ids, ",
             "set 'backgroundIDs = NULL' (current length: ", length(backgroundIDs), ", head: ", head(backgroundIDs), ")"))

    if(! length(foregroundIDs) > 0)
        stop(paste0("'foregroundIDs' cannot be empty (current length: ", length(foregroundIDs), ")"))

    if(! (is.data.frame(termTable) || is.tibble(termTable)))
        stop(paste0("'termTable' has to be a data.frame or a tibble.", class(termTable), ")"))

    if(! is.logical(dropIDs))
        stop(paste0("'dropIDs' parameter has to be either TRUE or FALSE (logical, currently ", class(dropIDs), ")"))

    # rename termTable column names
    colnames(termTable) <- c("ID", "term")

    if( class(termTable %>% .[["ID"]]) != class(foregroundIDs))
        stop(paste0("ID column of 'termTable' (", class(termTable[,"ID"]), ") ",
                    "has a different class than 'foregroundIDs' (", class(foregroundIDs),")"))

    if(!is.null(backgroundIDs) && class(termTable %>% .[["ID"]]) != class(backgroundIDs))
        stop(paste0("ID column of 'termTable' (", class(termTable[,"ID"]), ") ",
                    "has a different class than 'backgroundIDs' (", class(backgroundIDs),")"))

    # check that no rows are duplicated
    if( (termTable %>% duplicated %>% any)) {
        if( removeDuplicatedTerms ) {
            termTable <- termTable[ ! termTable %>% duplicated,  ]
        } else {
            stop("There cannot be duplicated rows in 'termTable'." )
        }
    }

    # if annotation not null it has to be a data frame or a tibble
    if(!is.null(annotation)) {
        if(! (is.data.frame(annotation) || is.tibble(annotation)))
            stop(paste0("'annotation' has to be either a data.frame or a tibble (currently: ", class(annotation),")"))
        if(! any(colnames(annotation) == "term"))
            stop("'annotation' has to have a column called term")
        if((annotation %>% .[["term"]] %>% sort %>% duplicated %>% any))
            stop("'annotation' cannot have any duplicated columns")
    }

    # if no background list was specified
    if( is.null(backgroundIDs) ) {
        # set all ids from the mapping table as background
        backgroundIDs = termTable %>% .[["ID"]] %>% sort %>% unique
    } else {
        # assert that all background IDs are found  in the mapping table
        if(! backgroundIDs %in% termTable$ID %>% all) {
            # if dropIDs is TRUE, drop all IDs not found in the termTable
            if(dropIDs) {
                # if we would lose all the IDs, raise error
                if( (! backgroundIDs %in% termTable$ID) %>% all )
                    stop("None of the background IDs were in the 'termTable' ID column.")

                cat(paste0(sum(! backgroundIDs %in% termTable$ID), " background IDs were dropped\n"))
                backgroundIDs <- backgroundIDs[ backgroundIDs %in% termTable$ID ]
            } else {
                stop(paste0("not all backgroundIDs are found in the ID columns of termTable",
                            " e.g.(", paste(backgroundIDs[ ! backgroundIDs %in% termTable$ID ] %>% head, collapse = " "),")"))
            }
        }
    }

    # assert that all foreground IDs are found in the mapping table
    if(! foregroundIDs %in% termTable$ID %>% all) {
        # if dropIDs is TRUE, drop all IDs not found in the termTable
        if( dropIDs) {
            cat(paste0(sum(! foregroundIDs %in% termTable$ID), " out of ", length(foregroundIDs),
                       " foreground IDs were dropped\n"))
            foregroundIDs <- foregroundIDs[ foregroundIDs %in% termTable$ID ]
        } else {
            stop(paste0("not all foregroundIDs are found in the ID columns of termTable",
                        " e.g.(", paste(foregroundIDs[ ! foregroundIDs %in% termTable$ID ] %>% head, collapse = " "),")"))
        }
    }

    # assert that all foreground IDs are found in the background list
    if(! foregroundIDs %in% backgroundIDs %>% all)
        stop(paste0("not all foregroundIDs are found in backgroundIDs",
                    "(", paste(foregroundIDs[! foregroundIDs %in% backgroundIDs] %>% head, collapse = " ") ,")"))

    # ----------- SETUP -------------
    # if terms is NULL all terms will be tested
    if(is.null(terms)) {
        # get a list of unique terms
        terms <- termTable %>% .[["term"]] %>% sort %>% unique
        terms <- terms [ terms != "" ]
    } else {
        if(! class(terms) != class(termTable[,"term"]))
            stop(paste0("'terms' variable has a different class (", class(terms) ,") ",
                  "from the column in the term table (", class(termTable[,"term"]), ")"))
        if(! terms %in% termTable$term %>% all)
            stop(paste0("'terms' variable has values not occuring",
                          " in the 'termTable' (e.g. ",  head(terms[! terms %in% termTable$term]) ,")"))
    }

    # generate a subset for background and foreground
    BACKGROUND <- termTable %>% dplyr::filter(ID %in% backgroundIDs)
    FOREGROUND <- termTable %>% dplyr::filter(ID %in% foregroundIDs)

    # parallelize
    if ( parallel ) {
        require( BiocParallel )
        LAPPLY = bplapply
    } else {
        LAPPLY = lapply
    }

    # ----------- ANALYSIS ----------
    if(!quiet)
        cat("calculating enrichment\n")

    # calculate Fisher's exact test for provided gene lists
    RESULTS <- .calculateEnrichment(FOREGROUND, BACKGROUND, terms, LAPPLY = LAPPLY,
                                    alternative = alternative,
                                    ...)

    # calculation of length of gene signature
    RESULTS <- RESULTS %>% mutate(signature.length = unlist(lapply(terms, function(x) {
        termTable %>% filter(term == x) %>% nrow
    })))

    # multiple testing correction
    if ( padj.method == "IHW" ) {
        ihwRes <- ihw(p.value ~ signature.length,  data = RESULTS, alpha = padj.cutoff)

        RESULTS <- RESULTS %>% mutate(p.adj = adj_pvalues(ihwRes),
                                      IHW.weights = weights(ihwRes))
    } else {
        RESULTS <- RESULTS %>% mutate(p.adj = p.adjust(p.value, method = padj.method))
    }

    # classification of significance using a cutoff
    RESULTS <- RESULTS %>% mutate(significant = p.adj < padj.cutoff)

    # order for significance and then for odds ratios
    RESULTS <- RESULTS[order(RESULTS %>% .[["significant"]], RESULTS %>% .[["oddsRatio"]], decreasing = TRUE),]

    # ----------- RANDOMIZED BACKGROUND LIST CORRECTION ----------
    if(!quiet)
        cat("calculating randomized background\n")

    if(correction) {
        # create list of randomly selected genes with the length of foreground list
        SAMPLED.LIST <- LAPPLY(1:permutations, function(x) sample_n(BACKGROUND, nrow(FOREGROUND), replace = FALSE))

        # calculate the significance and ranks
        if(!quiet)
            cat("calculate background enrichment\n")
        SAMPLED.RESULTLIST <- LAPPLY(SAMPLED.LIST, function(x) {
            if(!quiet)
                cat(".")
            .calculateEnrichment(x, BACKGROUND, terms, LAPPLY = LAPPLY, alternative = alternative, ...)
        })
        if(!quiet)
            cat("\n")

        # get a table with all the ranks
        SAMPLED.RANKS <- do.call(cbind, lapply(SAMPLED.RESULTLIST, function(x) x %>% select(p.value.rank)))

        # calculate row mean
        RESULTS[, "random.rank.mean"] <- rowMeans(SAMPLED.RANKS)

        # calculate standard deviation
        RESULTS[, "random.rank.sd"] <- apply(SAMPLED.RANKS, 1, sd, na.rm = TRUE)

        # calculate z-score
        if(!quiet)
            cat("calculate z-score")
        RESULTS <- RESULTS %>% mutate(z.score = ( p.value.rank - random.rank.mean ) / random.rank.sd)

        # calculate combined score
        RESULTS <- RESULTS %>% mutate(combined.score = p.value * z.score)
    }

    # ----------- POSTPROCESSING ----------
    if(!quiet)
        cat("adding annotation\n")

    # adding annotation
    if( !is.null(annotation) ) {
        RESULTS <- RESULTS %>% left_join(y = annotation, by = "term")
    }

    return(RESULTS)
}


# provide precalculated table for some annotations
# as argument

.calculateEnrichment <- function(FOREGROUND, UNIVERSE, terms, LAPPLY = lapply, alternative = "greater", ...) {
    # statistical testing
    TESTRESULTS <- LAPPLY(terms, function(x) {
        BACKGROUND <- UNIVERSE %>% filter(!ID %in% FOREGROUND[,"ID"])
        bg.in  <- BACKGROUND %>% filter(term == x) %>% nrow
        bg.out <- (BACKGROUND %>% nrow) - bg.in
        fg.in  <- FOREGROUND %>% filter(term == x) %>% nrow
        fg.out <- (FOREGROUND %>% nrow) - fg.in
        stopifnot(sum(c(bg.in, bg.out, fg.in, fg.out)) != nrow(UNIVERSE))

        ctable            <- matrix(c(fg.in, fg.out, bg.in, bg.out), ncol = 2, nrow = 2)
        significance.test <- fisher.test(ctable, alternative = alternative, ...)
        oddsRatio         <- as.numeric(significance.test$estimate)
        p.value           <- significance.test$p.value

        return(c("oddsRatio"      = oddsRatio,
                 "p.value"        = p.value,
                 "foreground.in"  = fg.in,
                 "foreground.out" = fg.out,
                 "background.in"  = bg.in,
                 "background.out" = bg.out))
    })

    # create a matrix
    TESTRESULTS <- do.call(rbind, TESTRESULTS)

    # add terms and convert to tibble
    TESTRESULTS <- data.frame(term = terms, TESTRESULTS, stringsAsFactors = FALSE) %>% as.tibble

    # calculating ranks based on p-values
    TESTRESULTS <- TESTRESULTS %>% mutate(p.value.rank = rank(p.value))

    TESTRESULTS
}


# calculate the propability that the gene set is pulled out by chance, dependent on the length
# so if a gene set is only 1
# Enrichr implements three approaches to compute enrichment. The first one is a
# standard method implemented within most enrichment analysis tools: the Fisher
# exact test. This is a proportion test that assumes a binomial distribution and
# independence for probability of any gene belonging to any set. The second test
# is a correction to the Fisher exact test that we developed based on intuition.
# We first compute enrichment using the Fisher exact test for many random input
# gene lists in order to compute a mean rank and standard deviation from the
# expected rank for each term in each gene-set library. Then, using a lookup
# table of expected ranks with their variances, we compute a z-score for
# deviation from this expected rank, this can be a new corrected score for
# ranking terms. Alternatively, we combined the p-value computed using the
# Fisher exact test with the z-score of the deviation from the expected rank
# by multiplying these two numbers as follows:
#    c=log(p)·z	(1)
# Where c is the combined score, p is the p-value computed using the Fisher
# exact test, and z is the z-score computed by assessing the deviation from
# the expected rank. Enrichr provides all three options for sorting enriched
# terms. In the results section, we show how we evaluated the quality of each
# of these three enrichment methods by examining how the methods rank terms
# that we know should be highly ranked.


# rank for random selection
# permutation test:
# 1. randomly select X genes from the background and do fisher exact test
# 2. calculate mean rank and standard devidation
# 3. calculate z-score to calculate deviation from expected ranks
# 4. Alternatively, we combined the p-value computed using the
# Fisher exact test with the z-score of the deviation from the expected rank
# by multiplying these two numbers as follows:
#    c=log(p)·z	(1)
#



