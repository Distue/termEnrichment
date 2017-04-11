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
#' @param padj.cutoff cutoff for p.adjusted value
#' @param padj.method p adjustment method for fisher exact test, see ?fisher.test
#' @param alternative alternative hypothesis for fisher exact test, see ?fisher.test
#' @param parallel logical the testing be parallised with BiocParallel
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
#' @importFrom stats reorder
#' @importFrom utils head
#' @return tibble annotated analysis results
#' @export
#' @usage termEnrichment(termTable, foregroundIDs, backgroundIDs, annotation, terms, padj.cutoff, padj.method, alternative, parallel)
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
                           annotation = NULL, terms = NULL,
                           padj.cutoff = 0.05, padj.method = "BH",
                           alternative = "greater", parallel = FALSE) {
    stopifnot(nrow(termTable) > 0 && ncol(termTable) == 2)
    stopifnot(is.null(backgroundIDs) || length(backgroundIDs) > 0)
    stopifnot(length(foregroundIDs) > 0)

    stopifnot(is.data.frame(termTable) || is.tibble(termTable))

    # rename termTable column names
    colnames(termTable) <- c("ID", "term")

    stopifnot(class(termTable) != class(foregroundIDs))
    stopifnot(class(termTable) != class(backgroundIDs))

    # check that no rows are duplicated
    stopifnot(! (termTable %>% duplicated %>% any))

    # if annotation not null it has to be a data frame or a tibble
    if(!is.null(annotation)) {
        stopifnot(is.data.frame(annotation) || is.tibble(annotation))
        stopifnot(any(colnames(annotation) == "term"))
        stopifnot(!(annotation %>% .[["term"]] %>% sort %>% duplicated %>% any))
    }

    # if no background list was specified
    if( is.null(backgroundIDs) ) {
        # set all ids from the mapping table as background
        backgroundIDs = termTable %>% .[["ID"]] %>% sort %>% unique
    } else {
        # assert that all background IDs are found  in the mapping table
        stopifnot(backgroundIDs %in% termTable$ID %>% all)
    }

    # assert that all foreground IDs are found in the mapping table
    stopifnot(foregroundIDs %in% termTable$ID %>% all)

    # assert that all foreground IDs are found in the background list
    stopifnot(foregroundIDs %in% backgroundIDs %>% all)


    # if terms is NULL all terms will be tested
    if(is.null(terms)) {
        # get a list of unique terms
        terms <- termTable %>% .[["term"]] %>% sort %>% unique
        terms <- terms [ terms != "" ]
    } else {
        stopifnot(class(terms) != class(termTable$term))
        stopifnot(terms %in% termTable$term %>% all)
    }

    # generate a subset for background and foreground
    BACKGROUND <- termTable %>% dplyr::filter(ID %in% backgroundIDs)
    FOREGROUND <- termTable %>% dplyr::filter(ID %in% foregroundIDs)

    # parallelize
    if ( parallel ) {
        require( BiocParallel )
        FUN = bplapply
    } else {
        FUN = lapply
    }

    # statistical testing
    RESULTS <- FUN(terms, function(x) {
        bg.in  <- BACKGROUND %>% filter(term == x) %>% nrow
        bg.out <- (backgroundIDs %>% length) - bg.in
        fg.in  <- FOREGROUND %>% filter(term == x) %>% nrow
        fg.out <- (foregroundIDs %>% length) - fg.in

        ctable            <- matrix(c(fg.in, fg.out, bg.in, bg.out), ncol = 2, nrow = 2)
        significance.test <- fisher.test(ctable, alternative = alternative)
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
    RESULTS <- do.call(rbind, RESULTS)

    # add terms and convert to tibble
    RESULTS <- data.frame(term = terms, RESULTS, stringsAsFactors = FALSE) %>% as.tibble

    # multiple testing correction
    RESULTS <- RESULTS %>% mutate(p.adj = p.adjust(p.value, method = padj.method))

    # classification of significance using a cutoff
    RESULTS <- RESULTS %>% mutate(significant = p.adj < padj.cutoff)

    # order for significance and then for odds ratios
    RESULTS <- RESULTS[order(RESULTS %>% .[["significant"]], RESULTS %>% .[["oddsRatio"]], decreasing = TRUE),]

    if( !is.null(annotation) ) {
        RESULTS <- RESULTS %>% left_join(y = annotation, by = "term")
    }

    return(RESULTS)
}


