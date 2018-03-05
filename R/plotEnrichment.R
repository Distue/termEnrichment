#' Displaying enrichment
#'
#' visualization Barplot of log2 odd ratios of significantly enriched RNA types
#'
#' @param X tibble. results table from the termEnrichment function
#' @param nameColumn character. Name of column for the plot to display
#' @param onlySignificant logical. display only signficant results if TRUE. significance is determined by column "significant"
#' @param maxDisplay numeric. maximum number of terms to display
#' @param verbose logical. hides messages if FALSE
#' @param gg logical. use ggplot2 if TRUE. default: TRUE
#' @param ratioType type of ratio to display. "oddsRatio", "ratio" - default "oddsRatio"
#' @param las numeric grapnical parameter las
#' @param oma vector graphical parameter oma
#' @param ... graphical parameters
#' @return ggplot object if gg = TRUE, else plots barchart
#' @author Thomas Schwarzl <schwarzl@embl.de>
#' @import tibble
#' @import dplyr
#' @importFrom stats reorder
#' @import ggplot2
#' @export
#' @usage plotEnrichment(X, nameColumn, onlySignificant,
#'                          maxDisplay, verbose,
#'                          gg, las, oma, ...)
#' @examples
#'
#' data(yeastGO)
#' data(backgroundIDs)
#' data(foregroundIDs)
#' data(yeastGOdesc)
#'
#' X <- termEnrichment(yeastGO, foregroundIDs, backgroundIDs, annotation = yeastGOdesc)
#' plotEnrichment(X)
plotEnrichment <- function(X, nameColumn = ifelse(any(colnames(X) == "name"), "name", "term"), onlySignificant = TRUE,
                           maxDisplay = 20, verbose = TRUE, gg = TRUE,  ratioType = "oddsRatio",
                           las = 2, oma = NULL, ...) {
    if(! (is.character(nameColumn) || is.null(nameColumn)) ) {
        stop("'nameColumn' has to be a character string or NULL if not defined")
    }

    if(is.null(nameColumn)) {
        nameColumn <- "term"
    }

    if(! (is.data.frame(X) || is.tibble(X)) ) {
        stop("'X' has to be a data.frame or a tibble")
    }

    if(! any(colnames(X) == nameColumn))
        stop(paste0("Column names of X do not contain nameColumn (currently: ", nameColumn, ")."))

    if(! ratioType %in% c("oddsRatio", "ratio"))
        stop("ratioType must be either \"oddsRatio\" or \"ratio\"")

    if(! all(c("significant", "p.adj") %in% colnames(X)))
        stop("X must contain columns \"significant\" and \"p.adj\"")

    if(ratioType == "oddsRatio" && ! "oddsRatio" %in% colnames(X))
        stop("X must contain columns \"oddsRatio\"")

    if(ratioType == "ratio" && ! "ratio" %in% colnames(X))
        stop("X must contain columns \"ratio\"")


    # replace all na's or empty spaces with the original term
    call = lazyeval::interp(~ ifelse( is.na(a), term, a), a = as.name(nameColumn))
    X <- X %>% mutate_(.dots = setNames(list(call), nameColumn))

    par.bkup <- par()

    # filters only significant terms
    if( onlySignificant ) {
        X <- X %>% dplyr::filter(significant == TRUE)
    }

    # if there are no terms, display a message if verbose and print an empty plot
    if( nrow(X) <= 0  ) {
        if( verbose ) message("No terms to display.")
        if( gg ) {
            ggplot(data = data.frame(label="No significant results", x=1, y=1), aes(x=x, y=y, label=label)) +
                geom_text() +
                theme_minimal() +
                xlab("") + ylab("") +
                theme(axis.ticks=element_blank(),
                      axis.text=element_blank(),
                      panel.border = element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
        } else {
            plot.new()
            text(x = 0.5, y = 0.5, labels = "No significant results")
        }

    } else {
        # if there are too many terms, just select the maximum to display
        if ( nrow(X) >= maxDisplay ) {
            if( verbose )
                message(paste0(nrow(X), "  to display, displaying the first ",
                               maxDisplay, ". Display more or less using the ",
                               "maxDisplay parameter."))

            X <- X %>% head(n = maxDisplay)
        }

        if( gg ) {
            if(ratioType == "oddsRatio") {
                p <- ggplot(X  %>% arrange(desc(oddsRatio)) %>% mutate(log2oddsRatio = log2(oddsRatio)),
                       aes_string(paste0("reorder(",nameColumn,", -log2oddsRatio)"), y = "log2oddsRatio", fill="p.adj"))
            } else if (ratioType == "ratio") {
                p <- ggplot(X  %>% arrange(desc(ratio)),
                            aes_string(paste0("reorder(",nameColumn,", -ratio)"), y = "ratio", fill="p.adj"))
            }
            p <- p +
                geom_bar(stat = "identity") +
                theme_minimal(base_size = 16) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                ylab(expression("log"[2]*" odds ratio")) +
                xlab("") +
                scale_fill_continuous(guide = guide_legend(title = "p-adjusted"))
            p
        } else {

            if( is.null(oma) ){
                oma = c(10,0,0,0)
            }

            par(lty = 0, oma = oma)

            if(ratioType == "oddsRatio") {
                barplot(log2(X %>% .[["oddsRatio"]]),
                        ylab = expression("log"[2]*" odds ratio"),
                        names = (X %>% .[[nameColumn]]), las = las, ...)
            } else if(ratioType == "ratio") {
                barplot(log2(X %>% .[["ratio"]]),
                        ylab = "ratio",
                        names = (X %>% .[[nameColumn]]), las = las, ...)
            }
            par(lty = par.bkup$lty, oma = par.bkup$oma)
        }
    }
}
