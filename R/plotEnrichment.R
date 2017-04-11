#' Displaying enrichment
#'
#' visualization Barplot of log2 odd ratios of significantly enriched RNA types
#'
#' @param X tibble. results table from the termEnrichment function
#' @param nameColumn character. Name of column for the plot to display
#' @param onlySignificant logical. display only signficant results if TRUE. significance is determined by column "significant"
#' @param maxDisplay numeric. maximum number of terms to display
#' @param verbose logical. hides messages if FALSE
#' @param gg logical. use ggplot2 if TRUE. dafault: TRUE
#' @param las numeric grapnical parameter las
#' @param oma vector graphical parameter oma
#' @param ... graphical parameters
#' @return ggplot object if gg = TRUE, else plots barchart
#' @author Thomas Schwarzl <schwarzl@embl.de>
#' @import tibble
#' @import dplyr
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
plotEnrichment <- function(X, nameColumn = "name", onlySignificant = TRUE,
                           maxDisplay = 20, verbose = TRUE, gg = TRUE,
                           las = 2, oma = NULL, ...) {
    stopifnot(is.character(nameColumn))
    stopifnot(is.data.frame(X) || is.tibble(X))
    stopifnot(any(colnames(X) == nameColumn))
    stopifnot(all(c("oddsRatio", "significant", "p.adj") %in% colnames(X)))

    par.bkup <- par()

    # filters only significant terms
    if( onlySignificant ) {
        X <- X %>% dplyr::filter(significant == TRUE)
    }

    # if there are no terms
    if( nrow(X) <= 0  ) {
        if( verbose ) message("No terms to display.")
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

            ggplot(X, aes(reorder(name, -oddsRatio), log2(oddsRatio), fill=p.adj)) +
                geom_bar(stat = "identity") +
                theme_minimal(base_size = 16) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                ylab(expression("log"[2]*" odds ratio")) +
                xlab("") +
                scale_fill_continuous(guide = guide_legend(title = "p-adjusted"))

        } else {

            if( is.null(oma) ){
                oma = c(10,0,0,0)
            }

            par(lty = 0, oma = oma)
            barplot(log2(X %>% .[["oddsRatio"]]),
                    ylab = expression("log"[2]*" odds ratio"),
                    names = (X %>% .[[nameColumn]]), las = las, ...)

            par(lty = par.bkup$lty, oma = par.bkup$oma)
        }
    }
}
