#' Helper function to look up variable information in data dictionary

varinfo <- function(string, dictionary = "dtd", attr = NULL) {
  require(crayon)

  if (!dictionary %in% c("dtd", "drv", "dtd"))
    messaging::emit_error("Not a data dictionary I know of...")

  if(is.null(attr)) {
    dc <- get(dictionary)
    sel <- dc[grep(string, dc$variable, ignore.case = TRUE),]
  }

  return(sel)
}
