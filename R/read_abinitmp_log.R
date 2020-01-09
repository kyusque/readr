#' Read common/combined log file into a tibble
#'
#' This is a fairly standard format for log files - it uses both quotes
#' and square brackets for quoting, and there may be literal quotes embedded
#' in a quoted string. The dash, "-", is used for missing values.
#'
#' @inheritParams read_delim
#' @export
#' @examples
#' read_log(readr_example("example.log"))
read_abinitmp_log <- function(file, skip = 0, n_max = Inf, progress = show_progress()) {
  tokenizer <- tokenizer_abinitmp_log()
  read_delimited(file, tokenizer, col_names = c("section", "body"), col_types = cols(section = col_character(), body = col_character()),
    skip = skip, n_max = n_max, progress = progress)
}
