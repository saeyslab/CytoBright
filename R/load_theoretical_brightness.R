#' Function to load theoretical brightness
#'
#' @returns A data frame which has a Fluorochrome and a Brightness column
#'
#' @importFrom readxl read_xlsx
#'
#' @export
theoretical_brightness <- function() {
  brightness <- readxl::read_xlsx(
    system.file("extdata/brightness_update210301.xlsx",
      package = "CytoBright"
    ),
    sheet = 1
  )
  brightness <- data.frame(brightness,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  brightness <- brightness[!is.na(brightness$Fluorochrome), ]
  brightness$Fluorochrome <- clean_fluorochromes(brightness$Fluorochrome)
  rownames(brightness) <- brightness$Fluorochrome
  brightness
}
