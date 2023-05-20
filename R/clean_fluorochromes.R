#' Adapt fluorochrome names to preferred short name
#' from https://github.com/ISAC-DSTF/ProbeTagDictionary/blob/master/ProbeTagDictionary.json
#' @param x Fluorochromes to clean
#' @export
#'
#' @importFrom stringr str_trim
#'
#' @examples
#' clean_fluorochromes(c("Bv-421", "Pacificblue"))
#'
clean_fluorochromes <- function(x) {
  x <- stringr::str_trim(x)

  # Remove prototype suffixes
  x <- gsub("-P2*$", "", x)
  x <- gsub("P2$", "", x)
  # Remove dots, except for Cy5.5 and Cy3.5
  x <- gsub("\\.", "", x)
  x <- gsub("Cy35", "Cy3.5", x)
  x <- gsub("Cy5-5", "Cy5.5", x)
  x <- gsub("Cy55", "Cy5.5", x)
  x <- gsub("Cy5-5", "Cy5.5", x)
  x <- gsub("Cy5_5", "Cy5.5", x)
  x <- gsub("PerCPCy5.5", "PerCP-Cy5.5", x)
  # Remove spaces
  x <- gsub(" ", "", x)
  # Use abbreviations with correct capitalisation
  x <- gsub("Bv", "BV", x)
  x <- gsub("BrilliantViolet", "BV", x)
  x <- gsub("Brilliantviolet", "BV", x)
  x <- gsub("brilliantviolet", "BV", x)

  x <- gsub("yellow", "Yellow", x)
  x <- gsub("RealYellow", "RY", x)

  x <- gsub("blue", "Blue", x)
  x <- gsub("Pacific", "Pac", x)
  x <- gsub("PacB$", "PacBlue", x)


  x <- gsub("Kiravia-Blue", "KiraviaBlue", x)
  x <- gsub("KIRAVIA-Blue", "KiraviaBlue", x)

  x <- gsub("AlexaFluor", "AF", x)
  x <- gsub("Ax", "AF", x)
  x <- gsub("Alexa", "AF", x)

  x <- gsub("Amcyan", "AmCyan", x)
  x <- gsub("AmCyan.*", "AmCyan", x)

  x <- gsub("eFLuor", "eFluor", x)

  x <- gsub("Texas Red", "TxRed", x)
  x <- gsub("TexasRed", "TxRed", x)

  x <- gsub("^Red *([0-9]+)", "R\\1", x)

  x <- gsub("P[Ee] *Cy", "PE-Cy", x)

  x <- gsub("SB", "Super Bright ", x)

  x <- gsub("Qdot *([0-9]+)", "QD\\1", x)

  x <- gsub("^APCF *([0-9]+)", "APC-Fire\\1", x)
  x <- gsub("^APCFire([0-9]+)", "APC-Fire\\1", x)
  x <- gsub("^APCCy([0-9]+)", "APC-Cy\\1", x)

  x <- gsub("Per[Cc]p", "PerCP", x)

  # If letters after number, bring letters to front
  x <- gsub("([A-z]+[0-9]+) *([A-z][A-z]+)", "\\2-\\1", x)

  # Remove trailing spaces
  x <- gsub(" *$", "", x)

  # Remove things between brackets
  x <- gsub("([^ ]) *\\(.*", "\\1", x)

  # Remove space or dash between fluor and number
  x <- gsub("(.*)[ -]([0-9]+)$", "\\1\\2", x)

  # Add dash for some specific dyes
  #x <- gsub("([A-z])-([0-9]*)$", "\\1\\2", x)
  x <- gsub("^(DY|SYTO|POPO|PO-PRO|TOTO|TO-PRO|YOYO|YO-PRO|Hoechst|Fluo|EtHD|Indo|JC)([0-9]*)$",
            "\\1-\\2", x)
  x
}
