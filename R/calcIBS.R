#' Calculate identity by state (IBS) between two or more lines in bins across the genome
#'
#' @param Hapmap - A hapmap file read into R. Needs to have chrom and pos columns named like that, but other names are flexible
#' @param ReferenceLine - The line to which IBS in the others will be compared
#' @param OtherLines - The line or lines that will be used for comparisons. If not specified, all other genotypes will be used
#' @param BinSize - the number of SNPs per bin
#' @param MinBinSize - At the ends of chromosomes or when there are many missing values, bins may be smaller than the number specified by `BinSize`. This number controls how small bins can be before NA is returned. (2 SNPs is probably not enough to give a good estimate of IBS, for example.)
#' @param Chromosomes - A vector of which chromosomes should be used for the calculations. All are used by default
#'
#' @return
#' @export
#'
make.IBS.table <- function(Hapmap, ReferenceLine, OtherLines, BinSize, MinBinSize=10, Chromosomes) {

  if (missing(Chromosomes)) {Chromosomes <- unique(Hapmap$chrom)}
  if (missing(OtherLines)) {OtherLines <- names(Hapmap %>% dplyr::select(-(dplyr::all_of(1:11)), -(dplyr::all_of(ReferenceLine))))}
  if (missing(BinSize)) {BinSize <- (nrow(Hapmap) / length(Chromosomes) / 20) %/% 1}

  AllGenotypes <- append(ReferenceLine, OtherLines)

  MergedGenos <- Hapmap %>%
    dplyr::arrange(chrom, pos) %>%
    dplyr::select(1:10, dplyr::all_of(AllGenotypes))
  LookupDataFrame <- data.frame(Chromosome=double(),
                                BinNumber=double(),
                                BinStartRow=double(),
                                BinEndRow=double(),
                                BinStartPos=double(),
                                BinEndPos=double(),
                                stringsAsFactors=FALSE)

  # I wish this wasn't a for loop, but it runs very fast and a little complicated, so I will leave it

  for (i in Chromosomes) {
    tempdf <- MergedGenos %>% dplyr::filter(chrom == i)
    MarkerNumber <- nrow(tempdf)
    NumberOfBins <- MarkerNumber %/% BinSize
    if ((MarkerNumber %% BinSize) > 0) {
      NumberOfBins <- NumberOfBins + 1
    }
    for (bin in 1:NumberOfBins) {
      StartRow <- (1 + ((bin - 1) * BinSize))
      EndRow <- bin * BinSize
      if (bin == NumberOfBins) {
        EndRow <- nrow(tempdf)
      }
      StartPos <- tempdf[StartRow, "pos"]
      EndPos <- tempdf[EndRow, "pos"]
      newdf <- data.frame(i, bin, StartRow, EndRow, StartPos, EndPos)
      names(newdf) <- c("Chromosome", "BinNumber", "BinStartRow", "BinEndRow", "BinStartPos", "BinEndPos")
      LookupDataFrame <- rbind(LookupDataFrame, newdf)
    }
  }

  GenotypesDataFrame <- data.frame(Genotype=factor(),
                                   Chromosome=double(),
                                   BinNumber=double(),
                                   BinStartPos=double(),
                                   BinEndPos=double(),
                                   ProportionIBS=double(),
                                   stringsAsFactors=FALSE)

  # Heterozygosity Lookup table!

  HetTable <- readr::read_csv(file="A, C, G, T, R, Y, S, W, K, M, N, +, -, 0
                             R, Y, R, Y, A, C, C, A, G, A, N, +, -, 0
                             W, S, S, W, G, T, G, T, T, C, N, +, -, 0
                             M, M, K, K, R, Y, S, W, K, M, N, +, -, 0
                             A, C, G, T, R, Y, S, W, K, M, N, +, -, 0", show_col_types = FALSE)
  HetTable <- as.data.frame(HetTable)

  # This is the slow part

  GenotypesDataFrame <- data.frame(Genotype=factor(),
                                        Chromosome=double(),
                                        BinNumber=double(),
                                        BinStartPos=double(),
                                        BinEndPos=double(),
                                        ProportionIBS=double(),
                                        stringsAsFactors=FALSE)


  for (row in 1:nrow(LookupDataFrame)) { # basically looping over bins
    tempdf <- MergedGenos %>% # Get just the genotypes from this bin
      dplyr::filter(chrom == LookupDataFrame[row, "Chromosome"]) %>%
      dplyr::filter(pos >= LookupDataFrame[row, "BinStartPos"]) %>%
      dplyr::filter(pos <= LookupDataFrame[row, "BinEndPos"])

    # data frame to store the output
    # Why do we reinitialize it for each bin?
    # We (sometimes) make the comparison for multiple taxa on the same bin (depends on the `OtherLines` argument) .
    # But does that require this structure?
    # We have to merge after finishing each individual and each chromosome. Could we use just one object?


    for (col in 12:ncol(tempdf)) { # Column 11 is the reference genotype, so we look over the comparison genotypes
      ProportionIBS <- c(NA)
      NumberOfMarkers <- 0
      for (marker in 1:nrow(tempdf)) {
        if (!is.na(tempdf[marker, col]) & !is.na(tempdf[marker, ReferenceLine]) & tempdf[marker, col] != "N" & tempdf[marker, ReferenceLine] != "N") { # If the SNP call is not missing in either parent
          if (tempdf[marker, ReferenceLine] %in% HetTable[, tempdf[marker, col]]) {
            ProportionIBS <- append(ProportionIBS, 1)
            NumberOfMarkers <- NumberOfMarkers + 1
            next()
          }
          if (tempdf[marker, col] != tempdf[marker, ReferenceLine]) {
            ProportionIBS <- append(ProportionIBS, 0)
            NumberOfMarkers <- NumberOfMarkers + 1
          }
        }
      }
      newdf <- data.frame(names(tempdf)[col], LookupDataFrame[row, "Chromosome"], LookupDataFrame[row, "BinNumber"], LookupDataFrame[row, "BinStartPos"], LookupDataFrame[row, "BinEndPos"], ifelse(NumberOfMarkers > MinBinSize, mean(ProportionIBS, na.rm=TRUE), NA))
      names(newdf) <- c("Genotype", "Chromosome", "BinNumber", "BinStartPos", "BinEndPos", "ProportionIBS")
      GenotypesDataFrame <- rbind(GenotypesDataFrame, newdf)
    }
  }

  return(GenotypesDataFrame)
}
