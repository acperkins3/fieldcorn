#' Make IBS Table
#'
#' @param Hapmap - A hapmap file read into R. Needs to have chrom and pos columns named like that, but other names are flexible
#' @param ReferenceLine - The line to which IBS in the others will be compared
#' @param OtherGenotypes - The line or lines that will be used for comparisons
#' @param BinSize
#'
#' @return
#' @export
#'
make.IBS.table <- function(Hapmap, ReferenceLine, OtherGenotypes, BinSize) {

  AllGenotypes <- append(ReferenceLine, OtherGenotypes)

  Hapmap %>%
    dplyr::arrange(chrom, pos) %>%
    dplyr::select(1:10, dplyr::all_of(AllGenotypes))
  LookupDataFrame <- data.frame(Chromosome=double(),
                                BinNumber=double(),
                                BinStartRow=double(),
                                BinEndRow=double(),
                                BinStartPos=double(),
                                BinEndPos=double(),
                                stringsAsFactors=FALSE)

  for (i in 1:10) {
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
                                   PercentIBS=double(),
                                   stringsAsFactors=FALSE)

  # Heterozygosity Lookup table!

  HetTable <- dplyr::read_csv(file="A, C, G, T, R, Y, S, W, K, M, N, +, -, 0
                             R, Y, R, Y, A, C, C, A, G, A, N, +, -, 0
                             W, S, S, W, G, T, G, T, T, C, N, +, -, 0
                             M, M, K, K, R, Y, S, W, K, M, N, +, -, 0
                             A, C, G, T, R, Y, S, W, K, M, N, +, -, 0")
  HetTable <- as.data.frame(HetTable)

  # This is the slow one
  #I'm gonna make it NA if there are fewer than 100 markers
  for (row in 1:nrow(LookupDataFrame)) {
    tempdf <- MergedGenos %>%
      dplyr::filter(chrom == LookupDataFrame[row, "Chromosome"]) %>%
      dplyr::filter(pos >= LookupDataFrame[row, "BinStartPos"]) %>%
      dplyr::filter(pos <= LookupDataFrame[row, "BinEndPos"])
    OtherGenotypesDataFrame <- data.frame(Genotype=factor(),
                                          Chromosome=double(),
                                          BinNumber=double(),
                                          BinStartPos=double(),
                                          BinEndPos=double(),
                                          PercentIBS=double(),
                                          stringsAsFactors=FALSE)
    for (col in 12:ncol(tempdf)) {
      PercentIBS <- c(NA)
      NumberOfMarkers <- 0
      for (marker in 1:nrow(tempdf)) {
        if (is.na(tempdf[marker, col])) {
          next()
        }
        if (is.na(tempdf[marker, ReferenceLine])) {
          next()
        }
        if (tempdf[marker, col] == "N") {
          next()
        }
        if (tempdf[marker, ReferenceLine] == "N") {
          next()
        }
        if (tempdf[marker, ReferenceLine] %in% HetTable[, tempdf[marker, col]]) {
          PercentIBS <- append(PercentIBS, 1)
          NumberOfMarkers <- NumberOfMarkers + 1
          next()
        }
        if (tempdf[marker, col] != tempdf[marker, ReferenceLine]) {
          PercentIBS <- append(PercentIBS, 0)
          NumberOfMarkers <- NumberOfMarkers + 1
        }
      }
      newdf <- data.frame(names(tempdf)[col], LookupDataFrame[row, "Chromosome"], LookupDataFrame[row, "BinNumber"], LookupDataFrame[row, "BinStartPos"], LookupDataFrame[row, "BinEndPos"], ifelse(NumberOfMarkers > 100, mean(PercentIBS, na.rm=TRUE), NA))
      names(newdf) <- c("Genotype", "Chromosome", "BinNumber", "BinStartPos", "BinEndPos", "PercentIBS")
      OtherGenotypesDataFrame <- rbind(OtherGenotypesDataFrame, newdf)
    }
    GenotypesDataFrame <- rbind(GenotypesDataFrame, OtherGenotypesDataFrame)
  }

  return(GenotypesDataFrame)
}
