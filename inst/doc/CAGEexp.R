## ----setup, echo = FALSE, results = "hide"------------------------------------
options(signif = 3, digits = 3)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE, fig.height = 5.5,
                      message = FALSE, error = FALSE, warning = TRUE)
set.seed(0xdada)

## ----CAGEprotocol, echo=FALSE, fig.align="center", fig.cap="Overview of CAGE experiment"----
knitr::include_graphics("images/CAGEprotocol.svg")

## ----load_CAGEr---------------------------------------------------------------
library(CAGEr)

## ----specify_input_files------------------------------------------------------
inputFiles = list.files( system.file("extdata", package = "CAGEr")
                       , "ctss$"
                       , full.names = TRUE)

## ----create_CAGEexp-----------------------------------------------------------
ce <- CAGEexp( genomeName     = "BSgenome.Drerio.UCSC.danRer7"
             , inputFiles     = inputFiles
             , inputFilesType = "ctss"
             , sampleLabels   = sub( ".chr17.ctss", "", basename(inputFiles))
)

## ----display_created_object---------------------------------------------------
ce

## ----load_data----------------------------------------------------------------
ce <- getCTSS(ce)
ce

## -----------------------------------------------------------------------------
CTSStagCountSE(ce)
CTSScoordinatesGR(ce)
CTSStagCountDF(ce)
CTSStagCountGR(ce, 1)  # GRanges for one sample with expression count.

## -----------------------------------------------------------------------------
sampleLabels(ce)

## -----------------------------------------------------------------------------
ce <- annotateCTSS(ce, exampleZv9_annot)

## -----------------------------------------------------------------------------
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
CTSScoordinatesGR(ce)

## -----------------------------------------------------------------------------
plotAnnot(ce, "counts")

## ----CorrelationScatterPlots, fig.cap="Correlation of raw CAGE tag counts per TSS"----
corr.m <- plotCorrelation2( ce, samples = "all"
                          , tagCountThreshold = 1, applyThresholdBoth = FALSE
                          , method = "pearson")

## -----------------------------------------------------------------------------
ce <- mergeSamples(ce, mergeIndex = c(3,2,4,4,1), 
                   mergedSampleLabels = c("Zf.unfertilized.egg", "Zf.high", "Zf.30p.dome", "Zf.prim6"))
ce <- annotateCTSS(ce, exampleZv9_annot)

## -----------------------------------------------------------------------------
librarySizes(ce)

## ----ReverseCumulatives, fig.cap="Reverse cumulative distribution of CAGE tags"----
plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)

## -----------------------------------------------------------------------------
ce <- normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)
ce[["tagCountMatrix"]]

## -----------------------------------------------------------------------------
ce <- clusterCTSS( ce
                 , threshold = 1
                 , thresholdIsTpm = TRUE
                 , nrPassThreshold = 1
                 , method = "distclu"
                 , maxDist = 20
                 , removeSingletons = TRUE
                 , keepSingletonsAbove = 5)

## -----------------------------------------------------------------------------
tagClustersGR(ce, sample = "Zf.unfertilized.egg")

## ----CumulativeDistribution, echo=FALSE, fig.align="center", fig.cap="Cumulative distribution of CAGE signal and definition of interquantile width"----
knitr::include_graphics("images/CumulativeDistributionAndQuantiles.svg")

## -----------------------------------------------------------------------------
ce <- cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
ce <- quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

## -----------------------------------------------------------------------------
tagClustersGR( ce, "Zf.unfertilized.egg"
             , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

## -----------------------------------------------------------------------------
plotInterquantileWidth(ce, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

## -----------------------------------------------------------------------------
ce <- aggregateTagClusters(ce, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
ce$outOfClusters / ce$librarySizes *100

## -----------------------------------------------------------------------------
consensusClustersGR(ce)

## -----------------------------------------------------------------------------
ce <- annotateConsensusClusters(ce, exampleZv9_annot)
ce <- cumulativeCTSSdistribution(ce, clusters = "consensusClusters", useMulticore = TRUE)
ce <- quantilePositions(ce, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)

## -----------------------------------------------------------------------------
consensusClustersGR( ce, sample = "Zf.unfertilized.egg"
		               , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

## -----------------------------------------------------------------------------
trk <- exportToTrack(CTSSnormalizedTpmGR(ce, "Zf.30p.dome"))
ce |> CTSSnormalizedTpmGR("all") |> exportToTrack(ce, oneTrack = FALSE)

## -----------------------------------------------------------------------------
split(trk, strand(trk), drop = TRUE)

## ----CTSSbedGraph, echo=FALSE, fig.cap="CAGE data bedGraph track visualized in the UCSC Genome Browser"----
knitr::include_graphics("images/CTSSbedGraph.svg")

## -----------------------------------------------------------------------------
iqtrack <- exportToTrack(ce, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneTrack = FALSE)
iqtrack
#rtracklayer::export.bed(iqtrack, "outputFileName.bed")

## ----tagClustersBed, echo=FALSE, fig.align="center", fig.cap="Tag clusters visualization in the genome browser"----
knitr::include_graphics("images/TagClustersBed.svg")

## -----------------------------------------------------------------------------
ce <- getExpressionProfiles(ce, what = "consensusClusters", tpmThreshold = 10, 
  nrPassThreshold = 1, method = "som", xDim = 4, yDim = 2)

consensusClustersGR(ce)$exprClass |> table(useNA = "always")

## ----ConsensusClustersExpressionProfiles--------------------------------------
plotExpressionProfiles(ce, what = "consensusClusters")

## ----ConsensusClustersExpressionProfiles_svg, echo=FALSE, fig.align="center", fig.cap="Expression clusters"----
knitr::include_graphics("images/ConsensusClustersExpressionProfiles.svg")

## -----------------------------------------------------------------------------
consensusClustersGR(ce) |> subset(consensusClustersGR(ce)$exprClass ==  "0_1")

## -----------------------------------------------------------------------------
# Not supported yet for CAGEexp objects, sorry.
# exportToBed(ce, what = "consensusClusters", 
# 		colorByExpressionProfile = TRUE)

## ----ConsensusClustersBed, echo=FALSE, fig.align="center", fig.cap="Consensus clusters colored by expression profile in the genome browser"----
knitr::include_graphics("images/ConsensusClustersBed.svg")

## -----------------------------------------------------------------------------
ce$group <- factor(c("a", "a", "b", "b"))
dds <- consensusClustersDESeq2(ce, ~group)

## -----------------------------------------------------------------------------
ce <- cumulativeCTSSdistribution(ce, clusters = "consensusClusters")

## -----------------------------------------------------------------------------
# Not supported yet for CAGEexp objects, sorry.
# scoreShift(ce, groupX = "Zf.unfertilized.egg", groupY = "zf_prim6",
# 		testKS = TRUE, useTpmKS = FALSE)

## ----ShiftingScore, echo=FALSE, fig.cap="Calculation of shifting score"-------
knitr::include_graphics("images/ShiftingScore.svg")

## -----------------------------------------------------------------------------
# Not supported yet for CAGEexp objects, sorry.
# shifting.promoters <- getShiftingPromoters(ce, 
# 		tpmThreshold = 5, scoreThreshold = 0.6,
# 		fdrThreshold = 0.01)
# head(shifting.promoters)

## ----ShiftingPromoter, echo=FALSE, fig.cap="Example of shifting promoter"-----
knitr::include_graphics("images/ShiftingPromoter.svg")

## ----create_df----------------------------------------------------------------
TSS.df <- read.table(system.file( "extdata/Zf.unfertilized.egg.chr17.ctss"
                                , package = "CAGEr"))
# make sure the column names are as required
colnames(TSS.df) <- c("chr", "pos", "strand", "Zf.unfertilized.egg")
# make sure the column classes are as required
TSS.df$chr <- as.character(TSS.df$chr)
TSS.df$pos <- as.integer(TSS.df$pos)
TSS.df$strand <- as.character(TSS.df$strand)
TSS.df$Zf.unfertilized.egg <- as.integer(TSS.df$Zf.unfertilized.egg)
head(TSS.df)

## ----coerce_to_df-------------------------------------------------------------
ce.coerced <- as(TSS.df, "CAGEexp")
ce.coerced

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

