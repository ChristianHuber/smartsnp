# smartsnp v.1
# Coding start date = 02/07/2020
# smartsnp::smart_permanova by Salvador Herrando-Perez (salherra@gmail.com) and Ray Tobler (tingalingx@gmail.com)
#
# Permutational Multivariate Analysis of Variance (PERMANOVA) of genotype data
# Computes global and pairwise MANOVA tests for differences in GROUP LOCATION for multiple types of proximity measures between samples
# Alpha values estimated by permutations
# Implements scaling by centering, z-scores and SMARTPCA controlling for genetic drift
# Original description of PERMANOVA at https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x / Anderson (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26, 32-46

#' @name smart_permanova
#' @title Smart Permutational Multivariate Analysis of Variance
#'
#' @description Computes Permutational Multivariate Analysis of Variance (PERMANOVA) for testing differences in group location using multivariate data. Variance partitioning computed on a sample-by-sample triangular matrix obtained from variable-by-sample data following Anderson (2001).
#' Calculates a range of inter-sample distances, similarities and dissimilarities.
#' Includes control for genetic drift for bi-allelic genetic markers such as single nucleotide polymorphisms (SNP) following Patterson, Price and Reich (2006) that can be combined with SMART Principal Component Analysis (PCA). Optimized to run fast matrix building and permutations for big datasets in ecological, evolutionary and genomic research.
#'
#' @details PERMANOVA is a form of linear modelling that partitions variation in a triangular matrix of inter-sample proximities obtained from variable-by-sample data.
#' Uses permutations to estimate the probability of observed group differences in SNP composition given a null hypothesis of no differences between groups (Anderson 2001).
#' Proximity between samples can be any type of distance, similarity or dissimilarity.
#' Original acronym \code{NPMANOVA} (Non-Parametric MANOVA) replaced with PERMANOVA (Anderson 2004, 2017).\cr
#'
#' Univariate ANOVA captures differences in mean and variance referred to as location and dispersion in PERMANOVA's multivariate context (Anderson & Walsh 2013, Warton, Wright and Wang 2012).
#' To attribute group differences to location (position of sample groups) and/or dispersion (spread of sample groups), PERMANOVA must be combined with PERMDISP as implemented through \code{smart_permdisp}.\cr
#'
#' Function \code{smart_permanova} uses \code{\link[vegan]{adonis}} to fit formula \code{snp_eucli ~ sample_group}, where \code{snp_eucli} is the sample-by-sample triangular matrix in Principal Coordinate Analysis (Gower 1966) space.
#' Current version restricted to one-way designs (one categorical predictor) though PERMANOVA can handle >1 crossed and/or nested factors (Anderson 2001) and continuous predictors (McArdle & Anderson 2001).
#' If >2 sample groups tested, \code{pairwise = TRUE} allows pairwise testing and correction for multiple testing by \code{holm (Holm)} [default], \code{hochberg (Hochberg)}, \code{hommel (Hommel)}, \code{bonferroni (Bonferroni)}, \code{BY (Benjamini-Yekuieli)}, \code{BH (Benjamini-Hochberg)} or \code{fdr (False Discovery Rate)}.\cr
#'
#' For big data, \code{\link[Rfast]{Dist}} builds sample-by-sample triangular matrix much faster than \code{\link[vegan]{vegdist}}.
#' \code{\link[Rfast]{Dist}} computes proximities \code{euclidean}, \code{manhattan}, \code{canberra1}, \code{canberra2}, \code{minimum}, \code{maximum}, \code{minkowski}, \code{bhattacharyya}, \code{hellinger}, \code{kullback_leibler} and \code{jensen_shannon}. \code{\link[vegan]{vegdist}} computes \code{manhattan}, \code{euclidean}, \code{canberra}, \code{clark}, \code{bray}, \code{kulczynski}, \code{jaccard}, \code{gower}, \code{altGower}, \code{morisita}, \code{horn}, \code{mountford}, \code{raup}, \code{binomial}, \code{chao}, \code{cao} and \code{mahalanobis}.
#' Euclidean distance required for SMARTPCA scaling.\cr
#'
#' \code{sample_remove} should include both samples removed from PCA and ancient samples projected onto PCA space (if any).\cr
#'
#' Data read from working directory with SNPs as rows and samples as columns.
#' Two alternative formats: (1) text file of SNPs by samples (file extension and column separators recognized automatically) read using \code{\link[data.table]{fread}}; or (2) duet of \code{EIGENSTRAT} files (see \url{https://reich.hms.harvard.edu/software}) using \code{\link[vroom]{vroom_fwf}}, including a genotype file of SNPs by samples (\code{*.geno}), and a sample file (\code{*.ind}) containing three vectors assigning individual samples to unique user-predefined groups (populations), sexes (or other user-defined descriptor) and alphanumeric identifiers.
#' For \code{EIGENSTRAT}, vector \code{sample_group} assigns samples to groups retrievable from column 3 of file \code{*.ind}.
#' SNPs with zero variance removed prior to SVD to optimize computation time and avoid undefined values if \code{scaling = "sd"} or \code{"drift"}.\cr
#'
#' Users can select subsets of samples or SNPs by introducing a vector including column numbers for samples (\code{sample_remove}) and/or row numbers for SNPs (\code{snp_remove}) to be removed from computations.
#' Function stops if the final number of SNPs is 1 or 2.
#' \code{EIGENSOFT} was conceived for the analysis of human genes and its SMARTPCA suite so accepts 22 (autosomal) chromosomes by default.
#' If >22 chromosomes are provided and the internal parameter \code{numchrom} is not set to the target number chromosomes of interest, SMARTPCA automatically subsets chromosomes 1 to 22.
#' In contrast, \code{smart_permanova} accepts any number of autosomes with or without the sex chromosomes from an \code{EIGENSTRAT} file.\cr
#'
#' @param snp_data File name read from working directory.
#' SNP = rows, samples = columns without row names or column headings.
#' SNP values must be count data (no decimals allowed).
#' File extension detected automatically whether text or \code{EIGENSTRAT}.
#' See details.
#' @param packed_data Logical value for \code{EIGENSTRAT}, irrelevant for text data.
#' Default \code{packed_data = FALSE} assumes uncompressed \code{EIGENSTRAT}.
#' \code{packed_data = TRUE} for compressed or binary \code{EIGENSTRAT} (\code{PACKENDANCESTRYMAP}).
#' @param sample_group Character or numeric vector assigning samples to groups. Coerced to factor.
#' @param sample_remove Logical \code{FALSE} or numeric vector indicating column numbers (samples) to be removed from computations.
#' Default \code{sample_remove =  FALSE} keeps all samples.
#' @param snp_remove Logical \code{FALSE} or numeric vector indicating row numbers (SNPs) to be removed from computations.
#' Default \code{snp_remove =  FALSE} keeps all SNPs. See details.
#' @param missing_value Number \code{9} or string \code{NA} indicating missing value.
#' Default \code{missing_value = 9} as in \code{EIGENSTRAT}.
#' If no missing values present, no effect on computation.
#' @param missing_impute String handling missing values.
#' Default \code{missing_impute = "remove"} removes SNPs with at least one missing value.
#' \code{missing_impute = "mean"} replaces missing values of each SNP by mean of non-missing values across samples.
#' If no missing values present, no effect on computation.
#' @param scaling String. Default \code{scaling = "drift"} scales SNPs to control for expected allele frequency dispersion caused by genetic drift (SMARTPCA).
#' \code{scaling = "center"} for \code{centering} (covariance-based PCA).
#' \code{scaling = "sd"} for \code{centered} SNPs divided by standard deviation (correlation-based PCA).
#' \code{scaling = "none"} for no scaling.
#' See details.
#' @param sample_distance Type of inter-sample proximity computed (distance, similarity, dissimilarity).
#' Default is \code{Euclidean distance}. See details.
#' @param program_distance A string value indicating R package to estimate proximities between pairs of samples.
#' Default \code{program_distance = "Rfast"} uses function \code{\link[Rfast]{Dist}}; \code{program_distance = "vegan"} uses \code{\link[vegan]{vegdist}}.
#' See details.
#' @param target_space String.
#' Default \code{target_space = "multidimensional"} applies PERMANOVA to sample-by-sample triangular matrix computed from variable-by-sample data, \code{pc_axes} has no effect on computation.
#' \code{target_space = "pca"} applies PERMANOVA to sample-by-sample data in PCA space, \code{pc_axes} determines number of PCA axes for testing.
#' @param pc_axes Number of PCA axes computed always starting with PCA axis 1. Default \code{pc_axes = 2} computes PCA axes 1 and 2 if \code{target_space = "pca"}.
#' No effect on computation if \code{target_space = "multidimensional"}.
#' @param pairwise Logical.
#' Default \code{pairwise = FALSE} computes global test.
#' \code{pairwise = TRUE} computes global and pairwise tests.
#' @param pairwise_method String specifying type of correction for multiple testing.
#' Default \code{"holm"}.
#' See details.
#' @param permutation_n Number of permutations resulting in PERMANOVA test \emph{p value}.
#' Default \code{9999}.
#' @param permutation_seed Number fixing random generator of permutations.
#' Default \code{1}.
#'
#' @return Returns a list containing the following elements:
#' \itemize{
#'   \item{permanova.samples}{Dataframe showing sample summary.
#'   Column \emph{Group} assigns samples to tested groups.
#'   Column \emph{Class} specifies if samples were used in, or removed from, testing.}
#'   \item{permanova.global_test}{List showing table with degrees of freedom, sum of squares, mean sum of squares, \emph{F} statistic, variance explained (\emph{R2}) and \emph{p} value.}
#'   \item{permanova.pairwise_test}{List showing table \emph{F} statistic, variance explained (\emph{R2}), \emph{p} value and corrected \emph{p} value per pair of groups.
#'   Obtained only if \code{pairwise = TRUE}.}
#'   \item{permanova.pairwise_correction}{String indicating type of correction for multiple testing.}
#'   \item{permanova.permutation_number}{Number of permutations applied to obtain the distribution of \emph{p value}.}
#'   \item{permanova.permutation_seed}{Number fixing random generator of permutations for reproducibility of results.}
#' }
#'
#' @examples
#' # Path to example genotype matrix "dataSNP"
#' pathToGenoFile = system.file("extdata", "dataSNP", package = "smartsnp")
#'
#' # Assign 50 samples to each of two groups
#' my_groups <- as.factor(c(rep("A", 50), rep("B", 50)))
#'
#' # Run PERMANOVA
#' permanovaR <- smart_permanova(snp_data = pathToGenoFile, sample_group = my_groups)
#'
#' # Extract summary table assigning samples to groups
#' permanovaR$permanova.samples
#'
#' # Extract PERMANOVA table
#' permanovaR$permanova.global_test
#'
#' # Plot means of squares per group
#' #run pca with truncated SVD (PCA 1 x PCA 2)
#' pcaR1 <- smart_pca(snp_data = pathToGenoFile, sample_group = my_groups)
#' #compute Euclidean inter-sample distances in PCA space (triangular matrix)
#' snp_eucli <- vegan::vegdist(pcaR1$pca.sample_coordinates[,c("PC1","PC2")], method = "euclidean")
#' #run PERMANOVA
#' permanova <- vegan::adonis(formula = snp_eucli ~ my_groups, permutations = 9999)
#' #extract meanSqs (groups versus residuals)
#' meanSqs <- as.matrix(t(permanova$aov.tab$MeanSqs[1:2]))
#' colnames(meanSqs) <- c("Groups", "Residuals")
#' #two horizontal plots
#' oldpar <- par(mfrow = c(2,1), oma = c(0,5,0.1,0.1), lwd = 2)
#' barplot(meanSqs, horiz = TRUE, main = "PERMANOVA mean of squares",
#'   cex.names = 2, cex.main = 2, col = c("grey40"))
#' #run ANOSIM
#' anosimD <- vegan::anosim(snp_eucli, my_groups, permutations = 999)
#' #remove outputs for clean plotting
#' #anosimD[2] <- ""; anosimD[5] <- ""
#' par(mar = c(5, 0.1, 3.5, 0.1))
#' plot(anosimD, xlab = "", ylab = "distance/similarity ranks",
#'   main = "Inter-sample proximity ranks", cex.main =2, cex.axis = 2,
#'   col = c("cyan", "red", "blue"))
#' par(oldpar)
#'
#' @references
#' Anderson, M. J. (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26, 32-46.\cr
#' Anderson, M. J. (2004). PERMANOVA_2factor: a FORTRAN computer program for permutational multivariate analysis of variance (for any two-factor ANOVA design) using permutation tests  (Department of Statistics, University of Auckland, New Zealand).\cr
#' Anderson, M. J. & D. C. I. Walsh (2013) PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological Monographs, 83, 557-574.\cr
#' Gower, J. C. (1966) Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika, 53, 325-338.\cr
#' McArdle, B. H. & M. J. Anderson (2001) Fitting multivariate models to community data: a comment on distance-based redundancy analysis. Ecology, 82, 290-297.\cr
#' Patterson, N., A. L. Price and D. Reich (2006) Population structure and eigenanalysis. PLoS Genetics, 2, e190.\cr
#' Warton, D. I., S. T. Wright and Y. Wang (2012) Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution, 3, 89-101.
#'
#' @seealso \code{\link[vegan]{adonis}} (package \bold{vegan}),
#' \code{\link[Rfast]{Dist}} (package \bold{Rfast}),
#' \code{\link[data.table]{fread}} (package \bold{data.table}),
#' \code{\link[vegan]{vegdist}} (package \bold{vegan}),
#' \code{\link[vroom]{vroom_fwf}} (package \bold{vroom})
#'
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom foreach %do%
#' @export
utils::globalVariables(c("i"))
smart_permanova <- function(snp_data, packed_data = FALSE,
                            sample_group, sample_remove = FALSE, snp_remove = FALSE,
                            missing_value = 9, missing_impute = "remove",
                            scaling = "drift",
                            sample_distance = "euclidean", program_distance = "Rfast",
                            target_space = "multidimensional", pc_axes = 2,
                            pairwise = FALSE, pairwise_method = "holm",
                            permutation_n = 9999, permutation_seed = 1) {

  ##--------------------------------------------------------------##
  # 1. Timing analytical steps
  ##--------------------------------------------------------------##

  # Overall starting time
  startT <- Sys.time()

  # cats cumulative run time at specific points
  get.time <- function(t) { # time conversion
    tt <- as.numeric(difftime(Sys.time(), t, units = "secs"))
    H <- tt %/% 3600
    rr <- tt %% 3600
    if(rr > 0) {
      M <- rr %/% 60
      rr2 <- rr %% 60
      if(rr2 > 0) {
        S <- round(rr2)
      }else{
        S <- 0
      }
    } else {
      M <- 0
      S <- 0
    }
    return(paste0(H, "h ", M, "m ", S, "s"))
  }

  ##--------------------------------------------------------------##
  # 2. Load data and filter samples and SNPs
  ##--------------------------------------------------------------##

  message("Checking argument options selected...")
  # Check options: convert to logical if needed
  if (!is.numeric(sample_remove)) sample_remove <- as.logical(sample_remove)
  # Print error and abort analysis
  if (missing_impute != "remove" & missing_impute != "mean") { # if users misspel missing_impute
    stop("Check spelling: missing_impute must be 'remove' or 'mean'")
  }
  if (missing_value != 9 & !is.na(missing_value)) { # if users assign missing values not equal to 9 or NA
    stop("Missing values must be coded as 9 or NA")
  }
  if (scaling != "none" & scaling != "center" & scaling != "sd" & scaling != "drift") { # # if users misspel scaling
    stop("Check spelling: scaling must be 'none', 'center', 'sd' or 'drift'")
  }
  if (program_distance != "Rfast" & program_distance != "vegan") { # if users misspel program_distance
    stop("Check spelling: program_distance must be 'Rfast' or 'vegan'")
  }
  if (target_space != "multidimensional" & target_space != "pca") { # if users misspel target_space
    stop("Check spelling: target_space must be 'multidimensional' or 'pca'")
  }
  message("Argument options are correct...")

  message("Loading data...")

  # Label samples
  samp_dat <- data.table::data.table(sample_group) # tabulate group labels
  data.table::setnames(samp_dat, "Group") # name column vector of group labels
  samp_dat[, I:= 1:.N] # cbind column vectors with samples listed from 1 to length(sample_group)
  sampN.full <- nrow(samp_dat) # total number of samples in dataset

  # Define samples used for PERMANOVA
  sample_PCA <- 1:sampN.full # assume all samples for PERMANOVA as default
  if (!isFALSE(sample_remove)) {
    sample_PCA <- setdiff(sample_PCA, sample_remove) # remove excluded samples
  }
  # Recalculate number of samples for PERMANOVA
  sampN <- length(sample_PCA)

  # Avoid error if users unfamiliar with data formats set packed_data = FALSE but their data is binary (packedancestrymap format) following https://stackoverflow.com/questions/16350164/native-method-in-r-to-test-if-file-is-ascii
  is.binary <- function(filepath, max = 1000) {
    f = file(filepath, "rb", raw = TRUE)
    b = readBin(f, "int", max, size = 1, signed = FALSE)
    close(f)
    return(max(b) > 128)
  }
  if (is.binary(snp_data) == TRUE) {
    packed_data = TRUE
    message("Data is binary (packedancestrymap)...")
  } else {
    packed_data = FALSE
  }

  # Read in data (columns = SNPs, rows = samples)
  if(length(grep("\\.geno$", snp_data)) == 1) { # eigenstrat format
    if(isFALSE(packed_data)) { # decompressed format
      snp_dat <- vroom::vroom_fwf(file = snp_data, # fast-loading data as tibble
                                  col_positions = vroom::fwf_widths(rep(1, sampN.full), col_names = NULL),
                                  col_types = paste(rep("i", sampN.full), collapse=""))
      snpN.full <- nrow(snp_dat) # number of SNP
      message(paste("Imported", snpN.full, "SNP by", sampN.full, "sample eigenstrat genotype matrix (decompressed or unpacked format)"))
      message(paste0("Time elapsed: ", get.time(startT)))
      message("Filtering data...")
      #convert to matrix
      snp_dat1 <- do.call(cbind, snp_dat[, sample_PCA])
    } else { # compressed binary format
      ss <- strsplit(snp_data, "\\.")[[1]]
      snp_data <- paste(ss[-length(ss)], collapse=".") # remove suffix
      snp_dat <- read_packedancestrymap(pref = snp_data)$geno
      attr(snp_dat, "dimnames") <- NULL # remove row and column name attributes
      snpN.full <- nrow(snp_dat) # number of SNP
      message(paste("Imported", snpN.full, "SNP by", sampN.full, "sample packed eigenstrat genotype matrix (compressed or packed format)"))
      message(paste0("Time elapsed: ", get.time(startT)))
      message("Filtering data...")
      snp_dat1 <- snp_dat[, sample_PCA]
      missing_value <- NA # reset missing value
    }
  } else { # generic input type (columns = samples, rows = SNPs)
    con <- file(snp_data,"r"); first_line <- readLines(con,n=1); close(con) # Read first line
    plink_traw_format_flag <- FALSE
    if (substr(first_line, 1, 8) == "CHR\tSNP\t") plink_traw_format_flag <- TRUE # Check for PLINK "traw" header line
    snp_dat <- data.table::fread(file = snp_data, header = plink_traw_format_flag)
    if (plink_traw_format_flag) snp_dat[, c("CHR","SNP", "(C)M", "POS", "COUNTED", "ALT"):=NULL] # If PLINK "traw" format, then remove non-genotype columns
    snpN.full <- nrow(snp_dat) # number of SNP
    message(paste("Imported", snpN.full, "SNP by", sampN.full, "sample genotype matrix"))
    message(paste0("Time elapsed: ", get.time(startT)))
    message("Filtering data...")
    snp_dat1 <- do.call(cbind, snp_dat[, sample_PCA, with = FALSE])
  }

  # Print error if users enter number of sample labels (sample_group) larger or smaller than number of samples in dataset (snp_dat)
  if (length(sample_group) != ncol(snp_dat)) {
    stop("length(sample_group) should be equal to number of samples in dataset: computation aborted")
  }

  # Remove original dataset from memory
  rm(snp_dat); gc()

  # Filter SNPS
  if (!isFALSE(snp_remove)) {
    snp.keep <- setdiff(1:snpN.full, snp_remove)
    snpN.full <- length(snp.keep) # update SNP count
    snp_dat1 <- snp_dat1[snp.keep, ] # subset SNPs by row number across modern samples
  }

  if (snpN.full < 3) {
    stop("Less than 3 SNPs remaining: computation aborted")
  }

  message(paste(snpN.full, "SNPs included in PERMANOVA computation"))
  message(paste(length(snp_remove), "SNPs omitted from PCA computation"))

  # Print number of samples used in PERMANOVA
  message(paste(length(sample_PCA), "samples included in PERMANOVA computation"))
  if (!isFALSE(sample_remove)) {
    message(paste(length(sample_remove), "samples ommitted from PERMANOVA computation"))
  }
  message("Completed data filtering")
  message(paste0("Time elapsed: ", get.time(startT)))

  ##--------------------------------------------------------------##
  # 3. Remove invariant SNPs
  ##--------------------------------------------------------------##

  message("Scanning for invariant SNPs...")

  # Identify invariant SNPs allowing for missing data
  if (is.na(missing_value)) { # missing values are NAs
    genoMean <- rowMeans(snp_dat1, na.rm = TRUE) # compute SNP means
    genoVar <- Rfast::rowVars(snp_dat1, na.rm = TRUE) # compute SNP variances
  } else { # missing values are numeric
    genoTab <- Rfast::rowTabulate(snp_dat1) # count frequency of values per SNP
    genoTab <- cbind(sampN - rowSums(genoTab), genoTab) # include 0's
    GTR <- setdiff(1:ncol(genoTab) - 1, missing_value)
    sumN <- rowSums(genoTab[, GTR+1])
    sumX <- foreach::foreach(i = GTR, .combine = "+") %do% {
      i * genoTab[, i+1]
    }
    sumX2 <- foreach::foreach(i = GTR, .combine = "+") %do% {
      i^2 * genoTab[, i+1]
    }
    genoMean <- sumX / sumN
    genoVar <- ((sumX2 / sumN) - genoMean ^ 2) * (sumN / (sumN - 1))
    rm(sumN, sumX, sumX2)
  }
  keepSNPs <-  which(genoVar != 0) # row number for SNPs with > 0 variance #SHP

  # Remove invariant SNPs
  if (length(keepSNPs > 0)) {
    snp_dat1 <- snp_dat1[keepSNPs, ] # remove invariant SNPs from modern samples
    genoMean <- genoMean[keepSNPs] # compute SNP means
    genoVar <- genoVar[keepSNPs] # compute SNP variances
    # if (!isFALSE(sample_project)) {
    #   snp_dat2 <- snp_dat2[keepSNPs, ] # remove invariant SNPs from ancient samples
    # }
  }
  rm(keepSNPs) # remove vector with indices for variant SNPs
  snpN <- nrow(snp_dat1) # number of SNPs with > 0 variance

  if (snpN == snpN.full) {
    message("Scan complete: no invariant SNPs found")
  } else {
    message(paste("Scan complete: removed", snpN.full - snpN, "invariant SNPs"))
  }
  message(paste0("Time elapsed: ", get.time(startT)))

  ##--------------------------------------------------------------##
  # 4. Deal with missing values (SNP removal or imputation)
  ##--------------------------------------------------------------##

  message("Checking for missing values...")

  # Index cells with missing value (counting sequence = top to bottom then left to right sequence)
  if (is.na(missing_value)) { # missing_value = NA
    missI <- which(is.na(snp_dat1))
  } else { # missing_value = numeric
    missI <- which(snp_dat1 == missing_value)
  }

  # Impute SNP with missing values with mean genotype
  if (length(missI) > 0) {
    # Determine SNPs with missing values
    Z <- missI %% snpN # label SNPs (rows) with missing values with 0
    Z[Z == 0] <- snpN
    snpMissI <- sort(unique(Z)) # row index for SNPs with missing values
    snpMissN <- length(snpMissI) # total SNPs with missing values
    message(paste(snpMissN, "SNPs contain missing values"))
    # Mean imputation
    if (missing_impute == "mean") {
      message("Imputing mean for SNPs with missing values...")
      snp_dat1[missI] <- genoMean[Z]
      message(paste("Imputation completed:", length(missI), "missing values imputed"))
    }
    # Removal of SNPs with missing data
    if (missing_impute == "remove") {
      message("Removing SNPs with missing data...")
      snp_dat1 <- snp_dat1[-snpMissI, ] #SNP
      genoMean <- genoMean[-snpMissI] # update genotype means
      genoVar <- genoVar[-snpMissI] # update genotype variance
      snpN <- nrow(snp_dat1)
      message(paste("Removal completed:", snpMissN, "SNPs removed"))
      message(paste(snpN, "SNPs remaining"))
    }
    rm(missI, Z)
  } else {
    message("Scan completed: no missing values found")
    rm(missI) # remove vector with cell number when missing value present
  }
  message(paste0("Time elapsed: ", get.time(startT)))

  ##--------------------------------------------------------------##
  # 5. Scale values by SNP
  ##--------------------------------------------------------------##

  if(scaling != "none") { # only run if standardization is requested
    message("Scaling values by SNP...")

    if (scaling == "drift") { # standardization controlling for genetic drift following Formula (3) by Patterson et al. 2006 (doi: 10.1371/journal.pgen.0020190)
      message("Centering and scaling by drift dispersion...")
      snp_drift <- sqrt((genoMean / 2) * (1 - genoMean/2))
      snp_dat1 <- (snp_dat1 - genoMean) / snp_drift
    }

    if (scaling == "sd") { # z-score standardization whereby all SNPs have mean = 0 and variance = 1
      message("Centering and standardizing (z-scores)...")
      if (missing_impute == "mean") {
        snp_sd <- Rfast::rowVars(snp_dat1, std = TRUE) # update variances to account for imputation
      } else { # calculate sd on precomputed variance after removing invariant snps
        snp_sd <- sqrt(genoVar)
      }
      snp_dat1 <- (snp_dat1 - genoMean) / snp_sd
    }

    if (scaling == "center") { # no standardization: this option runs a conventional PCA on raw (centered) data
      message("Centering...")
      snp_dat1 <- snp_dat1 - genoMean
    }
    message("Completed scaling")
    message(paste0("Time elapsed: ", get.time(startT)))
  } else {
    message("Note: SNP-based scaling not used")
  }

  ##--------------------------------------------------------------##
  # 6. Construct triangular inter-sample distance matrix
  ##--------------------------------------------------------------##

  # Compute sample by sample matrix of proximities using Rfast::Dist (fast) or vegan::vgdist (slow)
  if (target_space == "multidimensional") { # proximity matrix computed from raw data (SNP x sample)
    message("Construct triangular matrix of sample by sample proximities in multidimensional space...")
    if (program_distance == "Rfast") { # Rfast
      snp_eucli <- Rfast::Dist(t(snp_dat1), method = sample_distance) # default: method = "euclidean" (= Euclidean distance)
    }
    if (program_distance == "vegan") { # vegan
      snp_eucli <- vegan::vegdist(t(snp_dat1), method = sample_distance) # default: method = "euclidean" (= Euclidean distance)
    }
  }
  if (target_space == "pca") { # proximity matrix computed from PCA (principal components x sample)
    message("Construct triangular matrix of sample by sample proximities in PCA space...")
    snp_pca <- RSpectra::svds(t(snp_dat1), k = pc_axes) # single value decomposition for k axes
    snp_pc <- data.frame(snp_pca$u %*% diag(snp_pca$d)) # one column per PCA axis, one row per sample
    colnames(snp_pc) <- paste0("PC", c(1:ncol(snp_pc))) # rename columns with PC + number
    if (program_distance == "Rfast") { # Rfast
      snp_eucli <- Rfast::Dist(snp_pc, method = "euclidean") # default: method = "euclidean" (= Euclidean distance)
    }
    if (program_distance == "vegan") { # vegan
      snp_eucli <- vegan::vegdist(snp_pc, method = "euclidean") # default: method = "euclidean" (= Euclidean distance)
    }
  }

  message("Completed construction of triangular matrix of sample by sample proximities...")

  ##--------------------------------------------------------------##
  # 7. Compute variance partioning by PERMANOVA
  ##--------------------------------------------------------------##

  message(paste0("Computing variance partioning by PERMANOVA: global test..."))

  # Create group variable
  group <- factor(as.matrix(sample_group[sample_PCA]))

  # Set seed of random generator
  set.seed(permutation_seed)

  # Compute PERMANOVA (global test)
  pmanova <- vegan::adonis(formula = snp_eucli ~ group, permutations = permutation_n) # run test
  globalTable.anova <- pmanova$aov.tab[c(1:6)] # extract ANOVA table

  message("Completed PERMANOVA: global test")
  message(paste0("Time elapsed: ", get.time(startT)))

  # Compute PERMANOVA (pairwise tests) / if data contains two groups, pairwise test should mirror global test / Modified after https://rdrr.io/github/GuillemSalazar/EcolUtils/man/adonis.pair.html
  if (pairwise == TRUE & length(levels(group)) > 2) {
    message("Computing variance partioning by PERMANOVA: pairwise tests with/without", pairwise_method, "correction...")
    test.pair <- function(snp_eucli, group, nper = permutation_n, corr.method = pairwise_method) {
      comb.fact <- utils::combn(levels(group), 2)
      F.Model <- NULL
      R2 <- NULL
      pv.anova <- NULL
      for (i in 1:ncol(comb.fact)) { # iterate over group pair-wise comparisons
        message(paste0("Computing pairs ", paste(comb.fact[, i], collapse = " & "), "..."))
        GN <- group %in% comb.fact[, i]
        if(program_distance == "vegan"){
          model.temp <- vegan::adonis(as.matrix(snp_eucli, ncol = sampN)[GN, GN] ~ group[GN], permutations = permutation_n)
        } else {
          model.temp <- vegan::adonis(snp_eucli[GN, GN] ~ group[GN], permutations = permutation_n)
        }
        F.Model <- c(F.Model, model.temp$aov.tab$F.Model[1])
        R2 <- c(R2, model.temp$aov.tab$R2[1])
        pv.anova <- c(pv.anova, model.temp$aov.tab[[6]][1])
      }
      pv.corr.anova <- stats::p.adjust(pv.anova, method = corr.method) # compute FDR-corrected p values
      data.frame(GroupPair = paste(comb.fact[1, ], comb.fact[2, ], sep= "-"), F.Model = F.Model, R2 = R2, P.value = pv.anova, P.value.corrected = pv.corr.anova)
    }
    pairwiseTable <- test.pair(snp_eucli, group, nper = permutation_n) # run pairwise tests and tabulate ANOVA table
    message("Completed PERMANOVA: pairwise tests")
    message(paste0("Time elapsed: ", get.time(startT)))
  }
  if (pairwise == TRUE & length(levels(group)) == 2) {
    message(paste0("Pairwise tests not computed because number of groups is 2"))
  }

  ##--------------------------------------------------------------##
  # 8. Tabulate analytical summary
  ##--------------------------------------------------------------##

  message("Tabulating PERMANOVA outputs...")

  # Compile PERMANOVA results in a list
  sample.out <- merge(samp_dat, data.table::data.table(I = sample_PCA, Class = "PERMANOVA"), by = "I") # summary of samples whether included in PCA, projected and/or removed
  if (!isFALSE(sample_remove)) {
    sample.out <- rbind(sample.out, data.table::data.table(I = sample_remove, Group = sample_group[sample_remove], Class = "Removed"))
  }
  sample.out <- sample.out[order(I)] # retain original sample order
  sample.out[, I:=NULL] #remove index column
  if (pairwise == TRUE & length(levels(group)) > 2) {
    OUT <- list(permanova.samples = as.data.frame(sample.out), permanova.global_test = globalTable.anova, permanova.pairwise_test = pairwiseTable, permanova.pairwise_correction = pairwise_method, permanova.permutation_number = permutation_n, permanova.permutation_seed = permutation_seed) # create list
  } else {
    OUT <- list(permanova.samples = as.data.frame(sample.out), permanova.global_test = globalTable.anova, permanova.pairwise_test = "No pairwise testing implemented", permanova.permutation_number = permutation_n, permanova.permutation_seed = permutation_seed) # create list
  }
  message("Completed tabulations of PERMANOVA outputs")
  return(OUT)
}
# smartsnp v.1
# Coding end date = 08/09/2020
# smartsnp::smart_permanova by Salvador Herrando-Perez (salherra@gmail.com) and Ray Tobler (tingalingx@gmail.com)
