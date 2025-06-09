# smartsnp v.1
# Coding start date = 10/09/2020
# smartsnp::smart_mva by Salvador Herrando-Perez (salherra@gmail.com) and Ray Tobler (tingalingx@gmail.com)

# Wrapper integrating functionality of smart_pca, smart_permanova and smart_permdisp
# Principal Component Analysis (PCA), Permutational Multivariate Analysis of Variance (PERMANOVA) and Dispersion (PERMDISP) of genotype data
# Implements scaling by centering, z-scored and SMARTPCA controlling for genetic drift
# Projects ancient samples onto modern PCA space
# Computes global and pairwise MANOVA tests for differences in GROUP LOCATION and DISPERSION for multiple types of proximity measures between samples
# MANOVA alpha values estimated by permutations
# Singular Value Decomposition based on https://genomicsclass.github.io/book/pages/pca_svd.html
# For eigendecomposition see http://people.tamu.edu/~alawing/ESSM689Schedule.html and http://people.tamu.edu/~alawing/materials/ESSM689/pca.pdf
# Original description of SMARTPCA at https://doi.org/10.1371/journal.pgen.0020190 / Patterson, Price and Reich (2006) Population structure and eigenanalysis. PLoS Genetics, 2, e190
# Original description of PERMANOVA at https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x / Anderson (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26, 32-46
# Original description of PERMDISP at https://doi.org/10.1111/j.1461-0248.2006.00926.x / Anderson, Ellingsen & McArdle (2006) Multivariate dispersion as a measure of beta diversity. Ecology Letters, 9, 683-693

#' @name smart_mva
#' @title Smart Multivariate Analyses (wrapper of PCA, PERMANOVA and PERMDISP)
#'
#' @description Computes Principal Component Analysis (PCA) for variable x sample genotype data, such as Single Nucleotide Polymorphisms (SNP), in combination with Permutational Multivariate Analysis of Variance (PERMANOVA) and Permutational Multivariate Analysis of Dispersion (PERMDISP).
#' A wrapper of functions \code{smart_pca}, \code{smart_permanova} and \code{smart_permdisp}.
#' Genetic markers such as SNPs can be scaled by \code{centering}, z-scores and genetic drift-based dispersion.
#' The latter follows the SMARTPCA implementation of Patterson, Price and Reich (2006).
#' Optimized to run fast computation for big datasets.
#'
#' @details See details in other functions for conceptualization of PCA (\code{smart_pca}) (Hotelling 1993), SMARTPCA (Patterson, Price and Reich 2006), PERMANOVA (\code{smart_permanova}) (Anderson 2001) and PERMDISP (\code{smart_permdisp} (Anderson 2006), types of scaling, ancient projection, and correction for multiple testing.\cr
#'
#' Users can compute any combination of the three analyses by assigning \code{TRUE} or \code{FALSE} to \code{pca} and/or \code{permanova} and/or \code{permdisp}.\cr
#'
#' PERMANOVA and PERMDISP exclude samples (columns) specified in either \code{sample_remove} or \code{sample_project}.
#' Projected samples are not used for testing as their PCA coordinates are derived from, and therefore depend on, the coordinates of non-projected samples.\cr
#'
#' Data read from working directory with SNPs as rows and samples as columns. Two alternative formats: (1) text file of SNPs by samples (file extension and column separators recognized automatically) read using \code{\link[data.table]{fread}}; or (2) duet of \code{EIGENSTRAT} files (see \url{https://reich.hms.harvard.edu/software}) using \code{\link[vroom]{vroom_fwf}}, including a genotype file of SNPs by samples (\code{*.geno}), and a sample file (\code{*.ind}) containing three vectors assigning individual samples to unique user-predefined groups (populations), sexes (or other user-defined descriptor) and alphanumeric identifiers.
#' For \code{EIGENSTRAT}, vector \code{sample_group} assigns samples to groups retrievable from column 3 of file \code{*.ind}.
#' SNPs with zero variance removed prior to SVD to optimize computation time and avoid undefined values if \code{scaling = "sd"} or \code{"drift"}.\cr
#'
#' Users can select subsets of samples or SNPs by introducing a vector including column numbers for samples (\code{sample_remove}) and/or row numbers for SNPs (\code{snp_remove}) to be removed from computations.
#' Function stops if the final number of SNPs is 1 or 2.
#' \code{EIGENSOFT} was conceived for the analysis of human genes and its SMARTPCA suite so accepts 22 (autosomal) chromosomes by default.
#' If >22 chromosomes are provided and the internal parameter \code{numchrom} is not set to the target number chromosomes of interest, SMARTPCA automatically subsets chromosomes 1 to 22.
#' In contrast, \code{smart_mva} accepts any number of autosomes with or without the sex chromosomes from an \code{EIGENSTRAT} file.\cr
#'
#' @param snp_data snp_data}{File name read from working directory.
#' SNP = rows, samples = columns without row names or column headings.
#' SNP values must be count data (no decimals allowed).
#' File extension detected automatically whether text or \code{EIGENSTRAT}.
#' See details.
#' @param packed_data Logical value for \code{EIGENSTRAT}, irrelevant for text data.
#' Default \code{packed_data = FALSE} assumes uncompressed \code{EIGENSTRAT}.
#' \code{packed_data = TRUE} for compressed or binary \code{EIGENSTRAT} (\code{PACKENDANCESTRYMAP}).
#' @param sample_group Character or numeric vector assigning samples to groups.
#' Coerced to factor.
#' @param sample_remove Logical \code{FALSE} or numeric vector indicating column numbers (samples) to be removed from computations.
#' Default \code{sample_remove =  FALSE} keeps all samples.
#' @param snp_remove Logical \code{FALSE} or numeric vector indicating row numbers (SNPs) to be removed from computations.
#' Default \code{snp_remove =  FALSE} keeps all SNPs.
#' See details.
#' @param pca Logical indicating if PCA is computed.
#' Default \code{TRUE}.
#' @param permanova Logical indicating if PERMANOVA is computed.
#' Default \code{TRUE}
#' @param permdisp Logical indicating if PERMDISP is computed.
#' Default \code{TRUE}.
#' @param missing_value Number \code{9} or string \code{NA} indicating missing value.
#' Default \code{missing_value = 9} as in \code{EIGENSTRAT}.
#' If no missing values present, no effect on computation.
#' @param missing_impute String handling missing values.
#' Default \code{missing_impute = "mean"} replaces missing values of each SNP by mean of non-missing values across samples.
#' \code{missing_impute = "remove"} removes SNPs with at least one missing value.
#' If no missing values present, no effect on computation.
#' @param scaling String. Default \code{scaling = "drift"} scales SNPs to control for expected allele frequency dispersion caused by genetic drift (SMARTPCA).
#' \code{scaling = "center"} for \code{centering} (covariance-based PCA).
#' \code{scaling = "sd"} for \code{centered} SNPs divided by standard deviation (correlation-based PCA).
#' \code{scaling = "none"} for no scaling.
#' See details.
#' @param program_svd String indicating R package computing single value decomposition (SVD).
#' Default \code{program_svd = "Rspectra"} for \code{\link[RSpectra]{svds}}.
#' \code{program_svd = "bootSVD"} for \code{\link[bootSVD]{fastSVD}}.
#' See details.
#' @param sample_project Numeric vector indicating column numbers (ancient samples) projected onto (modern) PCA space.
#' Default \code{sample_project =  FALSE} implements no projection.
#' See details.
#' @param pc_project Numeric vector indicating the ranks of the PCA axes ancient samples are projected onto. Default \code{pc_ancient = c(1, 2)} for PCA axes 1 and 2.
#' If \code{program_svd = "RSpectra"}, \code{length(pc_ancient)} must be smaller than or equal to \code{pc_axes}.
#' No effect on computation, if no ancient samples present.
#' @param sample_distance Type of inter-sample proximity computed (distance, similarity, dissimilarity).
#' Default is \code{Euclidean distance}.
#' See details.
#' @param program_distance A string value indicating R package to estimate proximities between pairs of samples.
#' Default \code{program_distance = "Rfast"} uses function \code{\link[Rfast]{Dist}}; \code{program_distance = "vegan"} uses \code{\link[vegan]{vegdist}}.
#' See details.
#' @param target_space String.
#' Default \code{target_space = "multidimensional"} applies PERMANOVA and/or PERMDISP to sample-by-sample triangular matrix computed from variable-by-sample data, \code{pc_axes} has no effect on computation. \code{target_space = "pca"} applies PERMANOVA and/or PERMDISP to sample-by-sample data in PCA space, \code{pc_axes} determines number of PCA axes for testing.
#' @param pc_axes Number of PCA axes computed always starting with PCA axis 1.
#' Default \code{pc_axes = 2} computes PCA axes 1 and 2 if \code{target_space = "pca"}.
#' No effect on computation if \code{target_space = "multidimensional"}.
#' @param pairwise Logical.
#' Default \code{pairwise = FALSE} computes global test. \code{pairwise = TRUE} computes global and pairwise tests.
#' @param pairwise_method String specifying type of correction for multiple testing.
#' Default \code{"holm"}.
#' @param permutation_n Number of permutations resulting in PERMANOVA/PERMDISP test \emph{p value}.
#' Default \code{9999}.
#' @param permutation_seed Number fixing random generator of permutations.
#' Default \code{1}.
#' @param dispersion_type String indicating quantification of group dispersion whether relative to spatial \code{"median"} or \code{"centroid"} in PERMDISP.
#' Default \code{"median"}.
#' @param samplesize_bias Logical. \code{samplesize_bias = TRUE} for dispersion weighted by number of samples per group in PERMDISP.
#' Default \code{pairwise = FALSE} for no weighting.
#'
#' @return Returns a list containing the following elements:
#' \itemize{
#' \item{pca.snp_loadings}{Dataframe of principal coefficients of SNPs.
#' One set of coefficients per PCA axis computed.}
#' \item{pca.eigenvalues}{Dataframe of eigenvalues, variance and cumulative variance explained.
#' One eigenvalue per PCA axis computed.}
#' \item{pca_sample_coordinates}{Dataframe showing PCA sample summary. Column \emph{Group} assigns samples to groups. Column \emph{Class} specifies if samples "Removed" from PCA or "Projected" onto PCA space.
#' Sequence of additional columns shows principal components (coordinates) of samples in PCA space (1 column per PCA computed named PC1, PC2, ...).}
#' \item{test_samples}{Dataframe showing test sample summary.
#' Column \emph{Group} assigns samples to tested groups.
#' Column \emph{Class} specifies if samples were used in, or removed from, testing (PERMANOVA and/or PERMDISP).
#' Column \emph{Sample_dispersion} shows dispersion of individual samples relative to spatial \code{"median"} or \code{"centroid"} used in PERMDISP.}
#' \item{permanova.global_test}{List showing PERMANOVA table with degrees of freedom, sum of squares, mean sum of squares, \emph{F} statistic, variance explained (\emph{R2}) and \emph{p} value.}
#' \item{permanova.pairwise_test}{List showing PERMANOVA table with \emph{F} statistic, variance explained (\emph{R2}), \emph{p} value and corrected \emph{p} value per pair of groups.}
#' \item{permdisp.global_test}{List showing PERMDISP table with degrees of freedoms, sum of squares, mean sum of squares, \emph{F} statistic and \emph{p} value.}
#' \item{permdisp.pairwise_test}{List showing PERMDISP table with \emph{F} statistic, \emph{p} value and corrected \emph{p} value per pair of groups.
#' Obtained only if \code{pairwise = TRUE}.}
#' \item{permdisp.bias}{String indicating if PERMDISP dispersion corrected for number of samples per group.}
#' \item{permdisp.group_location}{Dataframe showing coordinates of spatial \code{"median"} or \code{"centroid"} per group in PERMDISP.}
#' \item{test.pairwise_correction}{String indicating type of correction for multiple testing in PERMANOVA and/or PERMDISP.}
#' \item{test.permutation_number}{Number of permutations applied to obtain the distribution of \emph{F} statistic of PERMANOVA and/or PERMDISP.}
#' \item{test.permutation_seed}{Number fixing random generator of permutations of PERMANOVA and/or PERMDISP for reproducibility of results.}
#' }
#'
#' @examples
#' # Path to example genotype matrix "dataSNP"
#' pathToGenoFile = system.file("extdata", "dataSNP", package = "smartsnp")
#'
#' # Assign 50 samples to each of two groups and colors
#' my_groups <- as.factor(c(rep("A", 50), rep("B", 50))); cols = c("red", "blue")
#'
#' # Run PCA, PERMANOVA and PERMDISP
#' mvaR <- smart_mva(snp_data = pathToGenoFile, sample_group = my_groups)
#' mvaR$pca$pca.eigenvalues # extract PCA eigenvalues
#' head(mvaR$pca$pca.snp_loadings) # extract principal coefficients (SNP loadings)
#' head(mvaR$pca$pca.sample_coordinates) # extract PCA PCs (sample position in PCA space)
#'
#' # plot PCA
#' plot(mvaR$pca$pca.sample_coordinates[,c("PC1","PC2")], cex = 2,
#'      pch = 19, col = cols[my_groups], main = "genotype smartpca")
#' legend("topleft", legend = levels(my_groups), cex = 1,
#'        pch = 19, col = cols, text.col = cols)
#'
#' # Extract PERMANOVA table
#' mvaR$test$permanova.global_test
#'
#' # Extract PERMDISP table
#' mvaR$test$permdisp.global_test # extract PERMDISP table
#'
#' # Extract sample summary and dispersion of individual samples used in PERMDISP
#' mvaR$test$test_samples
#'
#' @seealso \code{\link{smart_pca}},
#' \code{\link{smart_permanova}},
#' \code{\link{smart_permdisp}}
#'
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom data.table .SD
#' @importFrom foreach %do%
#' @export
utils::globalVariables(c("proj_i", "S", "i")) # assign non-binding global variables
smart_mva <- function(snp_data, packed_data = FALSE, sample_group,
                      sample_remove = FALSE, snp_remove = FALSE,
                      pca = TRUE, permanova = TRUE, permdisp = TRUE,
                      missing_value = 9, missing_impute = "mean",
                      scaling = "drift",
                      program_svd = "RSpectra",
                      sample_project = FALSE, pc_project = c(1:2), # specific to smart_pca
                      sample_distance = "euclidean", program_distance = "Rfast", # specific to smart_permanova and smart_permdisp
                      target_space = "multidimensional", pc_axes = 2, # specific to smart_permanova and smart_permdisp
                      pairwise = FALSE, pairwise_method = "holm", # specific to smart_permanova and smart_permdisp
                      permutation_n = 9999, permutation_seed = 1, # specific to smart_permanova and smart_permdisp
                      dispersion_type = "median", samplesize_bias = FALSE) # specific to smart_permdisp
{
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
  if(!is.numeric(sample_remove)) sample_remove <- as.logical(sample_remove)
  if(!is.numeric(sample_project)) sample_project <- as.logical(sample_project)

  # Print error and abort analysis if users enter number of sample labels (sample_group) larger or smaller than number of samples in dataset (snp_dat)
  if(isFALSE(sample_project)) pc_project <- 0 # force pc_project to 0 if no sample projection
  if (program_svd == "RSpectra" & max(pc_project) > pc_axes) {
    stop("Dimensionality of projected space (pc_project) must be equal to or smaller than dimensionality of PCA (pc_axes)\nComputation aborted")
  }
  if (missing_impute != "remove" & missing_impute != "mean") { # if users misspel missing_impute
    stop("Check spelling: missing_impute must be 'remove' or 'mean'")
  }
  if (missing_value != 9 & !is.na(missing_value)) { # if users assign missing values not equal to 9 or NA
    stop("Missing values must be coded as 9 or NA")
  }
  if (scaling != "none" & scaling != "center" & scaling != "sd" & scaling != "drift") { # # if users misspel scaling
    stop("Check spelling: scaling must be 'none', 'center', 'sd' or 'drift'")
  }
  if (program_svd != "RSpectra" & program_svd != "bootSVD") { # if users misspel missing_impute
    stop("Check spelling: program_svd must be 'RSpectra' or 'bootSVD'")
  }
  if (target_space != "multidimensional" & target_space != "pca") { # if users misspel target_space
    stop("Check spelling: target_space must be 'multidimensional' or 'pca'")
  }
  message("Argument options are correct...")


  message("Loading data...")

  # Label samples
  samp_dat <- data.table::data.table(sample_group) # tabulate group labels
  data.table::setnames(samp_dat, "Group") # name column vector of group labels
  samp_dat[, I:= 1:.N] # cbind column vectors with samples listed from 1 to to length(sample_group)
  sampN.full <- nrow(samp_dat) # total number of samples in dataset

  # Define samples used for PCA
  sample_test <- 1:sampN.full # assume all samples used as default
  sample_remove_perm <- c() # create empty vector to store ancient and/or removed samples for permanova and permdisp
  if (!isFALSE(sample_remove)) {
    sample_test <- setdiff(sample_test, sample_remove) # remove excluded samples
    sample_project <- setdiff(sample_project, sample_remove) # ensure no projected samples are amongst excluded samples
    sample_remove_perm <- sample_remove # create vector of excluded samples specific to permanova and permdisp functions
  }
  if (!isFALSE(sample_project)) {
    sample_test <- setdiff(sample_test, sample_project) # define samples used for PCA
    sample_remove_perm <- c(sample_project, sample_remove) # create vector specific to permanova and permdisp functions
  }

  # Recalculate number of samples for analyses
  sampN <- length(sample_test)

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
      snp_dat1 <- do.call(cbind, snp_dat[, sample_test])
      if (!isFALSE(sample_project)) {
        snp_dat2 <- do.call(cbind, snp_dat[, sample_project])
      }
    } else { # compressed binary format
      ss <- strsplit(snp_data, "\\.")[[1]]
      snp_data <- paste(ss[-length(ss)], collapse = ".") # remove suffix
      snp_dat <- read_packedancestrymap(pref = snp_data)$geno
      attr(snp_dat, "dimnames") <- NULL # remove row and column name attributes
      snpN.full <- nrow(snp_dat) # number of SNP
      message(paste("Imported", snpN.full, "SNP by", sampN.full, "sample packed eigenstrat genotype matrix (compressed or packed format)"))
      message(paste0("Time elapsed: ", get.time(startT)))
      message("Filtering data...")
      snp_dat1 <- snp_dat[, sample_test]
      if (!isFALSE(sample_project)) {
        snp_dat2 <- snp_dat[, sample_project]
      }
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
    snp_dat1 <- do.call(cbind, snp_dat[, sample_test, with = FALSE])
    if (!isFALSE(sample_project)) {
      snp_dat2 <- do.call(cbind, snp_dat[, sample_project, with = FALSE])
    } else {
      message(paste("No samples projected after PCA computation"))
    }
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
    if (!isFALSE(sample_project)) {
      snp_dat2 <- snp_dat2[snp.keep, ] # subset SNPs by row number across modern samples
    }
  }

  if (snpN.full < 3) {
    stop("Less than 3 SNPs remaining: computation aborted")
  }

  message(paste(snpN.full, "SNPs included in", paste(c("PCA", "PERMANOVA", "PERMDISP")[c(pca, permanova, permdisp)], collapse=" & "), "computations"))
  message(paste(length(snp_remove), "SNPs omitted from PCA computation"))

  # Print number of samples used in different tests
  message(paste(length(sample_test), "samples included in", paste(c("PCA", "PERMANOVA", "PERMDISP")[c(pca, permanova, permdisp)], collapse=" & "), "computations"))
  if (!isFALSE(sample_remove)) {
    if (pca == TRUE) {
      message(paste(length(sample_remove), "samples omitted from PCA computation"))
      if (!isFALSE(sample_project)) {
        message(paste(length(sample_project), "samples projected after PCA computation"))
      }
    }
  }
  if (permanova == TRUE | permdisp == TRUE) {
    message(paste(length(sample_remove_perm), "samples ommitted from", paste(c("PERMANOVA", "PERMDISP")[c(permanova, permdisp)], collapse= " & "), "analyses"))
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
    if (!isFALSE(sample_project)) {
      snp_dat2 <- snp_dat2[keepSNPs, ] # remove invariant SNPs from ancient samples
    }
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
      message("Imputing SNPs with missing values...")
      snp_dat1[missI] <- genoMean[Z]
      message(paste("Imputation with means completed:", length(missI), "missing values imputed"))
    }
    # Removal of SNPs with missing data
    if (missing_impute == "remove") {
      message("Removing SNPs with missing values...")
      snp_dat1 <- snp_dat1[-snpMissI, ] #SNP
      if (!isFALSE(sample_project)) {
        snp_dat2 <- snp_dat2[-snpMissI, ] # ensure snp_dat1 and snp_dat2 have same SNPs in projection
      }
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
      message("Centering data...")
      snp_dat1 <- snp_dat1 - genoMean
    }
    message(paste0("Completed scaling using ", scaling))
    message(paste0("Time elapsed: ", get.time(startT)))
  } else {
    message("Note: SNP-based scaling not used")
  }

  ##--------------------------------------------------------------##
  # 6. Run PCA
  ##--------------------------------------------------------------##

  if (pca == TRUE) {
    ##--------------------------------------------------------------##
    # 6a. Compute singular value decomposition
    ##--------------------------------------------------------------##

    message(paste0("Computing singular value decomposition using ", program_svd, "..."))

    # SVD computation
    if (program_svd == "RSpectra") {
      snp_pca <- RSpectra::svds(t(snp_dat1), k = pc_axes) # single value decomposition for k axes
    } else {
      snp_pca <- bootSVD::fastSVD(t(snp_dat1))
    }

    # Compute total variance in scaled data - necessary to account for % variance explained by each PCA axis
    if (missing_impute == "mean") {
      if (scaling == "standard") { # variance already recomputed
        snp_peigTotal <- sum(genoVar)
      } else {
        snp_peigTotal <- sum(Rfast::rowVars(snp_dat1))
      }
    } else { # no need to recompute variance
      snp_peigTotal <- sum(genoVar)
    }

    message("Completed singular value decomposition")
    message(paste0("Time elapsed: ", get.time(startT)))


    ##--------------------------------------------------------------##
    # 6b. Extract eigenvalues and eigenvectors
    ##--------------------------------------------------------------##

    message(paste0("Extracting eigenvalues and eigenvectors..."))

    # Extract principal components (sample coordinates in ordination) / computation based on https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
    snp_pc <- data.frame(snp_pca$u %*% diag(snp_pca$d)) # one column per PCA axis, one row per sample
    colnames(snp_pc) <- paste0("PC", c(1:ncol(snp_pc))) # rename columns with PC+number expressions

    # Extract principal coefficients (SNP loadings)
    snp_pcoef <- data.frame(snp_pca$v) # one column per PCA axis, one row per SNP
    colnames(snp_pcoef) <- paste0("PC", c(1:ncol(snp_pcoef))) # rename columns with PC+number expressions

    # Extract eigenvalues
    snp_peig <-  snp_pca$d ^ 2 / (sampN - 1) # one eigenvalue per PCA axis
    snp_peigN <- length(snp_peig) # number of eigenvalues

    # Calculate % raw and cumulative variance explained
    snp_eigVar <- snp_peig * 100 / snp_peigTotal # raw
    snp_eigVarCum <- cumsum(snp_eigVar) # cumulative

    # Store eigenvalues and variance explained
    snp_eig <- rbind(snp_peig, snp_eigVar, snp_eigVarCum)
    rownames(snp_eig) <- c("observed eigenvalues", "variance explained", "cumulative variance explained")
    colnames(snp_eig) <- paste("PC", c(1:snp_peigN), sep= "")

    message("Eigenvalues and eigenvectors extracted")
    message(paste0("Time elapsed: ", get.time(startT)))


    ##--------------------------------------------------------------##
    # 6c. Project ancient samples onto modern PCA space
    ##--------------------------------------------------------------##

    if (!isFALSE(sample_project)) {

      message("Projecting ancient samples onto PCA space")
      message("PCA space =", paste("PC", c(pc_project), sep = ""))
      sampN.anc <- ncol(snp_dat2)

      # Scale SNPS in ancient samples (snp_dat2) to SNP scale of variation in modern samples (snp_dat1)
      if (scaling == "drift") { # SMARTPCA scaling controlling for genetic drift following Formula (3) by Patterson et al. 2006 (doi: 10.1371/journal.pgen.0020190)
        snp_dat2 <- (snp_dat2 - genoMean) / snp_drift
        if(!is.na(missing_value)) {
          missing_value <- (missing_value - genoMean) / snp_drift # keep track of missing value
        }
        rm(genoMean, snp_drift)
      }
      if (scaling == "standard") { # z-score scaling whereby all SNPs have mean = 0 and variance = 1 = correlation-based PCA
        snp_dat2 <- (snp_dat2 - genoMean) / snp_sd
        if(!is.na(missing_value)) {
          missing_value <- (missing_value - genoMean) / snp_sd # keep track of missing value
        }
        rm(genoMean, snp_sd)
      }
      if (scaling == "center") {  # only centering applied = covariance-based PCA
        snp_dat2 <- snp_dat2 - genoMean
        if(!is.na(missing_value)) {
          missing_value <- missing_value - genoMean # keep track of missing value
        }
        rm(genoMean)
      }

      # Project ancient samples, loop over each sample
      proj_sp <- snp_pcoef[, pc_project] # get principal coefficients for projected PCA space
      aProj <- foreach::foreach(proj_i = 1:sampN.anc, .final = t, .combine = cbind) %do% {
        p_i <- data.table::data.table(proj_sp, S = snp_dat2[, proj_i]) # take one ancient sample
        if (any(is.na(missing_value))) { # exclude missing data
          NAout <- p_i[!is.na(S)]
        } else {
          NAout <- p_i[S != missing_value]
        }
        x <- as.matrix(NAout[, .SD, .SDcols = 1:(ncol(NAout) - 1)])
        solve(crossprod(x), crossprod(x, NAout$S))
      }
      rm(snp_dat2, proj_sp); gc()

      message("Completed ancient sample projection")
      message(paste(sampN.anc, "ancient samples projected"))
      message(paste0("Time elapsed: ", get.time(startT)))
    }

    ##--------------------------------------------------------------##
    # 6d. Compile PCA results
    ##--------------------------------------------------------------##

    sample.out <- merge(samp_dat, data.table::data.table(I = sample_test, Class = "PCA", snp_pc), by = "I")
    if (!isFALSE(sample_project)) {
      xx <- rep(NA, snp_peigN)
      names(xx) <- paste0("PC", 1:snp_peigN)
      proj_dt <- data.table::data.table(I = sample_project, Group = sample_group[sample_project], Class = "Projected", t(xx))
      for(i in 1:ncol(aProj)){
        data.table::set(proj_dt, j = colnames(aProj)[i], value = aProj[, i])
      }
      sample.out <- rbind(sample.out, proj_dt)
    }
    if (!isFALSE(sample_remove)) {
      xx <- rep(NA, snp_peigN)
      names(xx) <- paste0("PC", 1:snp_peigN)
      sample.out <- rbind(sample.out, data.table::data.table(I = sample_remove, Group = sample_group[sample_remove], Class= "Removed", t(xx)))
    }
    sample.out <- sample.out[order(I)] # retain original sample order
    sample.out[, I:=NULL] #remove index column
    pca.results <- list(pca.snp_loadings = snp_pcoef, pca.eigenvalues = snp_eig, pca.sample_coordinates = as.data.frame(sample.out))

    message("Completed PCA computation")
    message(paste0("Time elapsed: ", get.time(startT)))
  } else {
    pca.results <- NA
  }

  ##--------------------------------------------------------------##
  # 7. Run PERMANOVA and/or PERMDISP
  ##--------------------------------------------------------------##

  if (permanova == TRUE | permdisp == TRUE) {

    ##--------------------------------------------------------------##
    # 7a. Construct triangular inter-sample distance matrix
    ##--------------------------------------------------------------##

    message(paste("Performing", paste(c("PERMANOVA", "PERMDISP")[c(permanova, permdisp)], collapse=" & "), "testing..."))

    # Create group variable
    group <- factor(sample_group[sample_test])

    # Set seed of random generator
    set.seed(permutation_seed)

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
      if (pca == FALSE) { # if pca step used, can recycle SVD matrix
        snp_pca <- RSpectra::svds(t(snp_dat1), k = pc_axes) # single value decomposition for k axes
        snp_pc <- data.frame(snp_pca$u %*% diag(snp_pca$d)) # one column per PCA axis, one row per sample
        colnames(snp_pc) <- paste0("PC", c(1:ncol(snp_pc))) # rename columns with PC + number
      }
      if (program_distance == "Rfast") { # Rfast
        snp_eucli <- Rfast::Dist(snp_pc, method = "euclidean") # default: method = "euclidean" (= Euclidean distance)
      }
      if (program_distance == "vegan") { # vegan
        snp_eucli <- vegan::vegdist(snp_pc, method = "euclidean") # default: method = "euclidean" (= Euclidean distance)
      }
    }

    message("Completed construction of triangular matrix of sample by sample proximities")
    message(paste0("Time elapsed: ", get.time(startT)))

    ##--------------------------------------------------------------##
    # 7b. Compute inter- and intra-group dispersion by PERMANOVA
    ##--------------------------------------------------------------##

    if (permanova == TRUE) {
      message(paste0("Computing variance partioning by PERMANOVA: global test..."))

      # Compute PERMANOVA (global test)
      pmanova <- vegan::adonis2(formula = snp_eucli ~ group, permutations = permutation_n, by = "terms") # run test
      globalTable.anova <- pmanova[c(1:5)] # extract ANOVA table

      message("Completed PERMANOVA: global test")
      message(paste0("Time elapsed: ", get.time(startT)))
    } else {
      globalTable.anova <- "No PERMANOVA test implemented"
    }

    ##--------------------------------------------------------------##
    # 7c. Compute inter- and intra-group dispersion by PERMDISP
    ##--------------------------------------------------------------##

    if(permdisp==TRUE){
      message(paste0("Compute inter- and intra-group dispersion by PERMDISP: global test..."))

      # Compute PERMDISP (global test)
      dispCent <- vegan::betadisper(d = stats::as.dist(snp_eucli), group = group, type = dispersion_type, bias.adjust = samplesize_bias) # estimate group central points (median or centroid)
      globalTable.disp <- vegan::permutest(dispCent, pairwise = pairwise, permutations = permutation_n)[[1]][, -5] # run test and extract to ANOVA table

      # Extract group spatial median (or centroid) and sample distances to spatial median (or dentroid)
      dispGroup <- dispCent$centroids

      # Extract sample distances to spatial median (or dentroid)
      dispSample <- dispCent$distances

      if (dispersion_type == "median") {
        message("Completed PERMDISP: global test based on group spatial medians")
      } else {
        message("Completed PERMDISP: global test based on group centroids")
      }
      message(paste0("Time elapsed: ", get.time(startT)))
    } else {
      globalTable.disp <- "No PERMDISP test implemented"
      dispGroup <- "No PERMDISP test implemented"
      dispSample <- NA
    }

    if (samplesize_bias == TRUE) {
      disp_bias <- "Dispersion adjusted to number of samples per group"
    } else {
      disp_bias <- "Dispersion NOT adjusted to number of samples per group"
    }

    ##--------------------------------------------------------------##
    # 7c. Compute pairwise tests
    ##--------------------------------------------------------------##

    #  if data contains two groups, pairwise test should mirror global test
    #  Modified from pairdwise PERMANOVA after https://rdrr.io/github/GuillemSalazar/EcolUtils/man/adonis.pair.html

    if (pairwise == TRUE) {
      if (length(levels(group)) > 2) {
        message(paste("Computing", paste(c("PERMANOVA", "PERMDISP")[c(permanova, permdisp)], collapse = " & "), "pairwise tests; using", pairwise_method, "method to correct for multiple testing..."))

        test.pair <- function(snp_eucli, group, nper = permutation_n, corr.method = pairwise_method) {
          comb.fact <- utils::combn(levels(group), 2)
          F.Model.anova <- NULL; F.Model.disp <- NULL
          pv.anova <- NULL; pv.disp <- NULL;  R2 <- NULL
          for (i in 1:ncol(comb.fact)) { # iterate over group pair-wise comparisons
            GN <- group %in% comb.fact[, i]
            if (program_distance == "vegan") {
              test.mat <- as.matrix(snp_eucli, ncol=sampN)[GN, GN]
            } else {
              test.mat <- snp_eucli[GN, GN]
            }
            if(permanova == TRUE){
              model.temp1 <- vegan::adonis2(test.mat ~ group[GN], permutations = permutation_n, by = "terms")
              F.Model.anova <- c(F.Model.anova, model.temp1$F[1])
              pv.anova <- c(pv.anova, model.temp1[[5]][1])
              R2 <- c(R2, model.temp1$R2[1])
            }
            if(permdisp == TRUE){
              dispCent_pair <- vegan::betadisper(stats::as.dist(test.mat), group[GN], type = dispersion_type, bias.adjust = samplesize_bias)
              model.temp2 <- vegan::permutest(dispCent_pair, pairwise = FALSE, permutations = permutation_n)
              F.Model.disp <- c(F.Model.disp, model.temp2$tab[[4]][1])
              pv.disp <- c(pv.disp, model.temp2$tab[[6]][1])
            }
          }
          # compute FDR-corrected p values
          if (permanova == TRUE) {
            pv.corr.anova <- stats::p.adjust(pv.anova, method = pairwise_method)
            pair.anova.tab <- data.frame(GroupPair = paste(comb.fact[1, ], comb.fact[2, ], sep= "-"), F.Model = F.Model.anova, R2 = R2, P.value = pv.anova, P.value.corrected = pv.corr.anova)
          } else {
            pair.anova.tab <- "No PERMANOVA pairwise tests implemented"
          }
          if (permdisp == TRUE) {
            pv.corr.disp <- stats::p.adjust(pv.disp, method = pairwise_method)
            pair.disp.tab <- data.frame(GroupPair = paste(comb.fact[1, ], comb.fact[2, ], sep= "-"), F.Model = F.Model.disp, P.value = pv.disp, P.value.corrected = pv.corr.disp)
          } else {
            pair.disp.tab <- "No PERMDISP pairwise tests implemented"
          }
          list(permanova = pair.anova.tab, permdisp = pair.disp.tab)
        }

        pairwiseTable <- test.pair(snp_eucli, group, nper = permutation_n) # run pairwise tests

        message(paste("Completed", paste(c("PERMANOVA", "PERMDISP")[c(permanova, permdisp)], collapse = " & "), "pairwise tests"))
        message(paste0("Time elapsed: ", get.time(startT)))
      } else {
        message(paste0("Pairwise tests not computed because number of groups is 2"))
        pairwiseTable <- list(permdisp = "No PERMDISP pairwise tests implemented", permanova = "No PERMANOVA pairwise tests implemented")
      }
    } else {
      pairwiseTable <- list(permdisp = "No PERMDISP pairwise tests implemented", permanova = "No PERMANOVA pairwise tests implemented")
    }

    ##--------------------------------------------------------------##
    # 7d. Compile PERMANOVA and PERMDISP results
    ##--------------------------------------------------------------##

    # Summarize samples whether included or excluded in PERMANOVA and/or PERMDISP
    cc <- paste(c("PERMANOVA", "PERMDISP")[c(permanova, permdisp)], collapse = "/")
    sample.out <- merge(samp_dat, data.table::data.table(I = sample_test, Class = cc, Sample_dispersion = dispSample), by = "I")
    if (length(sample_remove_perm) > 0) {
      sample.out <- rbind(sample.out, data.table::data.table(I = sample_remove_perm, Group = sample_group[sample_remove_perm], Class = "Removed", Sample_dispersion = NA))
    }
    sample.out <- sample.out[order(I)] # retain original sample order
    sample.out[, I:=NULL] #remove index column

    perm.results <- list(test_samples = as.data.frame(sample.out),
                         permanova.global_test = globalTable.anova, permanova.pairwise_test = pairwiseTable$permanova,
                         permdisp.global_test = globalTable.disp, permdisp.pairwise_test = pairwiseTable$permdisp, permdisp.bias = disp_bias,
                         permdisp.group_location = dispGroup, test.pairwise_correction = pairwise_method,
                         test.permutation_number = permutation_n, test.permutation_seed = permutation_seed)

    message(paste("Completed", paste(c("PERMANOVA", "PERMDISP")[c(permanova, permdisp)], collapse = " & "), "computation"))
    message(paste0("Time elapsed: ", get.time(startT)))
  } else {
    perm.results <- "No PERMANOVA/PERMDISP tests implemented"
  }

  ##--------------------------------------------------------------##
  # 8. Tabulate analytical summary
  ##--------------------------------------------------------------##

  message(paste0("Complete."))
  return(list(data = unname(snp_dat1),  pca = pca.results, test = perm.results))
}
# smartsnp v.1
# Coding end date = 08/09/2020
# smartsnp::smart_mva by Salvador Herrando-Perez (salherra@gmail.com) and Ray Tobler (tingalingx@gmail.com)
