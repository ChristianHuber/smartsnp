# smartsnp v.1
# Coding start date = 13/07/2020
# smartsnp::smart_pca by Salvador Herrando-Perez (salherra@gmail.com) and Ray Tobler (tingalingx@gmail.com)

# Principal Component Analysis of genotype data
# Implements scaling by centering, z-scored and SMARTPCA controlling for genetic drift
# Projects ancient samples onto modern PCA space
# Singular Value Decomposition based on https://genomicsclass.github.io/book/pages/pca_svd.html
# For eigendecomposition see http://people.tamu.edu/~alawing/ESSM689Schedule.html and http://people.tamu.edu/~alawing/materials/ESSM689/pca.pdf
# Original description of SMARTPCA at https://doi.org/10.1371/journal.pgen.0020190 / Patterson, Price and Reich (2006) Population structure and eigenanalysis. PLoS Genetics, 2, e190

#' @name smart_pca
#'
#' @title Smart Principal Component Analysis
#'
#' @description Compute Principal Component Analysis (PCA) for variable x sample genotype data including covariance (\code{centered}), correlation (z-score) and SMARTPCA scaling,
#' and implements projection of ancient samples onto modern PCA space. SMARTPCA scaling controls for genetic drift when variables are bi-allelic genetic markers
#' such as single nucleotide polymorphisms (SNP) following Patterson, Price and Reich (2006).
#' Optimized to run fast single value decomposition for big datasets.
#'
#' @details PCA is a rigid rotation of a Cartesian coordinate system (samples = points, axes = variables or SNPs) that maximizes the dispersion of points along a new system of axes (Pearson 1901; Hotelling 1933; Jolliffe 2002).
#' In rotated space (ordination), axes are \code{principal axes} (PCA axes), \code{eigenvalues} measure variance explained, and \code{principal coefficients} measure importance of SNPs (eigenvectors), \code{principal components} are coordinates of samples (i.e., linear combinations of scaled variables weighted by eigenvectors).
#' Principal coefficients are direction cosines between original and PCA axes (Legendre & Legendre 2012). PCA can be computed by \code{eigenanalysis} or, as implemented here, single value decomposition (SVD). \cr
#'
#' SNPs can be scaled in four different ways prior to SVD: (1) no scaling; (2) covariance: SNPs \code{centered} such that \emph{M(i,j)} = \emph{C(i,j)} minus \emph{mean(j)}) where \emph{C(i,j)} is the number of variant alleles for SNP \emph{j} and sample \emph{i}, and \emph{M(i,j)} is the \code{centered} value of each data point; (3) correlation (z-scores): SNPs \code{centered} then divided by standard deviation \emph{sd(j)}, (4) SMARTPCA: SNPs \code{centered} then divided by \emph{sqrt(p(j)(1-p(j)))}, where \emph{p(j)} equals \emph{mean(j)} divided by \emph{2}, quantifies the underlying allele frequency (autosomal chromosomes) and conceptualizes that SNP frequency changes at rate proportional to \emph{sqrt(p(j)(1-p(j)))} per generation due to genetic drift (Patterson, Price and Reich 2006).
#' SMARTPCA standardization results in all SNPs that comply with Hardy-Weinberg equilibrium having identical variance.
#' SMARTPCA (Patterson, Price and Reich 2006) and \code{EIGENSTRAT} (Price, Patterson, Plenge, Weinblatt, Shadick and Reich 2006) are the computing suites of software \code{EIGENSOFT} (\url{https://reich.hms.harvard.edu/software}).\cr
#'
#' \code{\link[RSpectra]{svds}} runs single value decomposition much faster than \code{\link[bootSVD]{fastSVD}}. With \code{\link[RSpectra]{svds}}, \code{pc_axes} indicates number of eigenvalues and eigenvectors computed starting from PCA axis 1. \code{\link[bootSVD]{fastSVD}} computes all eigenvalues and eigenvectors. Eigenvalues calculated from singular values divided by number of samples minus 1. If number of samples equals number of SNPS, \code{\link[bootSVD]{fastSVD}} prints message alert that no computing efficiency is achieved for square matrices.\cr
#'
#' Ancient samples (with many missing values) can be projected onto modern PCA space derived from modern samples.
#' Following Nelson Taylor and MacGregor (1996), the projected coordinates of a given ancient sample equal the slope coefficient of linear fit through the origin of (scaled) non-missing SNP values of that sample (response) versus principal coefficients of same SNPs in modern samples.
#' Number of projected coordinates per ancient sample given by \code{length(pc_ancient)}.
#' With \code{\link[RSpectra]{svds}}, \code{pc_axes} must be larger or equal to \code{length(pc_ancient)}.\cr
#'
#' Data read from working directory with SNPs as rows and samples as columns.
#' Two alternative formats: (1) text file of SNPs by samples (file extension and column separators recognized automatically) read using \code{\link[data.table]{fread}}; or (2) duet of \code{EIGENSTRAT} files (see \url{https://reich.hms.harvard.edu/software}) using \code{\link[vroom]{vroom_fwf}}, including a genotype file of SNPs by samples (\code{*.geno}), and a sample file (\code{*.ind}) containing three vectors assigning individual samples to unique user-predefined groups (populations), sexes (or other user-defined descriptor) and alphanumeric identifiers.
#' For \code{EIGENSTRAT}, vector \code{sample_group} assigns samples to groups retrievable from column of file \code{*.ind}. SNPs with zero variance removed prior to SVD to optimize computation time and avoid undefined values if \code{scaling = "sd"} or \code{"drift"}.\cr
#'
#' Users can select subsets of samples or SNPs by introducing a vector including column numbers for samples (\code{sample_remove}) and/or row numbers for SNPs (\code{snp_remove}) to be removed from computations.
#' Function stops if the final number of SNPs is 1 or 2.
#' \code{EIGENSOFT} was conceived for the analysis of human genes and its SMARTPCA suite so accepts 22 (autosomal) chromosomes by default.
#' If >22 chromosomes are provided and the internal parameter \code{numchrom} is not set to the target number chromosomes of interest, SMARTPCA automatically subsets chromosomes 1 to 22.
#' In contrast, \code{smart_pca} accepts any number of autosomes with or without the sex chromosomes from an \code{EIGENSTRAT} file.\cr
#'
#' @param snp_data File name read from working directory.
#' SNP = rows, samples = columns without row names or column headings.
#' SNP values must be count data (no decimals allowed). File extension detected automatically whether text or \code{EIGENSTRAT}.
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
#' @param pc_axes A numeric value.
#' If \code{program_svd = "Rspectra"} this argument indicates number of PCA axes computed starting with PCA axis 1.
#' Default \code{pc_axes = 2} computes PCA axes 1 and 2.
#' No effect on computation if \code{program_svd = "bootSVD"} since all PCA axes are computed.
#' @param sample_project Numeric vector indicating column numbers (ancient samples) projected onto (modern) PCA space.
#' Default \code{sample_project =  FALSE} indicates no samples will be used for projection.
#' See details.
#' @param pc_project Numeric vector indicating the ranks of the PCA axes ancient samples are projected onto.
#' Default \code{pc_ancient = c(1, 2)} for PCA axes 1 and 2. If \code{program_svd = "RSpectra"}, \code{length(pc_ancient)} must be smaller than or equal to \code{pc_axes}.
#' No effect on computation, if no ancient samples present.
#'
#' @return Returns a list containing the following elements:
#' \itemize{
#'   \item {\code{pca.snp_loadings}} {Dataframe of principal coefficients of SNPs. One set of coefficients per PCA axis computed.}\cr
#'   \item {\code{pca.eigenvalues}} {Dataframe of eigenvalues, variance and cumulative variance explained. One eigenvalue per PCA axis computed.}\cr
#'   \item {\code{pca_sample_coordinates}} {Dataframe showing PCA sample summary.
#'   Column \emph{Group} assigns samples to groups.
#'   Column \emph{Class} specifies if samples "Removed" from PCA or "Projected" onto PCA space.
#'   Sequence of additional columns shows principal components (coordinates) of samples in PCA space (1 column per PCA computed named PC1, PC2, ...).}
#' }
#'
#' @examples
#' # Path to example genotype matrix "dataSNP"
#' pathToGenoFile = system.file("extdata", "dataSNP", package = "smartsnp")
#'
#' # Example 1: modern samples
#' #assign 50 samples to each of two groups and colors
#' my_groups <- c(rep("A", 50), rep("B", 50)); cols = c("red", "blue")
#' #run PCA with truncated SVD (PCA 1 x PCA 2)
#' pcaR1 <- smart_pca(snp_data = pathToGenoFile, sample_group = my_groups)
#' pcaR1$pca.eigenvalues # extract eigenvalues
#' head(pcaR1$pca.snp_loadings) # extract principal coefficients (SNP loadings)
#' head(pcaR1$pca.sample_coordinates) # extract principal components (sample position in PCA space)
#' #plot PCA
#' plot(pcaR1$pca.sample_coordinates[,c("PC1","PC2")], cex = 2,
#'      pch = 19, col = cols[as.factor(my_groups)], main = "genotype smartpca")
#' legend("topleft", legend = levels(as.factor(my_groups)), cex =1,
#'      pch = 19, col = cols, text.col = cols)
#'
#' # Example 2: modern and ancient samples (ancient samples projected onto modern PCA space)
#' #assign samples 1st to 10th per group to ancient
#' my_ancient <- c(1:10, 51:60)
#' #run PCA with truncated SVD (PCA 1 x PCA 2)
#' pcaR2 <- smart_pca(snp_data = pathToGenoFile, sample_group = my_groups, sample_project = my_ancient)
#' pcaR2$pca.eigenvalues # extract eigenvalues
#' head(pcaR2$pca.snp_loadings) # extract principal coefficients (SNP loading)
#' head(pcaR2$pca.sample_coordinates) # extract principal components (sample position in PCA space)
#' #assign samples to groups (A, ancient, B) and colors
#' my_groups[my_ancient] <- "ancient"; cols = c("red", "black", "blue")
#' #plot PCA
#' plot(pcaR2$pca.sample_coordinates[,c("PC1","PC2")],
#'      cex = 2, col = cols[as.factor(my_groups)], pch = 19, main = "genotype smartpca")
#' legend("topleft", legend = levels(as.factor(my_groups)),  cex = 1,
#'      pch = 19, col = cols, text.col = cols)
#'
#' @references
#' Hotelling, H. (1933) Analysis of a complex of statistical variables into principal components. Journal of Educational Psychology, 24, 417-441.\cr
#'
#' Jolliffe, I.T. (2002) Principal Component Analysis (Springer, New York, USA).\cr
#'
#' Legendre, P. & L. F. J. Legendre (2012). Numerical ecology.  Developments in environmental modelling (Elsevier, Oxford, UK).\cr
#'
#' Nelson, P.R.C., P.A. Taylor, and J.F. MacGregor (1996) Missing data methods in PCA and PLS: score calculations with incomplete observations. Chemometrics and Intelligent Laboratory Systems, 35, 45-65.\cr
#'
#' Patterson, N.J., A. L. Price and D. Reich (2006) Population structure and eigenanalysis. PLoS Genetics, 2, e190.\cr
#'
#' Pearson, K. (1901) On lines and planes of closest fit to systems of points in space. Philosophical Magazine, 2, 559-572.\cr
#'
#' Price, A.L., N.J. Patterson, R.M. Plenge, M.E. Weinblatt, N.A. Shadick and David Reich (2006). Principal components analysis corrects for stratification in genome-wide association studies. Nature Genetics, 38, 904-909.
#'
#' @seealso \code{\link[bootSVD]{fastSVD}} (package \bold{bootSVD}),
#' \code{\link[foreach]{foreach}} (package \bold{foreach}),
#' \code{\link[data.table]{fread}} (package \bold{data.table}),
#' \code{\link[Rfast]{rowVars}} (package \bold{Rfast}),
#' \code{\link[RSpectra]{svds}} (package \bold{RSpectra}),
#' \code{\link[vroom]{vroom_fwf}} (package \bold{vroom})
#'
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom data.table .SD
#' @importFrom foreach %do%
#'
#' @export
utils::globalVariables(c("proj_i", "S", "i")) # assign non-binding global variables
smart_pca <- function(snp_data, packed_data = FALSE,
                      sample_group, sample_remove = FALSE, snp_remove = FALSE,
                      missing_value = 9, missing_impute = "mean",
                      scaling = "drift",
                      program_svd = "RSpectra", pc_axes = 2,
                      sample_project = FALSE, pc_project = c(1:2)) {

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
  # Print error and abort analysis if users misspel missing_impute
  if (missing_impute != "remove" & missing_impute != "mean") {
    stop("Check spelling: missing_impute must equal 'remove' or 'mean'")
  }
  message("Argument options are correct...")

  message("Loading data...")

  # Label samples
  samp_dat <- data.table::data.table(sample_group) # tabulate group labels
  data.table::setnames(samp_dat, "Group") # name column vector of group labels
  samp_dat[, I:= 1:.N] # cbind column vectors with samples listed from 1 to length(sample_group)
  sampN.full <- nrow(samp_dat) # total number of samples in dataset

  # Define samples used for PCA
  sample_PCA <- 1:sampN.full # assume all samples for PCA as default
  if (!isFALSE(sample_remove)) {
    sample_PCA <- setdiff(sample_PCA, sample_remove) # remove excluded samples
    sample_project <- setdiff(sample_project, sample_remove) # ensure no projected samples are amongst excluded samples
  }
  if (!isFALSE(sample_project)) {
    sample_PCA <- setdiff(sample_PCA, sample_project) # define samples used for PCA
  }
  # Recalculate number of samples for PCA
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
      snp_dat1 <- snp_dat[, sample_PCA]
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
    snp_dat1 <- do.call(cbind, snp_dat[, sample_PCA, with = FALSE])
    if (!isFALSE(sample_project)) {
      snp_dat2 <- do.call(cbind, snp_dat[, sample_project, with=FALSE])
    } else {
      message(paste("No samples projected after PCA computation"))
    }
  }

  # Print error if users enter number of sample labels (sample_group) larger or smaller than number of samples in dataset (snp_dat)
  if (length(sample_group) != ncol(snp_dat)) {
    n_str <- paste0("n dataset = ",ncol(snp_dat),", but length(sample_group) is ",length(sample_group),"!")
    stop(n_str," length(sample_group) should be equal to number of samples in dataset: computation aborted")
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

  message(paste(snpN.full, "SNPs included in PCA", "computation"))
  if (length(snp_remove) != 1 | !is.logical(snp_remove)) message(paste(length(snp_remove), "SNPs omitted from PCA computation"))

  # Print number of samples used in, and projected onto, PCA
  message(paste(length(sample_PCA), "samples included in PCA computation"))
  if (!isFALSE(sample_remove)) {
    message(paste(length(sample_remove), "samples omitted from PCA computation"))
  }
  if (!isFALSE(sample_project)) {
    message(paste(length(sample_project), "samples projected after PCA computation"))
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
      snp_dat2 <- snp_dat2[keepSNPs, , drop = FALSE] # remove invariant SNPs from ancient samples
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
  # 6. Compute singular value decomposition
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
    if (scaling == "sd") { # variance already recomputed
      snp_peigTotal <- sum(genoVar)
    } else {
      snp_peigTotal <- sum(Rfast::rowVars(snp_dat1))
    }
  } else { # no need to recompute variance
    snp_peigTotal <- sum(genoVar)
  }

  # Remove PCA data matrix
  rm(snp_dat1); gc()

  message(paste0("Completed singular value decomposition using ", program_svd))
  message(paste0("Time elapsed: ", get.time(startT)))

  ##--------------------------------------------------------------##
  # 7. Extract eigenvalues and eigenvectors
  ##--------------------------------------------------------------##

  message(paste0("Extracting eigenvalues and eigenvectors..."))

  # Extract principal components (sample coordinates in ordination) / computation based on https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
  snp_pc <- data.frame(snp_pca$u %*% diag(snp_pca$d)) # one column per PCA axis, one row per sample
  colnames(snp_pc) <- paste0("PC", c(1:ncol(snp_pc))) # rename columns with PC + number

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

  rm(snp_pca)
  message("Eigenvalues and eigenvectors extracted")
  message(paste0("Time elapsed: ", get.time(startT)))

  ##--------------------------------------------------------------##
  # 8. Project ancient samples onto modern PCA space
  ##--------------------------------------------------------------##

  if (!isFALSE(sample_project)) {

    message("Projecting ancient samples onto PCA space")
    message("PCA space = ", paste("PC", c(pc_project), sep = ""))
    sampN.anc <- ncol(snp_dat2)

    # Scale SNPS in ancient samples (snp_dat2) to SNP scale of variation in modern samples (snp_dat1)
    if (scaling == "drift") { # SMARTPCA scaling controlling for genetic drift following Formula (3) by Patterson et al. 2006 (doi: 10.1371/journal.pgen.0020190)
      snp_dat2 <- (snp_dat2 - genoMean) / snp_drift
      if(!is.na(missing_value)) {
        missing_value <- (missing_value - genoMean) / snp_drift # keep track of missing value
      }
      rm(genoMean, snp_drift)
    }
    if (scaling == "sd") { # z-score scaling whereby all SNPs have mean = 0 and variance = 1 = correlation-based PCA
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
  # 9. Tabulate analytical summary
  ##--------------------------------------------------------------##

  message("Tabulating PCA outputs...")

  sample.out <- merge(samp_dat, data.table::data.table(I = sample_PCA, Class = "PCA", snp_pc), by = "I") # summary of samples whether included in PCA, projected and/or removed
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
    sample.out <- rbind(sample.out, data.table::data.table(I = sample_remove, Group = sample_group[sample_remove], Class = "Removed", t(xx)))
  }
  sample.out <- sample.out[order(I)] # retain original sample order
  sample.out[, I:=NULL] #remove index column
  OUT <- list(pca.snp_loadings = snp_pcoef, pca.eigenvalues = snp_eig, pca.sample_coordinates = as.data.frame(sample.out)) # output SNPS loadings, eigenvalues and sample coordinates in PCA space
  message("Completed tabulations of PCA outputs...")
  return(OUT)
}
##### smartsnp v.1
##### Coding end date = 08/09/2020
##### smartsnp::smart_pca by Salvador Herrando-Perez (salherra@gmail.com) and Ray Tobler (tingalingx@gmail.com)
