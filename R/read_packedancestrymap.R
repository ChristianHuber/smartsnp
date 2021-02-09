# smartsnp v.1
# Coding start date = 26/01/2021
# smartsnp::read_packedancestrymap by Robert Maier (rmaier@broadinstitute.org), modified by Christian D. Huber (christian.domitian.huber@gmail.com)

# Loading genotype data that is in packedancestrymap format standing for binary or compressed data

#' @name read_packedancestrymap
#'
#' @title Read Files in PACKEDANCESTRYMAP format
#'
#' @description This function loads genotype data in \code{PACKEDANCESTRYMAP} format (binary or compressed).
#'
#' @param pref The prefix of the file name that contains the genotype data (i.e., without the \code{*.geno}).
#' @return Returns a list containing a single element:
#' \itemize{
#'   \item {\code{geno}} {Genotype data as R matrix.}\cr
#' }
#'
#' @export
read_packedancestrymap = function (pref)
{
  pref = normalizePath(pref, mustWork = FALSE)
  fl = paste0(pref, ".geno")
  conn = file(fl, "rb")
  hd = strsplit(readBin(conn, "character", n = 1), " +")[[1]]
  close(conn)
  nindall = as.numeric(hd[2])
  nsnpall = as.numeric(hd[3])
  cat(paste0(basename(pref), ".geno has ", nindall,
             " samples and ", nsnpall, " SNPs.\n"))
  cat(paste0("Reading data for ", nindall, " samples and ",
             nsnpall, " SNPs\n"))
  cat(paste0("Expected size of genotype data: ",
             round((nsnpall * nindall * 8 + nsnpall * 112)/1e+06), " MB\n"))
  indvec = rep(1, nindall)
  geno = cpp_read_packedancestrymap(fl, nsnpall, nindall, indvec,
                                    first = 0, last = nsnpall, transpose = FALSE,
                                    verbose = TRUE)
  outlist = list(geno = geno)
  outlist
}
##### smartsnp v.1
##### Coding end date = 08/02/2021
##### smartsnp::read_packedancestrymap spelling checked by Salvador Herrando-Pérez (salherra@gmail.com)
