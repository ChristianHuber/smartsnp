#include <RcppArmadillo.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

#define PACK_DENSITY 4

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */


// [[Rcpp::export]]
NumericMatrix cpp_read_packedancestrymap(String genofile, int nsnp, int nind, IntegerVector indvec,
                                         int first, int last, bool transpose = false, bool verbose = false) {
  int val;
  long len, bytespersnp;
  int readsnps = last - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcerr << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long)in.tellg();
  bytespersnp = len/(nsnp+1);

  // size of packed data, in bytes, per SNP
  //np = (long)ceil((double)nind / PACK_DENSITY);

  int nindused = 0;
  int* blockused = new int[bytespersnp];
  for(int i = 0 ; i < bytespersnp; i++) {
    blockused[i] = 0;
  }
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] == 1) {
      nindused++;
      blockused[i/PACK_DENSITY] = 1;
    }
  }

  NumericMatrix geno(transpose?nindused:readsnps, transpose?readsnps:nindused);
  std::fill(geno.begin(), geno.end(), NA_REAL);

  // char* header = new char[bytespersnp];
  // in.seekg(0, std::ifstream::beg);
  // in.read((char*)header, bytespersnp);
  // Rcout << "header " << header << std::endl;

  in.seekg((first+1)*bytespersnp, std::ifstream::beg);
  char* tmp = new char[bytespersnp + 1];
  tmp[bytespersnp] = '\0';
  char tmpi;

  // Allocate more than the sample size since data must take up whole bytes
  char* tmp2 = new char[bytespersnp * PACK_DENSITY + 1];
  tmp2[bytespersnp * PACK_DENSITY] = '\0';

  int k;
  for(int j = 0 ; j < readsnps; j++) {
    //for(unsigned int j = 0 ; j < 3; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read((char*)tmp, sizeof(char) * bytespersnp);

    for(int l = 0; l < bytespersnp; l++) {
      if(!blockused[l]) continue;

      tmpi = tmp[l];
      k = PACK_DENSITY * l;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 11
       * which denotes missing.
       */
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
    }

    int c = 0;
    if(!transpose) {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = (double)tmp2[i];
        if(val != 3) geno(j, c) = val;
        c++;
      }
    } else {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = (double)tmp2[i];
        if(val != 3) geno(c, j) = val;
        c++;
      }
    }
  }
  if(verbose) Rcout << std::endl;

  delete[] tmp;
  delete[] tmp2;
  delete[] blockused;
  in.close();

  return geno;
}


