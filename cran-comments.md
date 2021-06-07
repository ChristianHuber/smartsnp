## Test environments
* local OS X installation, R 4.0.2
* fedora linux (devel)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

## Additional CRAN comments (Mar 3, 2021)
Thanks to Gregor Seyer for his constructive review of the first submission. I fixed all issues pointed out. Find below G. Seyer's remarks and my answers to each of them:

Please reduce the length of the title to less than 65 characters.

* We changed the title to "Fast Multivariate Analyses of Big Genomic Data"

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* There are no published references yet.

You write information messages to the console that cannot be easily
suppressed.
It is more R like to generate objects that can be used to extract the
information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning()  or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console.
(except for print, summary, interactive functions)

* I have replaced *cat()* and *print()* calls in all functions with *message()*. I have also turned off any messages to Rcout for the Rcpp function *cpp_read_packedancestrymap()*.

Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

* Fixed.
