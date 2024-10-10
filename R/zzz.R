#print package startup message
.onAttach <- function(libname, pkgname) {
  pkg.version <- utils::packageVersion("triangulaR")
  packageStartupMessage(
    c(paste0("This is triangulaR v.", pkg.version, "\n\n"),
      "          /\\
         /  \\
        /    \\
       /______\\", "\n\n",
      "Usage information is available at: https://github.com/omys-omics/triangulaR/ \n\n",
      "Please cite the following if you use triangulaR in a publication: \n\n",
      "Wiens, B. J., & Colella, J. P. (2024). triangulaR: an R package for identifying AIMs and building triangle plots using SNP data from hybrid zones. bioRxiv, 2024-03.")
  )
}
