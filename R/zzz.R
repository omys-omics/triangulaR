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
      "Wiens, B. J., DeCicco, L. H., & Colella, J. P. (2025). triangulaR: An R package for identifying AIMs and building triangle plots using SNP data from hybrid zones. Heredity, 1-12.")
  )
}
