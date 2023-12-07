#' alleleFreqDiff
#'
#' Generate a new vcfR object with only sites with an allele frequency difference above the given threshold
#'
#' @param vcfR vcfR object
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each sample must be assigned to a population
#' @param p1 (character) name of parental population 1
#' @param p2 (character) name of parental population 2
#' @param difference (numeric) allele frequency difference threshold, must be between 0 and 1
#'
#'
#' @return vcfR object
#' @export
#' @importFrom vcfR extract.gt is.het
#'
#'
#' @examples
#' #alleleFreqDiff(vcfR = example.vcfR, pm = example.popmap, p1 = 0, p2 = 20, difference = 0.5)
alleleFreqDiff <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)[-1])), sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  # Extract genotypes
  m <- extract.gt(vcfR)

  # recode alleles
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2

  # Filter and subset the genotypes for the two populations
  p1.gts <- m[, pm[pm$pop == p1,]]$id
  p2.gts <- m[, pm[pm$pop == p2,]]$id

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))

  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)

  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)

  # get names of loci with allele frequency difference above threshold
  loci <- rownames(af[af$diff >= difference,])

  # get names of all loci in vcfR object
  all.loci <- rownames(extract.gt(vcfR))

  # get indices in vcfR object of loci with difference above theshold
  loci.indices <- all.loci %in% loci
  loci.indices <- which(loci.indices)

  # subset vcfR by loci with allele frequency difference in parental pops above threshold
  vcfR.diff <- vcfR[loci.indices]

  # print statement
  s <- nrow(vcfR.diff@gt)
  print(paste0(s, " sites passed allele frequency difference threshold"))
  return(vcfR.diff)
}



#' hybridIndex
#'
#' Calculate hybrid index, heterozygosity, and percent missing data for each sample
#'
#' @param vcfR data in vcfR format
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each sample must be assigned to a population
#' @param p1 (character) name of parental population 1
#' @param p2 (character) name of parental population 2
#'
#' @return dataframe containing hybrid indices, heterozygosities, and percent missing data
#' @export
#'
#' @importFrom vcfR extract.gt is.het
#'
#' @examples
#' #hybridIndex(vcfR = example.vcfR, pm = example.popmap, p1 = 0, p2 = 20)
hybridIndex <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)[-1])), sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  # get number of differences above threshold
  d <- nrow(vcfR@gt)
  print(paste0("calculating hybrid indices and heterozygosities based on ", d, " sites"))

  m <- extract.gt(vcfR)

  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2

  # Filter and subset the genotypes for the two populations
  p1.gts <- m[, pm[pm$pop == p1, "id"]]
  p2.gts <- m[, pm[pm$pop == p2, "id"]]

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))

  # Determine p1 and p2 allele based on allele frequencies
  p1.allele <- ifelse(af_p1 > af_p2, 2, 0)
  p2.allele <- ifelse(af_p2 > af_p1, 2, 0)

  # Create a matrix to store hybrid index scores
  n <- matrix(nrow = nrow(m), ncol = ncol(m))

  # Compare genotypes and assign scores
  n[m == p1.allele] <- 0
  n[m == 1] <- 1
  n[m == p2.allele] <- 2
  n[is.na(m)] <- NA
  n[m == -9] <- NA

  colnames(n) <- colnames(m)
  rownames(n) <- rownames(m)

  # Count alleles and non-missing genotypes for each individual
  counts <- colSums(n, na.rm = TRUE)
  sites <- colSums(!apply(n, MARGIN = 2, is.na))

  # Calculate hybrid index
  hi <- counts / (sites * 2)

  # Create a dataframe for the results
  tri <- data.frame(
    id = names(hi),
    pop = pm[match(names(hi), pm$id), "pop"],
    hybrid.index = hi,
    heterozygosity = colSums(m == 1, na.rm = TRUE) / colSums(!is.na(m)),
    perc.missing = colSums(is.na(m)) / nrow(m)
  )

  return(tri)

}





#' AIMnames
#'
#' Get the names of ancestry informative markers that pass a given allele frequency difference threshold
#'
#' @param vcfR data in vcfR format
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each sample must be assigned to a population
#' @param p1 (character) name of parental population 1
#' @param p2 (character) name of parental population 2
#' @param difference (numeric) allele frequency difference threshold, must be between 0 and 1
#'
#' @return vector of locus names
#' @export
#'
#' @importFrom vcfR extract.gt
#'
#' @examples
#' #AIMnames(vcfR = example.vcfR, pm = example.popmap, p1 = 0, p2 = 20, difference = 1)
AIMnames <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)[-1])), sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  m <- extract.gt(vcfR)

  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2

  # Filter and subset the genotypes for the two populations
  p1.gts <- m[, pm[pm$pop == p1, "id"]]
  p2.gts <- m[, pm[pm$pop == p2, "id"]]

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))

  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)

  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)

  loci <- rownames(af[af$diff >= difference,])

  return(loci)
}






#' specFreqDiff
#'
#' Calculate allele frequency difference for loci above difference threshold. For all allele frequency differences, set difference at 0
#'
#' @param vcfR data in vcfR format
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each sample must be assigned to a population
#' @param p1 (character) name of parental population 1
#' @param p2 (character) name of parental population 2
#' @param difference (numeric) allele frequency difference threshold, must be between 0 and 1
#'
#' @return dataframe containing allele frequencies in parental pops and differences between them
#' @export
#'
#' @importFrom vcfR extract.gt
#'
#' @examples
#' #specFreqDiff(vcfR = example.vcfR, pm = example.popmap, p1 = 0, p2 = 20, difference = 1)
specFreqDiff <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = 0) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)[-1])), sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  m <- extract.gt(vcfR)

  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2

  # Filter and subset the genotypes for the two populations
  p1.gts <- m[, pm[pm$pop == p1, "id"]]
  p2.gts <- m[, pm[pm$pop == p2, "id"]]

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))

  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)

  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)

  # Filter by the difference threshold
  af <- af[af$diff >= difference, ]

  # Print statement
  s <- nrow(af)
  cat(s, "sites passed allele frequency difference threshold\n")

  return(af)
}




#' aimFreqDist
#'
#' Calculate frequencies of AIMs in parental pops. Allele frequency difference threshold must be >= 0.5
#'
#' @param vcfR data in vcfR format
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each sample must be assigned to a population
#' @param p1 (character) name of parental population 1
#' @param p2 (character) name of parental population 2
#' @param difference (numeric) allele frequency difference threshold, must be between 0 and 1
#'
#' @return dataframe containing allele frequency of AIMs in parental pops and differences between them
#' @export
#'
#' @importFrom vcfR extract.gt
#'
#' @examples
#' #aimFreqDist(vcfR = example.vcfR, pm = example.popmap, p1 = 0, p2 = 20, difference = 1)
aimFreqDist <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)[-1])), sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  if (difference < 0.5) {
    stop("Allele frequency difference threshold must be >= 0.5")
  }

  m <- extract.gt(vcfR)

  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2

  # Filter and subset the genotypes for the two populations
  p1.gts <- m[, pm[pm$pop == p1, "id"]]
  p2.gts <- m[, pm[pm$pop == p2, "id"]]

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))

  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)

  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)

  # only keep loci with allele frequency difference above threshold
  af <- af[af$diff >= difference,]

  # polarize allele freqs
  af[af$p1>0.5,"p1"] <- 1 - af[af$p1>0.5,"p1"]
  af[af$p2<0.5,"p2"] <- 1 - af[af$p2<0.5,"p2"]

  # print statement
  s <- nrow(af)
  print(paste0(s, " sites passed allele frequency difference threshold"))
  return(af)
}





#' triangle.plot
#'
#' Generate a triangle plot using the output of hybridIndex
#'
#' @param data Dataframe returned from hybridIndex function
#' @param colors (character) Colors to use for each population. Optional, if not supplied, default colors will be generated
#' @param outline (logical) Whether or not to draw possible triangle space as outline
#' @param cex (character) Size of points
#' @param alpha (numeric) Transparency of points
#' @param jitter (numeric) Amount by which to jitter points on plot (to facilitate visualization)
#'
#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' #triangle.plot(data = your.data, colors = your.colors)
triangle.plot <- function(data = NULL, colors = NULL, outline = T, cex = 2, alpha = 1, jitter = 0) {
  if(is.null(colors)) {
    color_ramp <- colorRampPalette(c("orange", "blue", "green", "red1", "yellow", "purple"))
    colors <- color_ramp(length(unique(data$pop)))
  }
  if(outline) {
    p <- ggplot(data, aes(x=hybrid.index, y=heterozygosity, color=as.factor(pop))) +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter)+
      geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
      geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
      stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab(paste("Hybrid Index"))+
      ylab(paste("Heterozygosity"))+
      labs(title = "") +
      scale_color_manual("pop", values = colors) +
      ylim(c(-0.05,1.05)) +
      xlim(c(-0.05,1.05)) +
      theme_classic()
  }
  if(!outline) {
    p <- ggplot(data, aes(x=hybrid.index, y=heterozygosity, color=as.factor(pop))) +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter)+
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab(paste("Hybrid Index"))+
      ylab(paste("Heterozygosity"))+
      labs(title = "") +
      scale_color_manual("pop", values = colors) +
      ylim(c(-0.05,1.05)) +
      xlim(c(-0.05,1.05)) +
      theme_classic()
  }
  return(p)
}


#' missing.plot
#'
#' Color the triangle plot by perecent missing data in each sample
#'
#' @param data Dataframe returned from hybridIndex function
#' @param outline (logical) Whether or not to draw possible triangle space as outline
#' @param cex (character) Size of points
#' @param alpha (numeric) Transparency of points
#' @param jitter (numeric) Amount by which to jitter points on plot (to facilitate visualization)
#'
#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' #missing.plot(data = your.data)
missing.plot <- function(data = NULL, outline = T, cex = 2, alpha = 1, jitter = 0) {
  if(outline) {
    p <- ggplot(data, aes(x=hybrid.index, y=heterozygosity, color=perc.missing)) +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter)+
      geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
      geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
      stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab(paste("Hybrid Index"))+
      ylab(paste("Heterozygosity"))+
      labs(title = "") +
      ylim(c(-0.05,1.05)) +
      xlim(c(-0.05,1.05)) +
      theme_classic()
  }
  if(!outline) {
    p <- ggplot(data, aes(x=hybrid.index, y=heterozygosity, color=perc.missing)) +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter)+
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab(paste("Hybrid Index"))+
      ylab(paste("Heterozygosity"))+
      labs(title = "") +
      ylim(c(-0.05,1.05)) +
      xlim(c(-0.05,1.05)) +
      theme_classic()
  }
  return(p)
}



