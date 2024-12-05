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
  p1.gts <- m[, pm[pm$pop == p1,]$id]
  p2.gts <- m[, pm[pm$pop == p2,]$id]

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
  p1.gts <- m[, pm[pm$pop == p1,]$id]
  p2.gts <- m[, pm[pm$pop == p2,]$id]

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
  p1.gts <- m[, pm[pm$pop == p1,]$id]
  p2.gts <- m[, pm[pm$pop == p2,]$id]

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
  p1.gts <- m[, pm[pm$pop == p1,]$id]
  p2.gts <- m[, pm[pm$pop == p2,]$id]

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
  p1.gts <- m[, pm[pm$pop == p1,]$id]
  p2.gts <- m[, pm[pm$pop == p2,]$id]

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
#' @param ind.labels (logical) Whether or not to label each individual on the triangle plot
#' @param cex (character) Size of points
#' @param alpha (numeric) Transparency of points
#' @param jitter (numeric) Amount by which to jitter points on plot (to facilitate visualization)
#' @param max.overlaps (numeric) Only necessary if labeling individuals. Increasing this will increase the number of individuals labeled, even if labels overlap.
#'
#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom ggrepel geom_label_repel
#'
#' @examples
#' #triangle.plot(data = your.data, colors = your.colors)
triangle.plot <- function(data = NULL, colors = NULL, outline = T, ind.labels = F, cex = 2, alpha = 1, jitter = 0, max.overlaps = 10) {
  if(is.null(colors)) {
    color_ramp <- colorRampPalette(c("orange", "blue", "green", "red1", "yellow", "purple"))
    colors <- color_ramp(length(unique(data$pop)))
  }
  if(outline) {
    p <- ggplot(data, aes(x=hybrid.index, y=heterozygosity, color=as.factor(pop))) +
      geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
      geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
      stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter)+
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab(paste("Hybrid Index"))+
      ylab(paste("Interclass Heterozygosity"))+
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
      ylab(paste("Interclass Heterozygosity"))+
      labs(title = "") +
      scale_color_manual("pop", values = colors) +
      ylim(c(-0.05,1.05)) +
      xlim(c(-0.05,1.05)) +
      theme_classic()
  }
  if(ind.labels) {
    p <- p + geom_label_repel(aes(label=id), size=2, max.overlaps = max.overlaps)
  }
  return(p)
}


#' missing.plot
#'
#' Color the triangle plot by percent missing data in each sample
#'
#' @param data Dataframe returned from hybridIndex function
#' @param outline (logical) Whether or not to draw possible triangle space as outline
#' @param ind.labels (logical) Whether or not to label each individual on the triangle plot
#' @param cex (character) Size of points
#' @param alpha (numeric) Transparency of points
#' @param jitter (numeric) Amount by which to jitter points on plot (to facilitate visualization)
#' @param max.overlaps (numeric) Only necessary if labeling individuals. Increasing this will increase the number of individuals labeled, even if labels overlap.

#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#'
#' @examples
#' #missing.plot(data = your.data)
missing.plot <- function(data = NULL, outline = T, ind.labels = F, cex = 2, alpha = 1, jitter = 0, max.overlaps = 10) {
  if(outline) {
    p <- ggplot(data, aes(x=hybrid.index, y=heterozygosity, color=perc.missing)) +
      geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
      geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
      stat_function(fun = function(hi) 2*hi*(1-hi), xlim = c(0,1), color = "black", linetype = "dashed") +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter)+
      guides(shape = guide_legend(override.aes = list(size = 5), order=2, label.theme= element_text(face="italic")))+
      xlab(paste("Hybrid Index"))+
      ylab(paste("Interclass Heterozygosity"))+
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
      ylab(paste("Interclass Heterozygosity"))+
      labs(title = "") +
      ylim(c(-0.05,1.05)) +
      xlim(c(-0.05,1.05)) +
      theme_classic()
  }
  if(ind.labels) {
    p <- p + geom_label_repel(aes(label=id), size=2, max.overlaps = max.overlaps)
  }
  return(p)
}


#' read.genotypes
#'
#' Read in genotype data from a tab-delimited file in matrix format (e.g. structure format) and convert it to a vcfR object
#'
#' @param file Name of the tab-delimited file. Individuals should be in rows and alleles/genotypes should be in columns.
#' @param data.type (character) Either "alleles" or "genotypes". Use alleles if each allele (e.g. A,C,G,T OR 1,2,3,4) is encoded separately in the file, either across rows or columns. Use genotypes if the alleles are encoded as a single genotype (e.g. 0=homozygous ref, 1=heterozygous, 2=homozygous alt). Default is alleles.
#' @param missing.value (character) Value used for missing data. Default is "N".
#' @param n.id.cols (numeric) Number of columns before genotypes start. Must be at least 1 with individual IDs.
#' @param n.id.rows (numeric) Number of rows before genotypes start. Usually there is 1 with SNP IDs. If not, artificial IDs will be assigned.
#' @param second.allele (numeric) Either "columns" or "rows". Use columns if second allele is in the next column. Use rows if second allele is in the next row. Not required if data are in genotype format.

#' @return vcfR object
#' @export
#'
#' @import vcfR
#'
#' @examples
#' #read.genotypes(file = alleles.str, data.type = "alleles", missing.value = -9, n.id.cols = 2, n.id.rows = 1, second.allele = "rows")
#' #read.genotypes(file = genotypes.txt, data.type = "genotypes", missing.value = "N", n.id.cols = 2, n.id.rows = 1)
read.genotypes <- function(file = NULL, data.type = "alleles", missing.value = "N", n.id.cols = NULL, n.id.rows = NULL, second.allele = "NA") {
  if (is.null(n.id.cols)) {
    stop( "Please indicate how many id columns there are. There should be at least one containing individual IDs")
  }
  if (is.null(n.id.rows)) {
    stop( "Please indicate how many id rows there are. There could be one containing snp IDs, but if not, indicate 0")
  }
  if (data.type == "alleles") {
    if (second.allele == "NA") {
      stop( "Please indicate whether the second allele occurs in the following column or the following row")
    }
  }

  if (n.id.cols < 1) {
    stop("There must be at least 1 column with individual IDs.")
  }

  # read in data
  g <- read.table(file = file, sep ="\t", header = F, colClasses = "character")
  g <- as.matrix(g)

  # assign SNP ids
  if (n.id.rows > 0) {
    ID <- unique(g[1,c(-1:(n.id.cols*-1))]) # retrieve SNP ids if they are included
  } else {
    if (second.allele == "columns") {
      ID <- 1:((ncol(g) - n.id.cols)/2) # assign sequentially by number of columns with genotypes, divided by 2
    }
    if (second.allele == "rows" || second.allele == "NA") {
      ID <- 1:(ncol(g) - n.id.cols) # assign sequentially by number of columns with genotypes
    }
  }

  # assign individual ids
  if (n.id.rows > 0) {
    inds <- unique(g[-1:(n.id.rows*-1),1])
  } else {
    inds <- unique(g[,1])
  }

  # remove ids from beginning rows and columns
  if (n.id.rows > 0) {
    g <- g[-1:(n.id.rows*-1),-1:(n.id.cols*-1)]
  } else {
    g <- g[,-1:(n.id.cols*-1)]
  }
  # set ids as column and row names
  if (data.type == "genotypes") {
    colnames(g) <- ID
    rownames(g) <- inds
    suppressWarnings(class(g) <- "numeric")
  } else {
    if (second.allele == "columns") {
      ID2 <- c()
      for (id in ID) {
        ID2 <- c(ID2, paste(id, 1, sep = "."))
        ID2 <- c(ID2, paste(id, 2, sep = "."))
      }
      colnames(g) <- ID2
    }
    if (second.allele == "rows") {
      colnames(g) <- ID
    }

    if (second.allele == "columns") {
      rownames(g) <- inds
    }
    if (second.allele == "rows") {
      inds2 <- c()
      for (ind in inds) {
        inds2 <- c(inds2, paste(ind, 1, sep = "."))
        inds2 <- c(inds2, paste(ind, 2, sep = "."))
      }
      rownames(g) <- inds2
    }
  }

  # Pivot matrix so second allele occurs in the following column instead of in the following row
  if (second.allele == "rows") {
    h <- matrix(nrow = nrow(g)/2, ncol = ncol(g)*2)
    for (a in 1:nrow(h)) {
      b <- a*2
      ig <- c()
      for (d in 1:ncol(g)) {
        ig <- c(ig, g[(b-1):b,d])
      }
      h[a,] <- ig
    }

    # rename rows with individual ids
    rownames(h) <- inds

    # rename columns with SNP ids
    ID2 <- c()
    for (id in ID) {
      ID2 <- c(ID2, paste(id, 1, sep = "."))
      ID2 <- c(ID2, paste(id, 2, sep = "."))
    }
    colnames(h) <- ID2
    g <- h
  }

  # Build vcfR fields if data is in allele format with either letters or numbers (e.g. A,C,G,T OR 1,2,3,4)
  if (data.type == "alleles") {
    # setup fix dataframe
    fix <- matrix(ncol=8, nrow=0)
    colnames(fix) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

    # setup gt dataframe
    gt <- matrix(ncol=nrow(g)+1, nrow = 0)
    colnames(gt) <- c("FORMAT", rownames(g))

    # Placeholders for each SNP's VCF metadata
    chrom <- "NA"  # Chromosome placeholder
    pos <- "NA"    # Position placeholder (you can modify this if you have position data)
    qual <- "."     # Quality placeholder
    filter <- "PASS"  # Placeholder for filter
    info <- "."     # Placeholder for info
    format <- "GT"  # Format for genotype

    # initialize counter for multiallelic sites
    multiallelic.sites <- 0

    for (snp_idx in 1:(ncol(g)/2)) {
      # SNP id for metadata
      snp_id <- sapply(strsplit(colnames(g)[2*snp_idx], "\\."), '[', 1)  # Unique SNP ID

      # Extract alleles for the current SNP across all samples (alleles are in columns 2*snp_idx-1 and 2*snp_idx)
      allele1 <- as.character(g[, (2 * snp_idx - 1)])  # Allele 1
      allele2 <- as.character(g[, (2 * snp_idx)])      # Allele 2

      # Combine the two allele columns for all samples and exclude "N"
      alleles <- c(allele1, allele2)
      alleles <- alleles[alleles != missing.value]  # Remove "N" (missing data)

      # Count allele frequencies
      allele_freq <- table(alleles)

      # Assign REF and ALT based on the frequency of alleles
      if (length(allele_freq) > 1) {
        if (length(allele_freq) == 2 && allele_freq[1] == allele_freq[2]) {
          # If frequencies are the same, randomly choose REF and ALT
          ref <- sample(names(allele_freq), 1)
          alt <- setdiff(names(allele_freq), ref)
        } else {
          # Otherwise, pick most frequent as REF and least frequent as ALT
          ref <- names(allele_freq)[which.max(allele_freq)]  # Most frequent allele
          alt <- names(allele_freq)[which.min(allele_freq)]  # Least frequent allele
        }
      } else {
        # If there's only one unique allele (monomorphic site)
        ref <- names(allele_freq)
        alt <- "N"  # No alternative allele for monomorphic sites
      }

      # Extract genotypes for each sample for this SNP
      genotypes <- sapply(1:nrow(g), function(i) {
        allele1 <- as.character(g[i, (2 * snp_idx - 1)])  # First allele
        allele2 <- as.character(g[i, (2 * snp_idx)])      # Second allele

        # Check for missing data or "N" alleles
        if (is.na(allele1) || is.na(allele2) || allele1 == missing.value || allele2 == missing.value) {
          NA  # Missing data format
        } else {
          # Convert genotypes to 0/1, 1/1, etc. based on REF/ALT
          if (allele1 == ref && allele2 == ref) {
            "0/0"
          } else if (allele1 == ref && allele2 == alt) {
            "0/1"
          } else if (allele1 == alt && allele2 == ref) {
            "1/0"
          } else if (allele1 == alt && allele2 == alt) {
            "1/1"
          } else {
            "XX"  # If there is any other combination
          }
        }
      })

      # Check that there are no "XX" which would indicate a multiallelic site
      if (!any(genotypes=="XX", na.rm = T)) {

        # Combine VCF fields into a row for the "fix" matrix
        fix_row <- c(chrom, pos, snp_id, ref, alt, qual, filter, info)
        fix <- rbind(fix, fix_row)

        # Add FORMAT and genotypes for each individual to "gt" matrix
        gt_row <- c(format, genotypes)
        gt <- rbind(gt, gt_row)

      } else {
        multiallelic.sites <- multiallelic.sites + 1
      }
    }

    if (multiallelic.sites > 0) {
      print(paste0(multiallelic.sites, " multiallelic sites were removed from the dataset"))
    }

    # rename rows
    rownames(fix) <- 1:nrow(fix)
    rownames(gt) <- 1:nrow(gt)
  }

  # Build vcfR fields if data is in number genotype format (e.g. 0,1,2, where 0 and 2 are homozygous states and 1 is heterozygous)
  if (data.type == "genotypes") {
    # setup fix dataframe
    fix <- matrix(ncol=8, nrow=0)
    colnames(fix) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

    # setup gt dataframe
    gt <- matrix(ncol=nrow(g)+1, nrow = 0)
    colnames(gt) <- c("FORMAT", rownames(g))

    # transpose the genotype matrix
    g <- t(g)

    # Placeholders for each SNP's VCF metadata
    chrom <- "NA"  # Chromosome placeholder
    pos <- "NA"    # Position placeholder (you can modify this if you have position data)
    ref <- "A"     # Assign arbitrary nucleotide as reference
    alt <- "T"     # Assign arbitrary nucleotide as alternate
    qual <- "."     # Quality placeholder
    filter <- "PASS"  # Placeholder for filter
    info <- "."     # Placeholder for info
    format <- "GT"  # Format for genotype

    # initialize counter for multiallelic sites
    multiallelic.sites <- 0

    for (snp_idx in 1:nrow(g)) {
      # SNP id for metadata
      snp_id <- rownames(g)[snp_idx]  # Unique SNP ID

      # Extract genotypes for each sample for this SNP
      genotypes <- sapply(1:ncol(g), function(i) {
        genotype <- as.character(g[snp_idx,i])
        # Check for missing data or "N" alleles
        if (is.na(genotype) || genotype == missing.value) {
          NA  # Missing data format
        } else {
          # Convert genotypes to 0/1, 1/1, etc. based on REF/ALT
          if (genotype == 0) {
            "0/0"
          } else if (genotype == 1) {
            "0/1"
          } else if (genotype == 2) {
            "1/1"
          } else {
            "XX"  # If there is any other combination
          }
        }
      })

      # Check that there are no "XX" which would indicate a multiallelic site
      if (!any(genotypes=="XX", na.rm = T)) {
        # Combine VCF fields into a row for the "fix" matrix
        fix_row <- c(chrom, pos, snp_id, ref, alt, qual, filter, info)
        fix <- rbind(fix, fix_row)

        # Add FORMAT and genotypes for each individual to "gt" matrix
        gt_row <- c(format, genotypes)
        gt <- rbind(gt, gt_row)

      } else {
        multiallelic.sites <- multiallelic.sites + 1
      }
    }

    if (multiallelic.sites > 0) {
      print(paste0(multiallelic.sites, " multiallelic sites were removed from the dataset"))
    }

    # rename rows
    rownames(fix) <- 1:nrow(fix)
    rownames(gt) <- 1:nrow(gt)
  }

  # make empty vcfR object
  v <- new("vcfR")

  # populate "meta" field
  v@meta <- c("##fileformat=VCFv4.2",
              "##source=triangulaR",
              "##INFO=<ID=.,Number=.,Type=.,Description=\"Placeholder for extra information\">",
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")

  # populate "fix" field
  v@fix <- fix

  # population "gt" field
  v@gt <- gt

  return(v)
}



