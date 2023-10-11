#' alleleFreqDiff
#'
#' @param vcfR vcfR object
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each individual must be assigned to a population
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
#' #alleleFreqDiff(vcfR = vcfR, pm = pm, p1 = p1, p2 = p2, difference = 0.5)
alleleFreqDiff <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)))[-1], sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  # make matrices of genotypes for each parental pop
  p1.gts <- data.frame(extract.gt(vcfR[samples = c(pm[pm$pop == p1,"id"])[[1]]]))
  p2.gts <- data.frame(extract.gt(vcfR[samples = c(pm[pm$pop == p2,"id"])[[1]]]))

  # make dataframe to keep track of frequency of "1" allele at each locus
  af <- data.frame(matrix(nrow = nrow(p1.gts), ncol = 2))
  colnames(af) <- c("p1", "p2")
  rownames(af) <- rownames(p1.gts)

  # recode missing data
  p1.gts[is.na(p1.gts)] <- -9
  p2.gts[is.na(p2.gts)] <- -9
  # recode alelles
  p1.gts[p1.gts=="0|0"] <- 0
  p1.gts[p1.gts=="0|1"] <- 1
  p1.gts[p1.gts=="1|0"] <- 1
  p1.gts[p1.gts=="1|1"] <- 2
  p1.gts[p1.gts=="0/0"] <- 0
  p1.gts[p1.gts=="0/1"] <- 1
  p1.gts[p1.gts=="1/0"] <- 1
  p1.gts[p1.gts=="1/1"] <- 2
  p2.gts[p2.gts=="0|0"] <- 0
  p2.gts[p2.gts=="0|1"] <- 1
  p2.gts[p2.gts=="1|0"] <- 1
  p2.gts[p2.gts=="1|1"] <- 2
  p2.gts[p2.gts=="0/0"] <- 0
  p2.gts[p2.gts=="0/1"] <- 1
  p2.gts[p2.gts=="1/0"] <- 1
  p2.gts[p2.gts=="1/1"] <- 2

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  for (i in 1:nrow(af)) {
    # p1:
    # get genotypes for locus, remove missing data
    loc.p1 <- p1.gts[i,]
    loc.p1 <- loc.p1[loc.p1!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p1"] <- sum(loc.p1)/(length(loc.p1)*2)

    #p2:
    # get genotypes for locus, remove missing data
    loc.p2 <- p2.gts[i,]
    loc.p2 <- loc.p2[loc.p2!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p2"] <- sum(loc.p2)/(length(loc.p2)*2)
  }
  # add allele frequency differences to af
  af$diff <- abs(af[,"p1"] - af[,"p2"])

  # get names of loci with allele frequency difference above threshold
  loci <- rownames(af[af$diff >= difference,])

  # get names of all loci in vcfR object
  all.loci <- rownames(extract.gt(vcfR))

  # get indices in vcfR object of loci with difference above theshold
  loci.indices <- all.loci %in% loci
  loci.indices <- which(loci.indices)

  # subset vcfR by loci with allele frequencies difference in parental pops above threshold
  vcfR.diff <- vcfR[loci.indices]

  # print statement
  s <- nrow(vcfR.diff@gt)
  print(paste0(s, " sites passed allele frequency difference threshold"))
  return(vcfR.diff)
}



#' hybridIndex
#'
#' @param vcfR data in vcfR format
#' @param pm data.frame containing two columns, "id" and "pop". The ids must match those in the vcfR object, but don't need to be in the same order. Each individual must be assigned to a population
#' @param p1 (character) name of parental population 1
#' @param p2 (character) name of parental population 2
#'
#' @return dataframe containing hybrid indices, heterozygosities, and percent missing data
#' @export
#'
#' @importFrom vcfR extract.gt is.het
#'
#' @examples
#' #hybridIndex(vcfR = vcfR, pm = pm, p1 = p1, p2 = p2)
hybridIndex <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }

  if (!(identical(sort(unique(colnames(vcfR@gt)))[-1], sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }

  # get number of differences above threshold
  d <- nrow(vcfR@gt)
  print(paste0("calculating hybrid indices and heterozygosities based on ", d, " sites"))

  m <- extract.gt(vcfR)
  # recode missing data
  m[is.na(m)] <- -9
  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2
  # make new matrix of same size as m
  n <- matrix(nrow = nrow(m), ncol = ncol(m))

  # make matrices of genotypes for each parental pop
  p1.gts <- data.frame(m[,c(pm[pm$pop == p1,"id"])[[1]]])
  p2.gts <- data.frame(m[,c(pm[pm$pop == p2,"id"])[[1]]])

  # make dataframe to keep track of frequency of "1" allele at each locus
  af <- data.frame(matrix(nrow = nrow(m), ncol = 2))
  colnames(af) <- c("p1", "p2")
  rownames(af) <- rownames(m)

  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)

  for (i in 1:nrow(m)) {

    # p1:
    # get genotypes for locus, remove missing data
    loc.p1 <- p1.gts[i,]
    loc.p1 <- loc.p1[loc.p1!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p1"] <- sum(loc.p1)/(length(loc.p1)*2)

    #p2:
    # get genotypes for locus, remove missing data
    loc.p2 <- p2.gts[i,]
    loc.p2 <- loc.p2[loc.p2!=-9]
    # add allele frequency of 1 to af dataframe
    af[i,"p2"] <- sum(loc.p2)/(length(loc.p2)*2)

    # for each locus, if parental pop 1 has a higher frequency of "2" allele, assign p1.allele as 2
    if (af[i,"p1"] > af[i,"p2"]) {
      p1.allele <- 2
      # else, if parental pop1 has a lower frequency of "2 allele", assign p1.allele as 0
    } else {
      p1.allele <- 0
    }

    # compare every individual to parental pop 1, giving scores of:
    # 0 = matching p1
    # 1 = being a het
    # 2 = not matching p1 (matching p2)
    # -9 = NA
    for (j in 1:ncol(m)) {
      if (m[i,j] == -9) {
        next
      }
      if (m[i,j] == p1.allele) {
        n[i,j] <- 0
      } else if (m[i,j] == 1) {
        n[i,j] <- 1
      } else {
        n[i,j] <- 2
      }
    }
  }
  colnames(n) <- colnames(m)
  rownames(n) <- rownames(m)
  # count alleles, removing NAs
  counts <- colSums(n, na.rm = T)
  # count nonmissing genotypes for each ind
  sites <- colSums(!apply(n, MARGIN = 2, is.na))
  # calculate hybrid index
  hi <- counts / (sites*2)
  # make dataframe
  tri <- data.frame(matrix(nrow = nrow(pm), ncol = 5))
  colnames(tri) <- c("id", "pop", "hybrid.index", "heterozygosity", "perc.missing")
  # add id column
  tri$id <- names(hi)
  tri$hybrid.index <- hi
  # add pop to hybrid index dataframe
  for (ind in pm$id) {
    tri[tri$id == ind, "pop"] <- pm[pm$id == ind, "pop"]
  }
  # get heterozygosity (don't count NAs as FALSE, leave as NA)
  het <- is.het(extract.gt(vcfR), na_is_false = F)
  for (ind in pm$id) {
    # divide the number of sites that are het by the number of sites that have data
    tri[tri$id == ind, "heterozygosity"] <- sum(het[,ind], na.rm = T)/sum(!is.na(het[,ind]))
  }
  #calculate missingness by individual
  ms <- colSums(is.na(vcfR@gt))/nrow(vcfR@gt)
  #add missingness to df
  tri$perc.missing <- ms[-1]

  return(tri)
}


#' triangle.plot
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



