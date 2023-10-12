
<!-- README.md is generated from README.Rmd. Please edit that file -->

# triangulaR

<!-- badges: start -->

<!-- badges: end -->

triangulaR is a package for calculating hybrid indices, heterozygosity,
and building triangle plots

## Installation

You can install the development version of triangulaR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("omys-omics/triangulaR")
```

## Citation

If you publish any work that uses triangulaR, please consider citing the
following paper:

## How to use this package:

triangulaR builds on the functionality of
[vcfR](https://knausb.github.io/vcfR_documentation/) to analyze SNP data
in R. The first step is to read in the data from a vcf file.

#### Step 1: Read in data

It is expected that the data have already filtered for quality
(e.g. setting genotype quality and depth thresholds, removing
individuals and sites with high missing data, etc.). For help with this
step, see [SNPfiltR](https://github.com/DevonDeRaad/SNPfiltR).

``` r
library(triangulaR)
#> Loading required package: ggplot2
#> Warning: package 'ggplot2' was built under R version 4.0.5
library(vcfR)
#> 
#>    *****       ***   vcfR   ***       *****
#>    This is vcfR 1.14.0 
#>      browseVignettes('vcfR') # Documentation
#>      citation('vcfR') # Citation
#>    *****       *****      *****       *****

# Read in data
data <- read.vcfR("../SecondaryContact/gen.19000.vcf", verbose = F)
data
#> ***** Object of Class vcfR *****
#> 420 samples
#> 1 CHROMs
#> 1,788 variants
#> Object size: 6.2 Mb
#> 0 percent missing data
#> *****        *****         *****

# Or, use example vcfR object
example.vcfR
#> ***** Object of Class vcfR *****
#> 420 samples
#> 1 CHROMs
#> 1,215 variants
#> Object size: 4.2 Mb
#> 0 percent missing data
#> *****        *****         *****
```

#### Step 2: Make a popmap

A popmap is a data.frame with two columns labeled “id” and “pop”. Each
name of each sample in the vcfR object should be included in the “id”
column. All individuals in the vcfR object should be included in the
popmap, and vice versa. Each individual needs to be assigned to a
population. Every individual must be assigned to a population, and there
can be any number of populations. IDs and pops should be character
strings.

``` r
# Here is an example of what a popmap should look like
print(head(example.popmap))
#>   id pop
#> 1 i0   0
#> 2 i1   0
#> 3 i2   0
#> 4 i3   0
#> 5 i4   0
#> 6 i5   0
print(tail(example.popmap))
#>       id pop
#> 415 i414  20
#> 416 i415  20
#> 417 i416  20
#> 418 i417  20
#> 419 i418  20
#> 420 i419  20
```

#### Step 3: Choose sites above an allele frequency difference threshold

Theoretically, the hybrid index of an individual represents the
proportion of ancestry received from each parental population. In
practice, one way to calculate hybrid indices is by identifying loci
with allele frequency differences above a chosen threshold in the
parental populations and scoring individuals by allele counts at those
loci. There is a balance between using a small amount of highly
diagnostic site (e.g. fixed differences) and a large amount of less
diagnostic sites. I recommend trying difference values for the allele
frequency difference threshold to see how this value affects results.

``` r
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
example.vcfR.diff <- alleleFreqDiff(vcfR = example.vcfR, pm = example.popmap, p1 = 0, p2 = 20, difference = 0.6)
#> [1] "216 sites passed allele frequency difference threshold"
```

#### Step 4: Calculate hybrid index and heterozygosity for each sample

Once sites above the allele frequency difference threshold have been
identified, hybrid index and heterozygosity for each sample can be
calculated.

``` r
# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
hi.het <- hybridIndex(vcfR = example.vcfR.diff, pm = example.popmap, p1 = 0, p2 = 20)
#> [1] "calculating hybrid indices and heterozygosities based on 216 sites"
```

#### Step 5: Visualize results as a triangle plot

``` r
# Generate colors (or leave blank to use default)
cols <- colorRampPalette(colors = c("#313695", "khaki2", "#a50026"))
# View triangle plot
triangle.plot(hi.het, colors = cols(21))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

#### Step 6: Color triangle plot by missing data

The data.frame returned by the hybridIndex function also contains the
percent of missing data in each individual. View the triangle plot with
samples colored by percent missing data to investigate its effect.

``` r
# There is no missing data in this data.set, so all samples have the same color
missing.plot(hi.het)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />
