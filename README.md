# fieldcorn
Functions useful to the fieldcorn lab

## Installation

```
devtools::install_github("acperkins3/fieldcorn")
```

## `make.IBS.table`

A function to calculate identify by state (IBS) between two or more lines from genotypes in [Hapmap format](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load). IBS is calculated in bins of a certain number of markers across the genome.

The current version has no issues I'm aware of but is slow. I will try vectorizing it at some point for speed.

### Example

A hapmap file can be read into `R` using `read.delim()`. We will use the example data contained in the package here

```
library(tidyverse)
library(fieldcorn)

table1 <- make.IBS.table(ExampleHapmap1, ReferenceLine = "PHN46", OtherLines = c("PHP02", "PHJ89"))

```

Here is what that table looks like

<p align="center"><img src="https://raw.githubusercontent.com/acperkins3/GDD-Plots-2022/main/ReadmeImages/IBSTable.png" /></p>

The `Genotype` column is the lines being compared with the `ReferenceLine`, PHN46.

We can plot that easily

```
table1 %>%
  mutate(BinStartPosMb = BinStartPos / 1000000) %>%
  ggplot(aes(BinStartPosMb, ProportionIBS, color = Genotype)) +
  geom_line() +
  facet_wrap(vars(Chromosome), scales = "free_x", nrow=2) +
  theme_bw() +
  theme(legend.title=element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(), legend.position="bottom", axis.text.x = element_text(size=8), axis.text.y = element_text(size=8)) +
  xlab("Position (Mb)") +
  ylab("Proportion IBS to PHN46")
```

<p align="center"><img src="https://raw.githubusercontent.com/acperkins3/GDD-Plots-2022/main/ReadmeImages/IBSTableExample.png" /></p>