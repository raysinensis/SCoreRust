# SCore-rust
Single Cell Pathway Scoring powered by Rust

Speed can be workflow altering. By reimplementing pathway (KEGG/Reactome) scoring algorithms for single cell sequencing, higher analysis throughput, richer insights, and better interactivity can be achieved.

Wrappers to call the faster implementations directly in R are designed as in-place replacements for SingleCellExperiment and Seurat workflows.

<img src="inst/bench1.png" width="600" align="center">

```
# testing on a small example
library(tidyverse)
library(Seurat)
so <- pbmc_small

# yielding similar scores, note that the random sampling isn't exactly the same
a <- ggplot(data.frame(Seurat = resa, SCore = resb), aes(x = Seurat, y = SCore)) + 
  geom_point(size = 1) +
  theme_classic() +
  coord_fixed()
  
# benchmark
res <- bench::mark(
  Seurat = AddModuleScore(so, features = list(c("LST1", "AIF1", "PSAP", "YWHAB", "MYO1G", "SAT1")), nbin = 10, ctrl = 10),
  SCore = calc_modulescore(as.matrix(so@assays$RNA@data),
                           c("LST1", "AIF1", "PSAP", "YWHAB", "MYO1G", "SAT1"),
                           rownames(as.matrix(so@assays$RNA@data))),
  iterations = 50,
  check = FALSE
)
b <- ggplot(res %>%
         mutate(expression = as.character(expression)) %>%
         select(expression, time) %>%
         unnest_longer(time), 
       aes(x = expression, y = time)) +
  geom_boxplot(aes(fill = expression), outlier.shape = NA) +
  geom_jitter(size = 0.5, alpha = 0.5) +
  theme_classic() +
  ylab("time (ms)") +
  NoLegend()
```
