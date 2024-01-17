library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(GGally)




metadata <-
  read_delim(
    "./meta.tsv",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

colnames(metadata)[1] <- "SampleID"
metadata <- metadata[metadata$body.site %in% c("gut", "tongue"), ]

kegg_abundance <-
  ko2kegg_abundance(
    "./pred_metagenome_unstrat.tsv"
  )


kegg_abundance2 <- kegg_abundance[, metadata$SampleID]

ko_AlDEx <-
  pathway_daa(
    abundance = kegg_abundance2,
    metadata = metadata,
    group = "body.site", 
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )
ko_AlDEx_df <-   ko_AlDEx[ko_AlDEx$method == "ALDEx2_Wilcoxon rank test", ]

ko_annotation <-pathway_annotation(pathway = "KO",
                                   daa_results_df = ko_AlDEx_df, 
                                   ko_to_kegg = TRUE)


Top30 <- ko_annotation %>% arrange(p_adjust) %>% top_n(30)
kegg_abundance_t30 <- kegg_abundance2[rownames(kegg_abundance2) %in% Top30$feature, ]


p <- pathway_errorbar(abundance = kegg_abundance_t30,
                      daa_results_df = Top30,
                      Group = metadata$body.site,
                      ko_to_kegg = T,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = T,
                      colors = NULL,
                      x_lab = NULL)

p




