library(phyloseq)   # 마이크로바이옴 데이터 분석 및 시각화
library(tidyverse)  # R의 데이터 핸들링 
library(readr)      # 파일 읽어오기 
library(ggpicrust2) # PICRUSt2결과 처리 및 시각화
library(ALDEx2)     # Aldex2 분석
library(pheatmap)   # heatmap 시각화


library(tidyverse)
library(psych)
library(pheatmap)
library(magrittr)
library(scico)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(GGally)


Group <- factor(metadata$group)
Level <- levels(Group)

## ALDEx2 분석 #### 
# 1) Aldex2분석을 위해 OTU read count 데이터를 CLR로 normalization 
ALDEx2_object <- ALDEx2::aldex.clr(round(path_counts), metadata$group)
# clr 변환은 각 샘플의 모든 taxa값이 0이며, taxa간의 상대적 비율 값이 그대로 보존된다

# 2) Aldex 통계 분석
ALDEx2_results <- ALDEx2::aldex.ttest(ALDEx2_object, paired.test = FALSE, verbose = FALSE)
ALDEx2_results <- bind_cols(notation_path, ALDEx2_results)
# t.test결과와 wilcoxon rank sum test결과를 반환한다 
# 각 taxa가 그룹내에서 유의한 수준으로 차이가 '있는지 없는지' 판별

# 3) Effect size계산
ALDEx2_effect <- ALDEx2::aldex.effect(ALDEx2_object)
ALDEx2_effect <- bind_cols(notation_path, ALDEx2_effect)
# 각 taxa가 차이가 있다면(유의하다면), '얼마나' 차이가 있는지 판별


# 4) 결과물 정리 
p_values_df <- data.frame(
  feature = rep(rownames(ALDEx2_results), 2), 
  description = rep(ALDEx2_results$description),
  method = c(rep("ALDEx2_Welch's t test", nrow(ALDEx2_results)), 
             rep("ALDEx2_Wilcoxon rank test", nrow(ALDEx2_results))),
  group1 = rep(Level[1], 2 * nrow(ALDEx2_results)),
  group2 = rep(Level[2], 2 * nrow(ALDEx2_results)),
  p_values = c(ALDEx2_results$we.ep, ALDEx2_results$wi.ep),
  effect = ALDEx2_effect$effect)

# 5) P- value값 보정 (BH 사용)
adjusted_p_values <- data.frame(
  feature = p_values_df$feature,
  p_adjust = p.adjust(p_values_df$p_values, method = "BH"))

# 6) 결과물 최종 정리
path_AlDEx_result <- cbind(p_values_df, p_adjust = adjusted_p_values$p_adjust)
path_AlDEx_result

path_ALDEx2 <- cbind(notation_path, ALDEx2_results, ALDEx2_effect, path_counts)




## Heatmap 그리기 ####

# "ALDEx2_Wilcoxon rank test"결과 중에서 p-value값 보정이 유의한 것만 추출 
path_AlDEx_result2 <- path_AlDEx_result[path_AlDEx_result$method == "ALDEx2_Wilcoxon rank test" &
                                        path_AlDEx_result$p_values < 0.05, ]
path_AlDEx_result2 <- path_AlDEx_result2[,-2] 
dim(path_AlDEx_result2)

# ggpicrust2의 annotation함수를 통해 이 kegg pathway데이터 가져오기
path_ann <-pathway_annotation(pathway = "MetaCyc",
                            daa_results_df = path_AlDEx_result2)

dim(path_ann)




### a) 전단계?  #######
# Pathway 1,2,3을 하나의 이름으로 표시하기
rownames(path_ann) <- path_ann$feature

# annotation이 완료된 kegg만 추출하기
pathway_abundance3 <- path_counts[path_ann$feature, ]


# kegg abundance를 relativer abundance로 바꾸기
col_sums <- colSums(pathway_abundance3)
path_relative_abundance <- t(t(pathway_abundance3) / col_sums)

# heatmap을 위한 데이터 생성
table <- merge(path_ann, path_relative_abundance, by = "row.names") 
rownames(table) <- table$feature


### b) heatmap 그리기 #####
library(pheatmap)
p <- pheatmap(table[, 10:21])
p


png("heatmap/Heatmap.png", width = 15, height = 20, units = "in", res = 600)
print(p)
dev.off()

rownames(table) <- table$description




### c) 보정하기 #######
annotation_row <- data.frame(
  row.names  = rownames(table),
  p_adjust = table$p_adjust,
  Effect = table$effect
)

annotation_row <- annotation_row[order(annotation_row$Effect, decreasing = F), ,
                                 drop = FALSE # rownames 사라지는걸 막음
]


annotation_col <- data.frame(
  row.names = metadata$sample,
  Group =  metadata$group
)

annotation_col$Group <- factor(annotation_col$Group , levels = c("CON", "AD") ) 

ann_colors = list(
  p_adjust = colorRampPalette(c("white", "green"))(100),
  Effect = colorRampPalette(c("red", "white", "blue"))(100),
  Group = c("CON" = "orange","AD"="purple")
)



df <- table[, 10:21]
colnames(df)
df <- df[rownames(annotation_row), rownames(annotation_col)]


pvalue <- table$p_values
pvalue <- round(pvalue, digits = 3)


plot <- pheatmap::pheatmap(mat = as.matrix(df),
                           color = colorRampPalette(c("blue", "white", "red"))(100),
                           annotation_col = annotation_col, 
                           annotation_row = annotation_row, 
                           annotation_colors = ann_colors,
                           scale = "row", # 행별로 정규화
                           cluster_rows = F,  
                           cluster_cols = F,
                           display_numbers = pvalue,
                           gaps_col = 6, # 8번째 열에서 heatmap 분리 
                           legend = T,
                           border_color=NA)


png("heatmap/Adjust Heatmap.png", width = 20, height = 20, units = "in", res = 600)
print(plot)
dev.off()







## Error bar 그리기 ######

### With KO ########

kegg_abundance <-
  ko2kegg_abundance(
    "./pred_metagenome_unstrat.tsv"
  )

ko_AlDEx <-
  pathway_daa(
    abundance = kegg_abundance,
    metadata = metadata,
    group = "group", 
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )

ko_AlDEx_df <-   ko_AlDEx[ko_AlDEx$method == "ALDEx2_Wilcoxon rank test", ]
write.csv(ko_AlDEx_df, "ref/KO.tsv")

ko_AlDEx_df2 <- read.csv("ref/KO.tsv", sep = "\t", header = T)
ko_annotation <-pathway_annotation(pathway = "KO",
                                   daa_results_df = ko_AlDEx_df2, 
                                   ko_to_kegg = TRUE)



Top30 <- ko_annotation %>% arrange(p_adjust) %>% top_n(30)
kegg_abundance_t30 <- kegg_abundance[rownames(kegg_abundance) %in% Top30$feature, ]


p <- pathway_errorbar(abundance = kegg_abundance_t30,
                      daa_results_df = Top30,
                      Group = metadata$group,
                      ko_to_kegg = T,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = T,
                      colors = NULL,
                      x_lab = NULL)






### With MetaCyc ########


path_AlDEx_df <-   path_AlDEx_result[path_AlDEx_result$method == "ALDEx2_Wilcoxon rank test" &
                                     path_AlDEx_result$p_values < 0.05, ]


Top <- path_AlDEx_df %>% arrange(p_values) %>% top_n(30)
write.csv(Top, "ref/Top.tsv")

Top <- read.csv("ref/Top.tsv", sep = "\t", header = T)

path_rela_t30 <- path_rela[rownames(path_rela) %in% Top$feature, ]


p <- pathway_errorbar(abundance = path_rela_t30,
                      daa_results_df = Top,
                      Group = metadata$group,
                      ko_to_kegg = F,
                      p_values_threshold = 0.05,
                      order = "description",
                      select = NULL,
                      p_value_bar = T,
                      colors = NULL,
                      x_lab = NULL)


