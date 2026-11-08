library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# 1. Importação 
dados <- read_excel("Maribacter - Parvome + Cazome.xlsx", 
                    sheet = "TABELA FINAL - PARVOMA", 
                    skip = 2)

# Renomeando na tabela provisória do R pela posição das colunas:
# 1 = Acess_number | 5 = ID_GenBank_RefSeq | 7 = Região | 8 = Tipo de Cluster
names(dados)[1] <- "Acess_number"
names(dados)[5] <- "Contig_ID"
names(dados)[7] <- "Regiao"
names(dados)[8] <- "Tipo_Cluster"

# 2. Tabela provisória corrigindo células mescladas e mantendo híbridos individuais
dados_heatmap <- dados %>%
  
  # Remove linhas que estejam completamente em branco na planilha
  filter(!if_all(everything(), is.na)) %>%
  
  # Limpa quebras de linha invisíveis que o Excel coloca em células mescladas (ex: "GCF_9475\n00985.1")
  mutate(
    Acess_number = str_remove_all(Acess_number, "\\s+"),
    Contig_ID    = str_trim(Contig_ID),
    Regiao       = as.character(str_trim(Regiao)),
    Tipo_Cluster = str_trim(Tipo_Cluster)
  ) %>%
  
  # Preenche de cima para baixo as células que vieram como NA por causa da mesclagem
  fill(Acess_number, Contig_ID, Regiao, .direction = "down") %>%
  
  # Filtra apenas linhas com genoma e cluster válidos
  filter(!is.na(Acess_number), !is.na(Tipo_Cluster)) %>%
  filter(!Tipo_Cluster %in% c("-", "N/A", "", "Tipo de Cluster")) %>%
  
  distinct(Acess_number, Contig_ID, Regiao, Tipo_Cluster) %>%
  
  # Conta a frequência correta de cada classe por genoma
  count(Acess_number, Tipo_Cluster)

# 3. Ordenação 
total_por_genoma <- dados_heatmap %>%
  group_by(Acess_number) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(total)

dados_heatmap$Acess_number <- factor(dados_heatmap$Acess_number, 
                                     levels = total_por_genoma$Acess_number)

# 4. Plot do Heatmap
plot_heatmap <- ggplot(dados_heatmap, aes(x = Tipo_Cluster, y = Acess_number, fill = n)) +
  
  geom_tile(color = "white", lwd = 0.1) + 
  
  scale_fill_gradient(
    low = "#41b6c4", 
    high = "#081d58", 
    name = "Count",
    breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)
  ) +
  
  # Títulos
  labs(
    title = "Global Distribution of BGCs in Maribacter",
    subtitle = "Heatmap: 128 analyzed genomes",
    x = "Cluster Family",
    y = "Genome"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(plot_heatmap)

# 5. Salvar em alta resolução
ggsave("Heatmap_Plot_NOVO.png", plot = plot_heatmap, width = 10, height = 14, dpi = 300)
