library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)

arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"

ler_aba_com_cabecalho <- function(arquivo, aba, coluna_ancora) {
  previa <- suppressMessages(read_excel(arquivo, sheet = aba, n_max = 15, col_names = FALSE, col_types = "text"))
  linha_cabecalho <- which(apply(previa, 1, function(linha) any(str_detect(str_squish(linha), fixed(coluna_ancora, ignore_case = TRUE)))))[1]
  if (is.na(linha_cabecalho)) linha_cabecalho <- 1
  df <- suppressMessages(read_excel(arquivo, sheet = aba, skip = linha_cabecalho - 1, col_types = "text"))
  colnames(df) <- str_squish(colnames(df))
  return(df)
}

# 1. Ler Metadados 2
meta <- ler_aba_com_cabecalho(arquivo_excel, "Metadados 2", "Assecion Number") %>%
  transmute(
    Accession    = str_extract(`Assecion Number`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    GTDB_Species = str_squish(`GTDB species`),
    Total_CDS    = as.numeric(str_replace_all(`CDS Count`, "[^0-9.]", ""))
  ) %>%
  filter(!is.na(Accession)) %>%
  mutate(
    GTDB_Species = if_else(is.na(GTDB_Species) | GTDB_Species %in% c("", "-", "NA"), "Espécie Indefinida", GTDB_Species)
  )

# 2. Ler aba CAZYME - dbCAN4
cazy_raw <- ler_aba_com_cabecalho(arquivo_excel, "CAZYME - dbCAN4", "GCF_|GCA_|Acess|Access|Genom|Gene")
col_acc_cazy <- colnames(cazy_raw)[apply(head(cazy_raw, 30), 2, function(col) any(str_detect(col, "(GCF|GCA)_[0-9]+\\.[0-9]+"), na.rm = TRUE))][1]

cazy_processado <- cazy_raw %>%
  mutate(Accession = str_extract(.data[[col_acc_cazy]], "(GCF|GCA)_[0-9]+\\.[0-9]+")) %>%
  fill(Accession, .direction = "down") %>%
  filter(!is.na(Accession))

col_total_cazy <- colnames(cazy_processado)[str_detect(colnames(cazy_processado), regex("^Total|Total_CAZ", ignore_case = TRUE))][1]
col_gene       <- colnames(cazy_processado)[str_detect(colnames(cazy_processado), regex("Gene|ID|ORF|Protein|Query|HMMER", ignore_case = TRUE))][1]

if (!is.na(col_total_cazy)) {
  contagem_cazymes <- cazy_processado %>%
    mutate(Val_Total = as.numeric(str_replace_all(.data[[col_total_cazy]], "[^0-9.]", ""))) %>%
    filter(!is.na(Val_Total)) %>%
    group_by(Accession) %>%
    summarise(Total_CAZymes = first(Val_Total), .groups = "drop")
} else {
  contagem_cazymes <- cazy_processado %>%
    filter(!is.na(.data[[col_gene]]) & !str_squish(.data[[col_gene]]) %in% c("", "-", "NA", col_gene)) %>%
    group_by(Accession) %>%
    summarise(Total_CAZymes = n(), .groups = "drop")
}

# 3. Cruzar dados e calcular o Conteúdo do CAZoma
df_cazoma <- meta %>%
  left_join(contagem_cazymes, by = "Accession") %>%
  replace_na(list(Total_CAZymes = 0)) %>%
  mutate(
    Razao_CAZoma       = round(Total_CAZymes / Total_CDS, 5),
    Porcentagem_CAZoma = round((Total_CAZymes / Total_CDS) * 100, 2)
  ) %>%
  arrange(desc(Porcentagem_CAZoma))

# ==============================================================================
# 4. GRÁFICO 1: VARIAÇÃO DO CONTEÚDO DO CAZOMA POR ESPÉCIE 
# ==============================================================================
ordem_especies_cazoma <- df_cazoma %>%
  group_by(GTDB_Species) %>%
  summarise(
    Mediana_Especie = median(Porcentagem_CAZoma, na.rm = TRUE),
    N_Especie = n(),
    .groups = "drop"
  ) %>%
  arrange(Mediana_Especie, N_Especie)

df_plot_especie <- df_cazoma %>%
  mutate(
    GTDB_Species = factor(GTDB_Species, levels = ordem_especies_cazoma$GTDB_Species),
    Indice_Y = as.numeric(GTDB_Species)
  )

p_cazoma_especie <- ggplot(df_plot_especie, aes(x = Porcentagem_CAZoma, y = Indice_Y, group = Indice_Y)) +
  geom_boxplot(fill = "#c6dbef", color = "#08519c", outlier.shape = NA, alpha = 0.7, width = 0.55) +
  geom_jitter(height = 0.12, width = 0, size = 1.8, color = "#08306b", alpha = 0.85) +
  scale_y_continuous(
    breaks = seq_along(ordem_especies_cazoma$GTDB_Species),
    labels = ordem_especies_cazoma$GTDB_Species,
    expand = expansion(mult = c(0.02, 0.02)),
    sec.axis = dup_axis(
      name = NULL,
      breaks = seq_along(ordem_especies_cazoma$GTDB_Species),
      labels = paste0("n = ", ordem_especies_cazoma$N_Especie)
    )
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Variação do Conteúdo do CAZoma por Espécie",
    subtitle = "Proporção de genes codificadores de CAZymes em relação ao total de CDS (%)",
    x = "Conteúdo do CAZoma (% do total de CDS)",
    y = expression(paste("Espécies de ", italic("Maribacter"), " (GTDB)"))
  ) +
  theme(
    # Colar título e subtítulo totalmente à esquerda
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0, color = "gray30", margin = margin(b = 12)),
    axis.text.y.left = element_text(face = "italic", color = "black", size = 9.5),
    axis.text.y.right = element_text(face = "bold", hjust = 0, color = "black", size = 9.5),
    axis.text.x = element_text(color = "black", size = 10.5),
    axis.title = element_text(color = "black", size = 11.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 15, r = 20, b = 15, l = 15)
  )

print(p_cazoma_especie)
# Altura aumentada (height = 10.5) para dar amplo espaço vertical entre os 40 nomes
ggsave("Variacao_CAZoma_Por_Especie.png", plot = p_cazoma_especie, width = 11.69, height = 10.5, dpi = 300, bg = "white")

# ==============================================================================
# 5. GRÁFICO 2: DISTRIBUIÇÃO GLOBAL DO CAZOMA 
# ==============================================================================
mediana_geral <- median(df_cazoma$Porcentagem_CAZoma, na.rm = TRUE)
cor_mediana   <- "#d95f02" 

p_cazoma_dist <- ggplot(df_cazoma, aes(x = Porcentagem_CAZoma)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "#6baed6", color = "white", alpha = 0.85) +
  geom_density(color = "#08306b", linewidth = 1) +
  geom_vline(xintercept = mediana_geral, linetype = "dashed", color = cor_mediana, linewidth = 1) +
  # Posicionado no topo 
  annotate("text", x = mediana_geral + 0.12, y = 1.35, 
           label = paste0("Mediana = ", mediana_geral, "%"), 
           color = cor_mediana, fontface = "bold", size = 4.2, hjust = 0) +
  theme_minimal(base_size = 12) +
  labs(
    title = bquote(bold("Distribuição Global do Conteúdo do CAZoma em " * italic("Maribacter") * " (n = 128)")),
    x = "Conteúdo do CAZoma (% de CAZymes / Total de CDS)",
    y = "Densidade"
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, size = 13, margin = margin(b = 12)),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 15, r = 20, b = 15, l = 15)
  )

print(p_cazoma_dist)
ggsave("Distribuicao_Global_CAZoma.png", plot = p_cazoma_dist, width = 9.5, height = 6, dpi = 300, bg = "white")

cat("\nGráficos atualizados salvos na sua pasta de trabalho!\n")