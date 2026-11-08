# ==============================================================================
# HEATMAP DE BGCs POR ESPÉCIE 
# ==============================================================================
library(tidyverse)
library(readxl)
library(ggplot2)

arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"

# Função para leitura correta dos cabeçalhos
ler_aba_com_cabecalho <- function(arquivo, aba, coluna_ancora) {
  previa <- suppressMessages(read_excel(arquivo, sheet = aba, n_max = 15, col_names = FALSE, col_types = "text"))
  linha_cabecalho <- which(apply(previa, 1, function(linha) any(str_detect(str_squish(linha), fixed(coluna_ancora, ignore_case = TRUE)))))[1]
  if (is.na(linha_cabecalho)) linha_cabecalho <- 1
  df <- suppressMessages(read_excel(arquivo, sheet = aba, skip = linha_cabecalho - 1, col_types = "text"))
  colnames(df) <- str_squish(colnames(df))
  return(df)
}

# 1. Leitura de Metadados e contagem de genomas por Espécie
meta <- ler_aba_com_cabecalho(arquivo_excel, "Metadados 2", "Assecion Number") %>%
  transmute(
    Accession    = str_extract(`Assecion Number`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    GTDB_Species = str_squish(`GTDB species`)
  ) %>%
  filter(!is.na(Accession)) %>%
  mutate(GTDB_Species = if_else(is.na(GTDB_Species) | GTDB_Species %in% c("", "-", "NA"), "Espécie Indefinida", GTDB_Species))

contagem_membros <- meta %>%
  group_by(GTDB_Species) %>%
  summarise(Total_Genomas = n_distinct(Accession), .groups = "drop")

# 2. Leitura do Parvoma
parvoma <- ler_aba_com_cabecalho(arquivo_excel, "TABELA FINAL - PARVOMA", "Acess_number") %>%
  mutate(Accession = str_extract(`Acess_number`, "(GCF|GCA)_[0-9]+\\.[0-9]+")) %>%
  fill(Accession, .direction = "down") %>%
  filter(!is.na(Accession))

df_cruzado <- parvoma %>%
  left_join(meta, by = "Accession") %>%
  filter(!is.na(GTDB_Species))

# 3. Extração dos BGCs e cálculo de porcentagens por espécie
bgcs_exatos <- df_cruzado %>%
  filter(!is.na(`Região`) & !str_squish(`Região`) %in% c("", "-", "NA", "0")) %>%
  mutate(Tipo_BGC = str_squish(`Tipo de Cluster`)) %>%
  filter(!is.na(Tipo_BGC) & Tipo_BGC != "" & Tipo_BGC != "-")

bgc_pct <- bgcs_exatos %>%
  group_by(GTDB_Species, Tipo_BGC) %>%
  summarise(Contagem = n(), .groups = "drop") %>%
  group_by(GTDB_Species) %>%
  mutate(
    Total_BGCs_Especie = sum(Contagem),
    Porcentagem = round((Contagem / Total_BGCs_Especie) * 100, 1)
  ) %>%
  ungroup()

# 4. Expansão da base para garantir as células brancas
todas_combinacoes <- expand_grid(
  GTDB_Species = unique(meta$GTDB_Species),
  Tipo_BGC     = unique(bgc_pct$Tipo_BGC)
)

bgc_heatmap_data <- todas_combinacoes %>%
  left_join(bgc_pct, by = c("GTDB_Species", "Tipo_BGC")) %>%
  replace_na(list(Contagem = 0, Porcentagem = 0)) %>%
  left_join(contagem_membros, by = "GTDB_Species")

# 5. Ordenação dos fatores 
ordem_especies_df <- bgc_heatmap_data %>%
  distinct(GTDB_Species, Total_Genomas) %>%
  arrange(Total_Genomas, GTDB_Species)

ordem_bgcs <- bgc_heatmap_data %>%
  group_by(Tipo_BGC) %>%
  summarise(Soma_Pct = sum(Porcentagem), .groups = "drop") %>%
  arrange(desc(Soma_Pct)) %>%
  pull(Tipo_BGC)

bgc_heatmap_data <- bgc_heatmap_data %>%
  mutate(
    GTDB_Species = factor(GTDB_Species, levels = ordem_especies_df$GTDB_Species),
    Indice_Y     = as.numeric(GTDB_Species),
    Tipo_BGC     = factor(Tipo_BGC, levels = ordem_bgcs),
    Indice_X     = as.numeric(Tipo_BGC)
  )

# 6. Plotagem do Heatmap com 'n = X' 
p_heatmap <- ggplot(bgc_heatmap_data, aes(x = Indice_X, y = Indice_Y, fill = Porcentagem)) +
  geom_tile(color = "gray80", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white", high = "#08306b", name = "Frequência (%)", 
    limits = c(0, 100), breaks = seq(0, 100, 25)
  ) +
  scale_x_continuous(
    breaks = seq_along(ordem_bgcs),
    labels = ordem_bgcs,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(ordem_especies_df$GTDB_Species),
    labels = ordem_especies_df$GTDB_Species,
    expand = c(0, 0),
    sec.axis = dup_axis(
      name = NULL,
      breaks = seq_along(ordem_especies_df$GTDB_Species),
      labels = paste0("n = ", ordem_especies_df$Total_Genomas)
    )
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Perfil Biossintético por Espécie",
    x = "Classe do Cluster (antiSMASH)",
    y = expression(paste("Espécies de ", italic("Maribacter"), " (GTDB)"))
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 9.5),
    axis.text.y.left = element_text(angle = 0, hjust = 1, face = "italic", color = "black", size = 9.5),
    axis.text.y.right = element_text(angle = 0, hjust = 0, face = "bold", color = "black", size = 9.5),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13, margin = margin(b = 12)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 15, r = 20, b = 30, l = 15)
  )

print(p_heatmap)

# 7. Salvar em alta resolução 
ggsave("Figura17.png", plot = p_heatmap, width = 10, height = 10.5, dpi = 300, bg = "white")

cat("\nFigura 'Figura17.png' gerada com sucesso!\n")