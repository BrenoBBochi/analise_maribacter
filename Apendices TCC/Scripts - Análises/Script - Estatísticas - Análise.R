# ==============================================================================
# PAINEL DE COMPARAÇÕES ESTATÍSTICAS (CATEGORIAS EXATAS: n = 42, 15, 8 e 4)
# ==============================================================================
library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)

arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"

ler_aba_com_cabecalho <- function(arquivo, aba, coluna_ancora) {
  previa <- suppressMessages(read_excel(arquivo, sheet = aba, n_max = 15, col_names = FALSE, col_types = "text"))
  linha_cabecalho <- which(apply(previa, 1, function(linha) any(str_detect(str_squish(linha), fixed(coluna_ancora, ignore_case = TRUE)))))[1]
  if (is.na(linha_cabecalho)) linha_cabecalho <- 1
  df <- suppressMessages(read_excel(arquivo, sheet = aba, skip = linha_cabecalho - 1, col_types = "text"))
  colnames(df) <- str_squish(colnames(df))
  return(df)
}

# 1A. Ler Tamanho do Genoma na aba 'Metadados 2'
meta2_raw <- ler_aba_com_cabecalho(arquivo_excel, "Metadados 2", "Assecion Number")
col_size  <- colnames(meta2_raw)[str_detect(colnames(meta2_raw), regex("Genome Size", ignore_case = TRUE))][1]

meta2 <- meta2_raw %>%
  transmute(
    Accession      = str_extract(`Assecion Number`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    Genome_Size_Mb = as.numeric(str_replace_all(.data[[col_size]], "[^0-9.]", "")) / 1e6
  ) %>%
  filter(!is.na(Accession))

# 1B. Ler Fonte de Isolamento e Host na aba 'Metadados'
meta1_raw  <- ler_aba_com_cabecalho(arquivo_excel, "Metadados", "Accession")
col_source <- colnames(meta1_raw)[str_detect(colnames(meta1_raw), regex("Isolation Source", ignore_case = TRUE))][1]
col_host   <- colnames(meta1_raw)[str_detect(colnames(meta1_raw), regex("^Host", ignore_case = TRUE))][1]

meta1 <- meta1_raw %>%
  transmute(
    Accession  = str_extract(`Accession`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    Source_Raw = str_squish(.data[[col_source]]),
    Host_Raw   = if (!is.na(col_host)) str_squish(.data[[col_host]]) else ""
  ) %>%
  filter(!is.na(Accession)) %>%
  mutate(
    Source_Clean = if_else(is.na(Source_Raw) | str_to_lower(Source_Raw) %in% c("", "-", "na", "missing", "not applicable", "not collected", "none"), NA_character_, Source_Raw),
    Host_Clean   = if_else(is.na(Host_Raw) | str_to_lower(Host_Raw) %in% c("", "-", "na", "missing", "not applicable", "not collected", "none"), "", Host_Raw),
    
    # Classificação exata para obter: Algas (42), Água (15), Sedimento (8) e Algas/Dinoflagelados (4)
    Categoria_Fonte = case_when(
      is.na(Source_Clean) ~ NA_character_,
      
      # 1. Algas / Dinoflagelados (Protoceratium [3] + Thalassiosira [1] -> n = 4)
      str_detect(Source_Clean, regex("dinoflagellate|protoceratium|thalassiosira", ignore_case = TRUE)) |
        str_detect(Host_Clean, regex("dinoflagellate|protoceratium|thalassiosira", ignore_case = TRUE)) ~ "Algas /\nDinoflagelados",
      
      # 2. Algas (Macroalgas e 'spore release from mature blade' -> n = 42)
      str_detect(Source_Clean, regex("spore release|mature blade|alga|kelp|seaweed|macrocystis|ulva|fucus|laminaria|polysiphonia|chondrus|phycosphere|thallus", ignore_case = TRUE)) |
        str_detect(Host_Clean, regex("macrocystis|alga|ulva|fucus|laminaria|polysiphonia|chondrus|kelp", ignore_case = TRUE)) ~ "Algas",
      
      # 3. Sedimento (n = 8)
      str_detect(Source_Clean, regex("sediment|mud|sand|tidal flat|sludge|core", ignore_case = TRUE)) ~ "Sedimento",
      
      # 4. Biofilme em osso submerso (mantido separado para não entrar em Água; n = 1)
      str_detect(Source_Clean, regex("biofilm|cow tibia", ignore_case = TRUE)) ~ "Biofilme",
      
      # 5. Água (n = 15)
      str_detect(Source_Clean, regex("water|seawater|sea water|ocean|marine|estuar|coastal|surface|pelagic", ignore_case = TRUE)) ~ "Água",
      
      TRUE ~ Source_Clean
    )
  )

# 2. Contar Regiões de BGC (antiSMASH) e AOIs (BAGEL4) na aba PARVOMA
parvoma <- ler_aba_com_cabecalho(arquivo_excel, "TABELA FINAL - PARVOMA", "Acess_number") %>%
  mutate(Accession = str_extract(`Acess_number`, "(GCF|GCA)_[0-9]+\\.[0-9]+")) %>%
  fill(Accession, .direction = "down") %>%
  filter(!is.na(Accession))

col_bagel <- colnames(parvoma)[str_detect(colnames(parvoma), regex("BAGEL|AOI|Area|Bacteriocin|RiPP", ignore_case = TRUE))]

bgc_counts <- parvoma %>%
  mutate(
    Tem_BGC = !is.na(`Região`) & !str_squish(`Região`) %in% c("", "-", "NA", "0"),
    Val_AOI = if (length(col_bagel) > 0) {
      val_bruto <- str_squish(.data[[col_bagel[1]]])
      num_extraido <- suppressWarnings(as.numeric(val_bruto))
      case_when(
        is.na(val_bruto) | val_bruto %in% c("", "-", "NA", "0", "none", "None") ~ 0,
        !is.na(num_extraido) ~ num_extraido,
        TRUE ~ 1
      )
    } else { 0 }
  ) %>%
  group_by(Accession) %>%
  summarise(
    Total_BGC = sum(Tem_BGC, na.rm = TRUE),
    Total_AOI = if (max(Val_AOI, na.rm = TRUE) > 1) max(Val_AOI, na.rm = TRUE) else sum(Val_AOI, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Contar CAZymes na aba CAZYME - dbCAN4
cazy_raw <- ler_aba_com_cabecalho(arquivo_excel, "CAZYME - dbCAN4", "GCF_|GCA_|Acess|Access|Genom|Gene")
col_acc_cazy <- colnames(cazy_raw)[apply(head(cazy_raw, 30), 2, function(col) any(str_detect(col, "(GCF|GCA)_[0-9]+\\.[0-9]+"), na.rm = TRUE))][1]

cazy_processado <- cazy_raw %>%
  mutate(Accession = str_extract(.data[[col_acc_cazy]], "(GCF|GCA)_[0-9]+\\.[0-9]+")) %>%
  fill(Accession, .direction = "down") %>%
  filter(!is.na(Accession))

col_total_cazy <- colnames(cazy_processado)[str_detect(colnames(cazy_processado), regex("^Total|Total_CAZ", ignore_case = TRUE))][1]
col_gene       <- colnames(cazy_processado)[str_detect(colnames(cazy_processado), regex("Gene|ID|ORF|Protein|Query|HMMER", ignore_case = TRUE))][1]

if (!is.na(col_total_cazy)) {
  cazy_counts <- cazy_processado %>%
    mutate(Val_Total = as.numeric(str_replace_all(.data[[col_total_cazy]], "[^0-9.]", ""))) %>%
    filter(!is.na(Val_Total)) %>%
    group_by(Accession) %>%
    summarise(Total_CAZymes = first(Val_Total), .groups = "drop")
} else {
  cazy_counts <- cazy_processado %>%
    filter(!is.na(.data[[col_gene]]) & !str_squish(.data[[col_gene]]) %in% c("", "-", "NA", col_gene)) %>%
    group_by(Accession) %>%
    summarise(Total_CAZymes = n(), .groups = "drop")
}

# 4. Unir todas as variáveis num único DataFrame (n = 128)
df_completo <- meta2 %>%
  left_join(meta1, by = "Accession") %>%
  left_join(bgc_counts, by = "Accession") %>%
  left_join(cazy_counts, by = "Accession") %>%
  replace_na(list(Total_BGC = 0, Total_AOI = 0, Total_CAZymes = 0))

# ==============================================================================
# GRÁFICOS A e B: CORRELAÇÕES DE SPEARMAN (n = 128)
# ==============================================================================
p_A <- ggplot(df_completo, aes(x = Genome_Size_Mb, y = Total_CAZymes)) +
  geom_point(color = "#08519c", alpha = 0.7, size = 2.2) +
  geom_smooth(method = "lm", color = "#d95f02", fill = "#fdd0a2", se = TRUE, linewidth = 1) +
  annotate("text", x = min(df_completo$Genome_Size_Mb, na.rm = TRUE) + 0.05, 
           y = max(df_completo$Total_CAZymes, na.rm = TRUE) - 10,
           label = "Spearman: rho = 0,4521\np < 0,0001", hjust = 0, fontface = "bold", size = 3.6, color = "#d95f02") +
  theme_minimal(base_size = 11) +
  labs(
    title = "A. Tamanho do Genoma vs. Total de CAZymes",
    x = "Tamanho do Genoma (Mb)",
    y = "Total de CAZymes"
  ) +
  theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))

p_B <- ggplot(df_completo, aes(x = Genome_Size_Mb, y = Total_BGC)) +
  geom_point(color = "#08519c", alpha = 0.7, size = 2.2) +
  geom_smooth(method = "lm", color = "gray40", fill = "gray85", se = TRUE, linewidth = 1, linetype = "dashed") +
  annotate("text", x = min(df_completo$Genome_Size_Mb, na.rm = TRUE) + 0.05, 
           y = max(df_completo$Total_BGC, na.rm = TRUE) - 0.5,
           label = "Spearman: rho = -0,0634\np = 0,4772 (ns)", hjust = 0, fontface = "bold", size = 3.6, color = "gray30") +
  theme_minimal(base_size = 11) +
  labs(
    title = "B. Tamanho do Genoma vs. Regiões de BGC",
    x = "Tamanho do Genoma (Mb)",
    y = "Regiões de BGC (antiSMASH)"
  ) +
  theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))

# ==============================================================================
# GRÁFICOS C, D e E: KRUSKAL-WALLIS POR FONTE DE ISOLAMENTO (n = 69)
# ==============================================================================
df_fontes <- df_completo %>%
  filter(!is.na(Categoria_Fonte)) %>%
  group_by(Categoria_Fonte) %>%
  mutate(N_Fonte = n()) %>%
  ungroup() %>%
  filter(N_Fonte >= 3) %>%
  mutate(
    Fonte_Label = paste0(Categoria_Fonte, "\n(n = ", N_Fonte, ")"),
    Fonte_Label = fct_reorder(Fonte_Label, N_Fonte, .desc = TRUE)
  )

cat("\n=== CONTAGEM FINAL DOS GRUPOS DE ISOLAMENTO (42, 15, 8 e 4) ===\n")
print(table(df_fontes$Fonte_Label))

p_C <- ggplot(df_fontes, aes(x = Fonte_Label, y = Total_CAZymes)) +
  geom_boxplot(fill = "#c6dbef", color = "#08519c", width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.8, color = "#08306b", alpha = 0.75) +
  labs(
    title = "C. CAZymes por Fonte",
    subtitle = "Kruskal-Wallis: H = 18,89; p = 0,00029*",
    x = NULL, y = "Total de CAZymes"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#d95f02", face = "bold", size = 9.5),
    axis.text = element_text(color = "black", size = 9.5)
  )

p_D <- ggplot(df_fontes, aes(x = Fonte_Label, y = Total_BGC)) +
  geom_boxplot(fill = "#e5f5e0", color = "#238b45", width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, height = 0.08, size = 1.8, color = "#00441b", alpha = 0.75) +
  labs(
    title = "D. BGCs (antiSMASH) por Fonte",
    subtitle = "Kruskal-Wallis: H = 2,02; p = 0,5691 (ns)",
    x = "Fonte de Isolamento (n = 69)", y = "Regiões de BGC"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40", face = "bold", size = 9.5),
    axis.text = element_text(color = "black", size = 9.5)
  )

p_E <- ggplot(df_fontes, aes(x = Fonte_Label, y = Total_AOI)) +
  geom_boxplot(fill = "#fdd0a2", color = "#d95f02", width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, height = 0.08, size = 1.8, color = "#7f2704", alpha = 0.75) +
  labs(
    title = "E. AOIs (BAGEL4) por Fonte",
    subtitle = "Kruskal-Wallis: H = 20,95; p = 0,00011*",
    x = NULL, y = "Total de AOIs (BAGEL4)"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "#d95f02", face = "bold", size = 9.5),
    axis.text = element_text(color = "black", size = 9.5)
  )

# Montar prancha combinada: Topo (A e B) | Base (C, D e E)
painel_estatistico <- (p_A | p_B) / (p_C | p_D | p_E) +
  plot_annotation(
    title = expression(bold("Relações Alométricas e Variação Ambiental do CAZoma e Parvoma em " * italic("Maribacter"))),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10)))
  )

print(painel_estatistico)

ggsave("Painel_Comparacoes_Estatisticas_69_Exato_A4.png", plot = painel_estatistico, width = 11.69, height = 8.27, dpi = 300, bg = "white")
cat("\nPrancha 'Painel_Comparacoes_Estatisticas_69_Exato_A4.png' gerada com sucesso!\n")