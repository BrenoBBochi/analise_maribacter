# ==============================================================================
# CARACTERIZAÇÃO COMPLETA DOS METADADOS 
# ==============================================================================
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"

# Função para leitura segura do cabeçalho
ler_aba_com_cabecalho <- function(arquivo, aba, coluna_ancora) {
  previa <- suppressMessages(read_excel(arquivo, sheet = aba, n_max = 15, col_names = FALSE, col_types = "text"))
  linha_cabecalho <- which(apply(previa, 1, function(linha) any(str_detect(str_squish(linha), fixed(coluna_ancora, ignore_case = TRUE)))))[1]
  if (is.na(linha_cabecalho)) linha_cabecalho <- 1
  df <- suppressMessages(read_excel(arquivo, sheet = aba, skip = linha_cabecalho - 1, col_types = "text"))
  colnames(df) <- str_squish(colnames(df))
  return(df)
}

# 1. Ler Metadados 2 
meta_raw <- ler_aba_com_cabecalho(arquivo_excel, "Metadados 2", "Assecion Number")

col_gtdb   <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^GTDB", ignore_case = TRUE))][1]
col_comp   <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("CheckM2 Compl", ignore_case = TRUE))][1]
col_cont   <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("CheckM2 Conta", ignore_case = TRUE))][1]
col_gc     <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^GC Percentage", ignore_case = TRUE))][1]
col_ass    <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^Assembly Level", ignore_case = TRUE))][1]
col_contig <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^Contig Count", ignore_case = TRUE))][1]
col_cat    <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^Genome Categor", ignore_case = TRUE))][1]
col_source <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^Isolation source", ignore_case = TRUE))][1]
col_countr <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^Country", ignore_case = TRUE))][1]
col_reldat <- colnames(meta_raw)[str_detect(colnames(meta_raw), regex("^Release Date", ignore_case = TRUE))][1]

df_meta <- meta_raw %>%
  transmute(
    Accession        = str_extract(`Assecion Number`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    GTDB_Species     = str_squish(.data[[col_gtdb]]),
    Completeness     = as.numeric(str_replace_all(.data[[col_comp]], "[^0-9.]", "")),
    Contamination    = as.numeric(str_replace_all(.data[[col_cont]], "[^0-9.]", "")),
    GC_Percent       = as.numeric(str_replace_all(.data[[col_gc]], "[^0-9.]", "")),
    Contig_Count     = as.numeric(str_replace_all(.data[[col_contig]], "[^0-9.]", "")),
    Assembly_Level   = str_squish(.data[[col_ass]]),
    Genome_Category  = str_squish(.data[[col_cat]]),
    Isolation_Source = str_squish(.data[[col_source]]),
    Country_Raw      = str_squish(.data[[col_countr]]),
    Release_Date_Raw = str_squish(.data[[col_reldat]])
  ) %>%
  filter(!is.na(Accession)) %>%
  mutate(
    Assembly_Level   = if_else(is.na(Assembly_Level) | Assembly_Level %in% c("", "-", "NA"), "Não Informado", Assembly_Level),
    Genome_Category  = if_else(is.na(Genome_Category) | Genome_Category %in% c("", "-", "NA"), "Não Informado", Genome_Category),
    Isolation_Source = if_else(is.na(Isolation_Source) | Isolation_Source %in% c("", "-", "NA"), "Não Informado", Isolation_Source),
    Country_Clean    = if_else(is.na(Country_Raw) | Country_Raw %in% c("", "-", "NA", "missing", "not applicable"), 
                               "Não Informado", 
                               str_squish(str_remove(Country_Raw, ":.*$"))),
      Release_Year     = case_when(
      str_detect(Release_Date_Raw, "^[0-9]{5}(\\.[0-9]+)?$") ~ 
        as.numeric(format(as.Date(as.numeric(Release_Date_Raw), origin = "1899-12-30"), "%Y")),
      str_detect(Release_Date_Raw, "(19|20)[0-9]{2}") ~ 
        as.numeric(str_extract(Release_Date_Raw, "(19|20)[0-9]{2}")),
      TRUE ~ NA_real_
    )
  )



cat("\n=== DADOS PARA PREENCHER OS 'XXXXX' DO PARÁGRAFO ===\n")
cat("Total de Espécies GTDB:", n_distinct(df_meta$GTDB_Species), "\n")
cat("\n-> Menor Completude:\n")
print(df_meta %>% arrange(Completeness, Contamination) %>% select(Accession, GTDB_Species, Completeness, Contamination) %>% head(2))
cat("\n-> Maior Completude:\n")
print(df_meta %>% arrange(desc(Completeness), desc(Contamination)) %>% select(Accession, GTDB_Species, Completeness, Contamination) %>% head(2))
cat("\n-> Menor e Maior Contaminação:\n")
print(df_meta %>% arrange(Contamination) %>% select(Accession, GTDB_Species, Completeness, Contamination) %>% head(1))
print(df_meta %>% arrange(desc(Contamination)) %>% select(Accession, GTDB_Species, Completeness, Contamination) %>% head(1))

# ==============================================================================
# 2. RESUMOS ESTATÍSTICOS PARA TABELAS (COLUNAS G, H, I, K, O, R, S)
# ==============================================================================
resumo_continuas <- df_meta %>%
  pivot_longer(cols = c(GC_Percent, Contig_Count), names_to = "Variavel", values_to = "Valor") %>%
  group_by(Variavel) %>%
  summarise(
    N_Validos = sum(!is.na(Valor)),
    Media     = round(mean(Valor, na.rm = TRUE), 2),
    DP        = round(sd(Valor, na.rm = TRUE), 2),
    Mediana   = round(median(Valor, na.rm = TRUE), 2),
    IQR       = round(IQR(Valor, na.rm = TRUE), 2),
    Minimo    = round(min(Valor, na.rm = TRUE), 2),
    Maximo    = round(max(Valor, na.rm = TRUE), 2),
    Shapiro_p = signif(shapiro.test(Valor)$p.value, 4),
    .groups   = "drop"
  )

gerar_freq <- function(df, coluna) {
  df %>%
    count(.data[[coluna]], name = "Contagem") %>%
    mutate(Porcentagem = round((Contagem / sum(Contagem)) * 100, 2)) %>%
    arrange(desc(Contagem))
}

freq_assembly <- gerar_freq(df_meta, "Assembly_Level")
freq_category <- gerar_freq(df_meta, "Genome_Category")
freq_source   <- gerar_freq(df_meta, "Isolation_Source")
freq_country  <- gerar_freq(df_meta, "Country_Clean")
freq_year     <- df_meta %>%
  count(Release_Year, name = "Genomas_Ano") %>%
  arrange(Release_Year) %>%
  mutate(
    Porcentagem = round((Genomas_Ano / sum(Genomas_Ano)) * 100, 2),
    Acumulado   = cumsum(Genomas_Ano)
  )

cat("\n=== EVOLUÇÃO DE DEPÓSITOS POR ANO (CORRIGIDO) ===\n")
print(freq_year, n = 30)

write_xlsx(
  list(
    "G_I_Continuas"      = resumo_continuas,
    "H_Assembly_Level"   = freq_assembly,
    "K_Genome_Category"  = freq_category,
    "O_Isolation_Source" = freq_source,
    "R_Country_Limpo"    = freq_country,
    "S_Release_Year"     = freq_year
  ),
  "Resumo_Completo_Metadados_GHIKORS.xlsx"
)

# ==============================================================================
# 3. CONSTRUÇÃO DO PAINEL VISUAL MULTIPARTES (A, B, C, D)
# ==============================================================================

# Painel A: Fragmentação por Nível de Montagem
p_A <- ggplot(df_meta, aes(x = reorder(Assembly_Level, Contig_Count, FUN = median), y = Contig_Count)) +
  geom_boxplot(fill = "#c6dbef", color = "#08519c", width = 0.45, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 1.8, color = "#08306b", alpha = 0.75) +
  scale_y_log10() +
  theme_minimal(base_size = 11) +
  labs(
    title = "A. Fragmentação por Nível de Montagem",
    x = "Nível de Montagem (Assembly Level)",
    y = "Nº de Contigs (Escala Log10)"
  ) +
  theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))

# Painel B: Distribuição do Conteúdo GC (%)
med_gc <- median(df_meta$GC_Percent, na.rm = TRUE)
p_B <- ggplot(df_meta, aes(x = GC_Percent)) +
  geom_histogram(aes(y = after_stat(density)), bins = 18, fill = "#6baed6", color = "white", alpha = 0.85) +
  geom_density(color = "#08306b", linewidth = 0.9) +
  geom_vline(xintercept = med_gc, linetype = "dashed", color = "#d95f02", linewidth = 0.9) +
  annotate("text", x = med_gc + 0.3, y = 0.25, label = paste0("Mediana = ", round(med_gc, 2), "%"),
           color = "#d95f02", fontface = "bold", hjust = 0, size = 3.5) +
  theme_minimal(base_size = 11) +
  labs(
    title = "B. Distribuição do Conteúdo GC (%)",
    x = "Conteúdo GC (%)",
    y = "Densidade"
  ) +
  theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"))

# Painel C: Origem Geográfica
df_country_top <- df_meta %>%
  mutate(Country_Plot = fct_lump_n(Country_Clean, n = 8, other_level = "Outros Países/Regiões")) %>%
  count(Country_Plot, name = "Contagem") %>%
  mutate(Pct = round((Contagem / sum(Contagem)) * 100, 1))

p_C <- ggplot(df_country_top, aes(x = Contagem, y = reorder(Country_Plot, Contagem))) +
  geom_col(fill = "#2171b5", width = 0.65) +
  geom_text(aes(label = paste0(Contagem, " (", Pct, "%)")), hjust = -0.1, size = 3.3, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  theme_minimal(base_size = 11) +
  labs(
    title = "C. Distribuição Geográfica (Country)",
    x = "Número de Genomas",
    y = NULL
  ) +
  theme(plot.title = element_text(face = "bold"), axis.text = element_text(color = "black"), panel.grid.minor = element_blank())

# Painel D: Evolução Temporal de Depósitos (Anos Reais Corrigidos)
p_D <- ggplot(freq_year %>% filter(!is.na(Release_Year)), aes(x = factor(Release_Year), y = Genomas_Ano)) +
  geom_col(fill = "#4292c6", width = 0.7) +
  geom_text(aes(label = Genomas_Ano), vjust = -0.4, size = 3.5, fontface = "bold", color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  theme_minimal(base_size = 11) +
  labs(
    title = "D. Evolução Anual de Depósitos no NCBI",
    x = "Ano de Liberação (Release Year)",
    y = "Genomas Depositados"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

# Combina os 4 gráficos em uma única prancha A4 Paisagem
painel_metadados <- (p_A | p_B) / (p_C | p_D)
print(painel_metadados)

ggsave("Painel_Metadados_GHIKORS_A4.png", plot = painel_metadados, width = 11.69, height = 8.27, dpi = 300, bg = "white")
cat("\nPrancha 'Painel_Metadados_GHIKORS_A4.png' e planilha 'Resumo_Completo_Metadados_GHIKORS.xlsx' geradas com sucesso!\n")