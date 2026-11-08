# ==============================================================================
# ANÁLISE ESTATÍSTICA ESTRITA DOS DADOS BRUTOS (SEM CATEGORIZAÇÃO ARTIFICIAL)
# ==============================================================================
library(tidyverse)
library(readxl)
library(rstatix)
library(writexl)

arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"

# Função para localizar a linha exata do cabeçalho sem perder dados
ler_aba_com_cabecalho <- function(arquivo, aba, coluna_ancora) {
  previa <- suppressMessages(read_excel(arquivo, sheet = aba, n_max = 15, col_names = FALSE, col_types = "text"))
  linha_cabecalho <- which(apply(previa, 1, function(linha) any(str_detect(str_squish(linha), fixed(coluna_ancora, ignore_case = TRUE)))))[1]
  if (is.na(linha_cabecalho)) linha_cabecalho <- 1
  
  df <- suppressMessages(read_excel(arquivo, sheet = aba, skip = linha_cabecalho - 1, col_types = "text"))
  colnames(df) <- str_squish(colnames(df))
  return(df)
}

# Função auxiliar para limpar strings vazias ou marcadores de ausência 
limpar_vazio <- function(x) {
  x_limpo <- str_squish(as.character(x))
  if_else(x_limpo %in% c("", "-", "NA", "N/A", "na", "n/a", "missing", "not applicable", "not collected"), NA_character_, x_limpo)
}

# ==============================================================================
# 1. LEITURA EXATA DAS 3 ABAS
# ==============================================================================

# --- 1.1 Aba: Metadados 2 ---
meta_raw <- ler_aba_com_cabecalho(arquivo_excel, "Metadados 2", "Assecion Number")

metadados <- meta_raw %>%
  transmute(
    Accession          = str_extract(`Assecion Number`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    GTDB_Species       = limpar_vazio(`GTDB species`),
    NCBI_Strain_ID     = limpar_vazio(`NCBI Strain ID`),
    Completeness       = as.numeric(str_replace(`CheckM2 Completeness (%)`, ",", ".")),
    Contamination      = as.numeric(str_replace(`CheckM2 Contamination (%)`, ",", ".")),
    Genome_Size_bp     = as.numeric(str_replace_all(`Genome Size (bp)`, "[^0-9\\.]", "")),
    Genome_Size_Mb     = Genome_Size_bp / 1e6,
    GC_Percent         = as.numeric(str_replace(`GC Percentage (%)`, ",", ".")),
    Assembly_Level     = limpar_vazio(`Assembly Level`),
    Contig_Count       = as.numeric(str_replace_all(`Contig Count`, "[^0-9]", "")),
    N50_bp             = as.numeric(str_replace_all(`N50 (bp)`, "[^0-9\\.]", "")),
    Genome_Category    = limpar_vazio(`Genome Category`),
    Coding_Density     = as.numeric(str_replace(`Coding Density (%)`, ",", ".")),
    CDS_Count          = as.numeric(str_replace_all(`CDS Count`, "[^0-9]", "")),
    Protein_Count      = as.numeric(str_replace_all(`Protein Count`, "[^0-9]", "")),
    Host_Meta2         = limpar_vazio(`Host`),
    Lat_Lon            = limpar_vazio(`Latitude and Longitude`),
    Country            = limpar_vazio(`Country`),
    Collection_Date    = limpar_vazio(`Collection Date`),
    LPSN_Species       = limpar_vazio(`LPSN species validation`)
  ) %>%
  filter(!is.na(Accession)) %>%
  distinct(Accession, .keep_all = TRUE)

cat("-> Genomas lidos em 'Metadados 2':", nrow(metadados), "\n")

# --- 1.2 Aba: TABELA FINAL - PARVOMA ---
parvoma_raw <- ler_aba_com_cabecalho(arquivo_excel, "TABELA FINAL - PARVOMA", "Acess_number")

parvoma_por_genoma <- parvoma_raw %>%
  mutate(Accession = str_extract(`Acess_number`, "(GCF|GCA)_[0-9]+\\.[0-9]+")) %>%
  # Preenche para baixo APENAS o código do genoma (célula mesclada na Coluna A)
  fill(Accession, .direction = "down") %>%
  filter(!is.na(Accession)) %>%
  group_by(Accession) %>%
  summarise(
    NCBI_Organism_Parvoma  = first(na.omit(limpar_vazio(`NCBI_organism_name`))),
    Isolation_Source       = first(na.omit(limpar_vazio(`Isolation Source`))),
    Host_Parvoma           = first(na.omit(limpar_vazio(`Host`))),
    Total_Regioes_BGC      = sum(!is.na(limpar_vazio(`Região`)) & limpar_vazio(`Região`) != "0"),
    Soma_Protoclusters_BGC = sum(as.numeric(str_squish(`Total`)), na.rm = TRUE),
    Total_AOIs_BAGEL4      = sum(!is.na(limpar_vazio(`AOI`)) & limpar_vazio(`AOI`) != "0"),
    .groups = "drop"
  )

cat("-> Genomas lidos em 'TABELA FINAL - PARVOMA':", nrow(parvoma_por_genoma), "\n")

# --- 1.3 Aba: CAZYME - dbCAN4 ---
cazy_raw <- ler_aba_com_cabecalho(arquivo_excel, "CAZYME - dbCAN4", "Recommend Results")

contagem_cazymes <- cazy_raw %>%
  mutate(
    Accession = str_extract(`Amostra`, "(GCF|GCA)_[0-9]+\\.[0-9]+"),
    Rec_Res   = limpar_vazio(`Recommend Results`)
  ) %>%
  filter(!is.na(Accession), !is.na(Rec_Res)) %>%
  group_by(Accession) %>%
  summarise(
    Total_CAZymes = n(),
    Total_GH      = sum(str_detect(Rec_Res, "GH[0-9]+")),
    Total_GT      = sum(str_detect(Rec_Res, "GT[0-9]+")),
    Total_PL      = sum(str_detect(Rec_Res, "PL[0-9]+")),
    Total_CE      = sum(str_detect(Rec_Res, "CE[0-9]+")),
    Total_AA      = sum(str_detect(Rec_Res, "AA[0-9]+")),
    Total_CBM     = sum(str_detect(Rec_Res, "CBM[0-9]+")),
    .groups = "drop"
  )

cat("-> Genomas lidos em 'CAZYME - dbCAN4':", nrow(contagem_cazymes), "\n\n")

# ==============================================================================
# 2. JUNÇÃO PURA DAS TABELAS 
# ==============================================================================
base_bruta <- metadados %>%
  left_join(parvoma_por_genoma, by = "Accession") %>%
  left_join(contagem_cazymes, by = "Accession")

# ==============================================================================
# 3. AUDITORIA DE PREENCHIMENTO DOS METADADOS BRUTOS
# ==============================================================================
cat("====================================================================\n")
cat("1. AUDITORIA DE COMPLETUDE DOS CAMPOS BRUTOS (N =", nrow(base_bruta), "GENOMAS)\n")
cat("====================================================================\n")

campos_auditoria <- c("Isolation_Source", "Host_Parvoma", "Host_Meta2", "Country", "Lat_Lon", "Collection_Date", "GTDB_Species", "LPSN_Species")

tabela_auditoria <- map_dfr(campos_auditoria, function(col) {
  n_total <- nrow(base_bruta)
  n_preenchido <- sum(!is.na(base_bruta[[col]]))
  n_vazio <- sum(is.na(base_bruta[[col]]))
  tibble(
    Coluna_Original  = col,
    Preenchidos_N    = n_preenchido,
    Preenchidos_Pct  = round((n_preenchido / n_total) * 100, 2),
    Vazios_N         = n_vazio,
    Vazios_Pct       = round((n_vazio / n_total) * 100, 2)
  )
})
print(tabela_auditoria, width = Inf)

cat("\n--- Valores Brutos Exatos encontrados em 'Isolation Source' (TABELA FINAL - PARVOMA) ---\n")
contagem_iso_bruto <- base_bruta %>%
  count(Isolation_Source, sort = TRUE, name = "N_Genomas") %>%
  mutate(Porcentagem = round((N_Genomas / nrow(base_bruta)) * 100, 2))
print(contagem_iso_bruto, n = 100)

cat("\n--- Valores Brutos Exatos encontrados em 'Host' (TABELA FINAL - PARVOMA) ---\n")
contagem_host_bruto <- base_bruta %>%
  count(Host_Parvoma, sort = TRUE, name = "N_Genomas") %>%
  mutate(Porcentagem = round((N_Genomas / nrow(base_bruta)) * 100, 2))
print(contagem_host_bruto, n = 100)

cat("\n--- Valores Brutos Exatos encontrados em 'Country' (Metadados 2) ---\n")
contagem_country_bruto <- base_bruta %>%
  count(Country, sort = TRUE, name = "N_Genomas") %>%
  mutate(Porcentagem = round((N_Genomas / nrow(base_bruta)) * 100, 2))
print(contagem_country_bruto, n = 100)

# ==============================================================================
# 4. TESTES DE NORMALIDADE (SHAPIRO-WILK) NAS VARIÁVEIS NUMÉRICAS ORIGINAIS
# ==============================================================================
cat("\n====================================================================\n")
cat("2. TESTES DE NORMALIDADE DE SHAPIRO-WILK (Variáveis Numéricas)\n")
cat("====================================================================\n")

variaveis_numericas <- c(
  "Genome_Size_bp", "GC_Percent", "Coding_Density", "CDS_Count", "Protein_Count",
  "Total_CAZymes", "Total_GH", "Total_GT", "Total_PL", "Total_CE", "Total_AA", "Total_CBM",
  "Total_Regioes_BGC", "Soma_Protoclusters_BGC", "Total_AOIs_BAGEL4"
)

tabela_normalidade <- map_dfr(variaveis_numericas, function(v) {
  vals <- na.omit(base_bruta[[v]])
  if (length(vals) >= 3 && length(unique(vals)) > 1) {
    sw <- shapiro.test(vals)
    tibble(
      Variavel     = v,
      N            = length(vals),
      Media        = round(mean(vals), 2),
      DP           = round(sd(vals), 2),
      Mediana      = round(median(vals), 2),
      IQR          = round(IQR(vals), 2),
      Min_Max      = paste0(round(min(vals), 2), " - ", round(max(vals), 2)),
      Shapiro_W    = round(unname(sw$statistic), 4),
      Shapiro_p    = sw$p.value,
      Distribuicao = if_else(sw$p.value > 0.05, "Normal (Paramétrica)", "Não-Normal (Não-Paramétrica)")
    )
  } else {
    tibble(Variavel = v, N = length(vals), Distribuicao = "Dados insuficientes/constantes")
  }
})
print(tabela_normalidade, width = Inf)

# ==============================================================================
# 5. TESTES DE CORRELAÇÃO (PEARSON OU SPEARMAN CONFORME SHAPIRO-WILK)
# ==============================================================================
cat("\n====================================================================\n")
cat("3. TESTES DE CORRELAÇÃO ENTRE VARIÁVEIS NUMÉRICAS\n")
cat("====================================================================\n")

pares_correlacao <- list(
  c("Genome_Size_bp",    "Total_CAZymes"),
  c("CDS_Count",         "Total_CAZymes"),
  c("Protein_Count",     "Total_CAZymes"),
  c("Genome_Size_bp",    "Total_Regioes_BGC"),
  c("Genome_Size_bp",    "Soma_Protoclusters_BGC"),
  c("Genome_Size_bp",    "Total_AOIs_BAGEL4"),
  c("CDS_Count",         "Total_Regioes_BGC"),
  c("Total_CAZymes",     "Total_Regioes_BGC"),
  c("Total_CAZymes",     "Soma_Protoclusters_BGC"),
  c("Total_CAZymes",     "Total_AOIs_BAGEL4"),
  c("GC_Percent",        "Total_CAZymes"),
  c("GC_Percent",        "Total_Regioes_BGC")
)

resultados_correlacao <- map_dfr(pares_correlacao, function(par) {
  var_x <- par[1]
  var_y <- par[2]
  df_par <- base_bruta %>% select(all_of(c(var_x, var_y))) %>% drop_na()
  
  norm_x <- shapiro.test(df_par[[var_x]])$p.value > 0.05
  norm_y <- shapiro.test(df_par[[var_y]])$p.value > 0.05
  metodo <- if (norm_x && norm_y) "pearson" else "spearman"
  
  teste <- suppressWarnings(cor.test(df_par[[var_x]], df_par[[var_y]], method = metodo))
  
  tibble(
    Variavel_X    = var_x,
    Variavel_Y    = var_y,
    N             = nrow(df_par),
    Metodo_Usado  = str_to_title(metodo),
    Coeficiente   = round(unname(teste$estimate), 4),
    p_valor       = teste$p.value,
    Significancia = case_when(
      teste$p.value < 0.001 ~ "*** (p < 0.001)",
      teste$p.value < 0.01  ~ "** (p < 0.01)",
      teste$p.value < 0.05  ~ "* (p < 0.05)",
      TRUE                  ~ "ns (Não significativo)"
    )
  )
})
print(resultados_correlacao, width = Inf)

# ==============================================================================
# 6. COMPARAÇÃO ENTRE GRUPOS BRUTOS (SEM AGRUPAMENTO ARTIFICIAL, N >= 3)
# ==============================================================================
cat("\n====================================================================\n")
cat("4. TESTES DE COMPARAÇÃO DE GRUPOS NAS COLUNAS ORIGINAIS (N >= 3 por grupo)\n")
cat("====================================================================\n")

testar_coluna_categorica_bruta <- function(df, col_grupo, vars_resposta) {
  df_valido <- df %>%
    filter(!is.na(.data[[col_grupo]])) %>%
    group_by(.data[[col_grupo]]) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  n_grupos <- n_distinct(df_valido[[col_grupo]])
  
  cat("\n--------------------------------------------------------------------\n")
  cat("Variável Agrupadora Bruta:", col_grupo, "| Grupos com N >= 3:", n_grupos, "\n")
  
  if (n_grupos < 2) {
    cat("-> Insuficiente para teste de hipótese: menos de 2 grupos possuem N >= 3 observações idênticas.\n")
    return(tibble())
  }
  
  print(df_valido %>% count(.data[[col_grupo]], sort = TRUE, name = "N_Amostras"))
  
  res_grupos <- map_dfr(vars_resposta, function(var_y) {
    df_sub <- df_valido %>% select(all_of(c(col_grupo, var_y))) %>% drop_na()
    if (n_distinct(df_sub[[col_grupo]]) < 2) return(tibble())
    
    formula_teste <- as.formula(paste(var_y, "~", col_grupo))
    
    modelo_lm <- lm(formula_teste, data = df_sub)
    p_shapiro_res <- tryCatch(shapiro.test(residuals(modelo_lm))$p.value, error = function(e) 0)
    p_levene <- tryCatch(df_sub %>% levene_test(formula_teste) %>% pull(p), error = function(e) 0)
    
    usa_parametrico <- (p_shapiro_res > 0.05) && (p_levene > 0.05)
    
    if (usa_parametrico) {
      teste_global <- df_sub %>% anova_test(formula_teste)
      p_global <- as.numeric(teste_global$p)
      est_global <- as.numeric(teste_global$F)
      nome_teste <- "ANOVA One-Way"
    } else {
      teste_global <- df_sub %>% kruskal_test(formula_teste)
      p_global <- as.numeric(teste_global$p)
      est_global <- as.numeric(teste_global$statistic)
      nome_teste <- "Kruskal-Wallis"
    }
    
    tibble(
      Coluna_Grupo       = col_grupo,
      Variavel_Resposta  = var_y,
      N_Total_Testado    = nrow(df_sub),
      N_Grupos_Testados  = n_distinct(df_sub[[col_grupo]]),
      Shapiro_Residuos_p = signif(p_shapiro_res, 3),
      Levene_p           = signif(p_levene, 3),
      Teste_Aplicado     = nome_teste,
      Estatistica        = round(est_global, 3),
      p_valor            = p_global,
      Significancia      = if_else(p_global < 0.05, "SIGNIFICATIVO (p < 0.05)", "Não significativo (p >= 0.05)")
    )
  })
  
  print(res_grupos, width = Inf)
  
  # Se houver resultado significativo, roda Post-hoc automaticamente
  vars_sig <- res_grupos %>% filter(p_valor < 0.05) %>% pull(Variavel_Resposta)
  if (length(vars_sig) > 0) {
    for (v_s in vars_sig) {
      cat("\n-> Post-hoc de Dunn (ajuste BH) para", v_s, "por", col_grupo, "(apenas pares com p.adj < 0.05):\n")
      f_post <- as.formula(paste(v_s, "~", col_grupo))
      pos <- df_valido %>% dunn_test(f_post, p.adjust.method = "BH") %>% filter(p.adj < 0.05)
      if (nrow(pos) > 0) print(pos, width = Inf) else cat("   Nenhum par isolado atingiu p.adj < 0.05 após correção para múltiplos testes.\n")
    }
  }
  
  return(res_grupos)
}

vars_alvo <- c("Total_CAZymes", "Total_Regioes_BGC", "Soma_Protoclusters_BGC", "Total_AOIs_BAGEL4", "Genome_Size_bp")

res_isolation <- testar_coluna_categorica_bruta(base_bruta, "Isolation_Source", vars_alvo)
res_host      <- testar_coluna_categorica_bruta(base_bruta, "Host_Parvoma", vars_alvo)
res_country   <- testar_coluna_categorica_bruta(base_bruta, "Country", vars_alvo)
res_gtdb      <- testar_coluna_categorica_bruta(base_bruta, "GTDB_Species", vars_alvo)
res_lpsn      <- testar_coluna_categorica_bruta(base_bruta, "LPSN_Species", vars_alvo)

todos_testes_grupos <- bind_rows(res_isolation, res_host, res_country, res_gtdb, res_lpsn)

# Exporta os resultados estritos
write_xlsx(
  list(
    "1_Base_Bruta_Unificada"  = base_bruta,
    "2_Auditoria_Metadados"   = tabela_auditoria,
    "3_Contagem_Iso_Source"   = contagem_iso_bruto,
    "4_Normalidade_Shapiro"   = tabela_normalidade,
    "5_Correlacoes"           = resultados_correlacao,
    "6_Testes_Grupos_Brutos"  = todos_testes_grupos
  ),
  "Relatorio_Estatistico_Bruto_Maribacter.xlsx"
)

cat("\n====================================================================\n")
cat("CONCLUÍDO! Arquivo 'Relatorio_Estatistico_Bruto_Maribacter.xlsx' gerado.\n")
cat("====================================================================\n")