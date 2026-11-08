# ==============================================================================
# GRÁFICO DE AFINIDADE POR SUBSTRATO 
# ==============================================================================
library(tidyverse)
library(readxl)
library(ggplot2)
library(RColorBrewer)

arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"

# 1. Encontra automaticamente a aba que contém "Recommend Results" e "Substrate"
abas <- excel_sheets(arquivo_excel)
aba_correta <- NULL
for (aba in abas) {
  cols <- colnames(read_excel(arquivo_excel, sheet = aba, n_max = 1, .name_repair = "minimal"))
  if ("Recommend Results" %in% cols && "Substrate" %in% cols) {
    aba_correta <- aba
    break
  }
}

cat("-> A ler a folha:", aba_correta, "\n")
cazymes <- read_excel(arquivo_excel, sheet = aba_correta)

# 2. Limpeza da coluna "Substrate" 
dados_contagem <- cazymes %>%
  select(CAZyme_Bruta = `Recommend Results`, Substrato_Bruto = `Substrate`) %>%
  mutate(
    CAZyme_Bruta    = as.character(CAZyme_Bruta),
    Substrato_Bruto = as.character(Substrato_Bruto),
    ID_Linha        = row_number()
  ) %>%
  # Remove linhas sem CAZyme válida ou sem substrato
  filter(
    !is.na(CAZyme_Bruta), CAZyme_Bruta != "-", CAZyme_Bruta != "",
    !is.na(Substrato_Bruto), Substrato_Bruto != "-", Substrato_Bruto != ""
  ) %>%
  # Extrai a classe principal da CAZyme (GH, GT, PL, CE, AA, CBM)
  mutate(Classe = str_extract(CAZyme_Bruta, "(GH|GT|PL|CE|AA|CBM)")) %>%
  filter(!is.na(Classe)) %>%
  separate_rows(Substrato_Bruto, sep = "\\s*[;,\\|\\+]\\s*") %>%
  mutate(
    Substrato = str_trim(str_to_lower(Substrato_Bruto)),
    Substrato = str_to_sentence(Substrato)
  ) %>%
  # Remove os "-" que sobraram após separar casos como "-;host glycan"
  filter(!Substrato %in% c("-", "", "Na", "None", "N/a")) %>%
  distinct(ID_Linha, Classe, Substrato)

# 3. Ranking completo dos substratos mais frequentes
ranking_substratos <- dados_contagem %>%
  count(Substrato, sort = TRUE, name = "Total_Enzimas")

cat("\n=== RANKING DE SUBSTRATOS PREDITOS EM MARIBACTER ===\n")
print(ranking_substratos, n = 25)

# Seleciona os Top 20 para o gráfico
top20_substratos <- ranking_substratos %>%
  slice_head(n = 20)

dados_plot <- dados_contagem %>%
  filter(Substrato %in% top20_substratos$Substrato) %>%
  count(Substrato, Classe, name = "Quantidade") %>%
  left_join(top20_substratos, by = "Substrato") %>%
  mutate(Substrato = fct_reorder(Substrato, Total_Enzimas))

cores_base <- brewer.pal(n = 9, name = "Set1")
cores_base[6] <- "#8B4513"
classes_ordem <- c("AA", "CBM", "CE", "GH", "GT", "PL")
cores_classes <- setNames(cores_base[1:6], classes_ordem)

# 4. Gera o Gráfico de Barras Horizontais Empilhadas por Classe
p_substratos <- ggplot(dados_plot, aes(x = Quantidade, y = Substrato, fill = Classe)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.78) +
  geom_text(
    data = top20_substratos %>% mutate(Substrato = fct_reorder(Substrato, Total_Enzimas)),
    aes(x = Total_Enzimas, y = Substrato, label = Total_Enzimas),
    inherit.aes = FALSE,
    hjust = -0.2,
    size = 3.6,
    fontface = "bold",
    color = "gray20"
  ) +
  scale_fill_manual(values = cores_classes, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Substratos alvo preditos em relação as CAZymes de Maribacter",
    subtitle = "Os 20 substratos mais detectados pelo dbCAN-sub",
    x = "Número de CAZymes",
    y = "Substrato predito",
    fill = "CAZyme"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "gray15", size = 10.5),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_substratos)

ggsave("Figura_Substratos_Maribacter_dbCANsub.png", p_substratos, width = 10.5, height = 7, dpi = 300, bg = "white")
write_tsv(ranking_substratos, "Ranking_Substratos_Maribacter.tsv")