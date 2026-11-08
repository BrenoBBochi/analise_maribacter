library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ------------------------------------------------------------------------------
# 1. LEITURA DIRETA DA PLANILHA EXCEL 
# ------------------------------------------------------------------------------
arquivo_excel <- "Maribacter - Parvome + Cazome.xlsx"
nome_aba      <- "TABELA FINAL - PARVOMA"

dados_brutos <- read_excel(arquivo_excel, sheet = nome_aba, skip = 2)

# Nomeando apenas as colunas necessárias pela posição
names(dados_brutos)[1]  <- "Acess_number"
names(dados_brutos)[16] <- "Codigo_AOI_P" # Coluna P
names(dados_brutos)[19] <- "Classe_BAGEL" # Coluna S

# Preenchendo células mescladas do Acess_number em tabela provisória
dados_limpos <- dados_brutos %>%
  filter(!if_all(everything(), is.na)) %>%
  mutate(Acess_number = str_remove_all(Acess_number, "\\s+")) %>%
  fill(Acess_number, .direction = "down")

# Lista completa dos 128 genomas
todos_genomas <- dados_limpos %>%
  distinct(Acess_number)

# Extração das 171 AOIs válidas
aois_validas <- dados_limpos %>%
  filter(!is.na(Codigo_AOI_P), !Codigo_AOI_P %in% c("-", "0", "N/A", "")) %>%
  mutate(
    Codigo_AOI_P = str_trim(Codigo_AOI_P),
    Classe_BAGEL = str_to_title(str_trim(Classe_BAGEL))
  ) %>%
  select(Acess_number, Codigo_AOI_P, Classe_BAGEL) %>%
  distinct()

# ------------------------------------------------------------------------------
# 2. DADOS PARA O PAINEL A (Distribuição nos 128 Genomas)
# ------------------------------------------------------------------------------
dados_painel_A <- todos_genomas %>%
  left_join(
    aois_validas %>% count(Acess_number, name = "Qtd_AOIs"),
    by = "Acess_number"
  ) %>%
  mutate(Qtd_AOIs = replace_na(Qtd_AOIs, 0)) %>%
  count(Qtd_AOIs, name = "Qtd_Genomas") %>%
  mutate(
    Porcentagem = (Qtd_Genomas / sum(Qtd_Genomas)) * 100,
    Categoria_AOI = factor(
      paste0(Qtd_AOIs, if_else(Qtd_AOIs == 1, " AOI", " AOIs")),
      levels = c("0 AOIs", "1 AOI", "2 AOIs")
    ),
    Rotulo = paste0(Qtd_Genomas, "\n(", round(Porcentagem, 1), "%)")
  )

# ------------------------------------------------------------------------------
# 3. DADOS PARA O PAINEL B (Classes e Sobreposição nas 171 AOIs)
# ------------------------------------------------------------------------------
dados_painel_B <- aois_validas %>%
  count(Classe_BAGEL, name = "Qtd_AOIs") %>%
  mutate(
    Porcentagem = (Qtd_AOIs / sum(Qtd_AOIs)) * 100,
    Classe_Formatada = case_when(
      Classe_BAGEL == "Sactipeptides" ~ "Sactipeptides",
      Classe_BAGEL == "Lanthipeptide_class_i" ~ "Lanthipeptide\n(Class I)",
      TRUE ~ Classe_BAGEL
    ),
    Status_Deteccao = if_else(
      Classe_BAGEL == "Lanthipeptide_class_i",
      "Detectado por ambos (BAGEL4 + antiSMASH)",
      "Exclusivo do BAGEL4"
    ),
    Rotulo = paste0(Qtd_AOIs, " (", round(Porcentagem, 1), "%)")
  )

# ------------------------------------------------------------------------------
# 4. CONSTRUÇÃO DOS GRÁFICOS A e B
# ------------------------------------------------------------------------------
grafico_A <- ggplot(dados_painel_A, aes(x = Categoria_AOI, y = Qtd_Genomas, fill = Categoria_AOI)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.4) +
  geom_text(aes(label = Rotulo), vjust = -0.3, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = c("0 AOIs" = "#BDBDBD", "1 AOI" = "#4DBBD5", "2 AOIs" = "#00A087")) +
  scale_y_continuous(limits = c(0, 78), expand = c(0, 0)) +
  labs(
    title = "A. Distribuição de AOIs por genoma",
    subtitle = "Total = 128 genomas de Maribacter analisados",
    x = "Quantidade de AOIs preditas por genoma",
    y = "Número de genomas"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray30", size = 10),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold")
  )

grafico_B <- ggplot(dados_painel_B, aes(x = Qtd_AOIs, y = reorder(Classe_Formatada, Qtd_AOIs), fill = Status_Deteccao)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.4) +
  geom_text(aes(label = Rotulo), hjust = -0.1, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = c(
    "Exclusivo do BAGEL4" = "#3C5488",
    "Detectado por ambos (BAGEL4 + antiSMASH)" = "#E64B35"
  )) +
  scale_x_continuous(limits = c(0, 215), expand = c(0, 0)) +
  labs(
    title = "B. Classes de RiPPs e complementaridade",
    subtitle = "Total = 171 AOIs identificadas",
    x = "Número de Áreas de Interesse (AOIs)",
    y = "Classe predita (BAGEL4)",
    fill = "Detecção entre ferramentas:"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9.5),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray30", size = 10),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold")
  )

# ------------------------------------------------------------------------------
# 5. EXPORTAÇÃO DA FIGURA FINAL
# ------------------------------------------------------------------------------
painel_final_bagel <- grafico_A + grafico_B + plot_layout(widths = c(1, 1.15))

print(painel_final_bagel)

ggsave(
  filename = "Figura_Resultados_BAGEL4.png",
  plot     = painel_final_bagel,
  width    = 11,
  height   = 5.2,
  dpi      = 300,
  bg       = "white"
)