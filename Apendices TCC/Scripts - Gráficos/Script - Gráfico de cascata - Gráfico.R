library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)

# 1. Carregar os dados
caminho_ficheiro <- "Maribacter - Parvome + Cazome.xlsx"
dados <- read_excel(caminho_ficheiro, sheet = "Metadados 2")

# 2. Preparar e ordenar os dados
dados_grafico <- dados %>%
  arrange(desc(`N50 (bp)`)) %>% 
  mutate(Amostra_ID = factor(`Assecion Number`, levels = `Assecion Number`))

# 3. Gráfico 1: N50 
grafico_n50 <- ggplot(dados_grafico, aes(x = Amostra_ID, y = `N50 (bp)`)) +
  geom_col(fill = "#2c7fb8", width = 1) +
  labs(title = "N50 por Amostra", x = NULL, y = "N50 (bp)") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

# 4. Gráfico 2: Contig Count 
grafico_contigs <- ggplot(dados_grafico, aes(x = Amostra_ID, y = `Contig Count`)) +
  geom_col(fill = "#e34a33", width = 1) +
  labs(title = "Número de Contigs por Amostra", x = "Amostras (Ordenadas pelo maior N50)", y = "Nº de Contigs") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

# 5. Unir os gráficos na vertical 
grafico_final <- grafico_n50 / grafico_contigs + 
  plot_annotation(
    title = "Perfil de Qualidade das Montagens - Maribacter",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

# 6. Exibir e Exportar
print(grafico_final)
ggsave("Graficos_Amostras_Cascata.png", plot = grafico_final, width = 10, height = 8, dpi = 300)
ggsave("Graficos_Amostras_Cascata.pdf", plot = grafico_final, width = 10, height = 8)