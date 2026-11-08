# 1. Cria uma matriz de Presença/Ausência e Contagem para todos os 128 genomas
analise_coocorrencia <- dados_heatmap %>%
  # Garante que todas as classes apareçam como colunas (preenchendo ausências com 0)
  pivot_wider(names_from = Tipo_Cluster, values_from = n, values_fill = 0) %>%
  # Seleciona as colunas de interesse (incluindo terpene para contexto completo)
  select(Acess_number, `terpene-precursor`, T3PKS, terpene) %>%
  # Cria categorias claras de presença/ausência
  mutate(
    Status_Precursor = if_else(`terpene-precursor` > 0, "Com terpene-precursor", "Sem terpene-precursor"),
    Status_T3PKS     = if_else(T3PKS > 0, "Com T3PKS", "Sem T3PKS")
  )

# 2. Tabela 1: Porcentagem global (sobre os 128 genomas) e por combinação exata de cópias
resumo_combinacoes <- analise_coocorrencia %>%
  count(Status_Precursor, Status_T3PKS, `terpene-precursor`, T3PKS, terpene, name = "Qtd_Genomas") %>%
  mutate(
    Pct_Total_128 = round((Qtd_Genomas / sum(Qtd_Genomas)) * 100, 2)
  ) %>%
  group_by(Status_Precursor) %>%
  mutate(
    # Calcula a % DENTRO do grupo (ex: Entre os que têm precursor, quantos % NÃO têm T3PKS?)
    Pct_Dentro_do_Grupo = round((Qtd_Genomas / sum(Qtd_Genomas)) * 100, 2)
  ) %>%
  ungroup()

print("=== RESUMO DETALHADO DAS COMBINAÇÕES ===")
print(resumo_combinacoes)

# 3. Tabela 2: Cruzamento direto 2x2 (Presença vs. Ausência)
tabela_cruzada <- analise_coocorrencia %>%
  count(Status_Precursor, Status_T3PKS, name = "Qtd_Genomas") %>%
  group_by(Status_Precursor) %>%
  mutate(
    Total_Grupo = sum(Qtd_Genomas),
    Porcentagem_no_Grupo = round((Qtd_Genomas / Total_Grupo) * 100, 2)
  )

print("=== PORCENTAGEM DE T3PKS DENTRO DE CADA GRUPO ===")
print(tabela_cruzada)

# 4.Teste Exato de Fisher
# Verifica se a ausência de T3PKS quando há terpene-precursor é estatisticamente significativa (p < 0.05)
matriz_fisher <- table(analise_coocorrencia$Status_Precursor, analise_coocorrencia$Status_T3PKS)
teste_fisher <- fisher.test(matriz_fisher)

cat("\n=== TESTE EXATO DE FISHER (Mutuamente Exclusivos?) ===\n")
print(matriz_fisher)
cat("p-valor:", format.pval(teste_fisher$p.value, digits = 4), "\n")