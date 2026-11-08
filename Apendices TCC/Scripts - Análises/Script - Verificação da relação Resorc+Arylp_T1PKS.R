# 1. Matriz de presença/ausência para T1PKS, resorcinol e arylpolyene
analise_t1pks <- dados_heatmap %>%
  pivot_wider(names_from = Tipo_Cluster, values_from = n, values_fill = 0) %>%
  select(Acess_number, T1PKS, resorcinol, arylpolyene) %>%
  mutate(
    Status_T1PKS = if_else(T1PKS > 0, "Com T1PKS", "Sem T1PKS"),
    Status_Res   = if_else(resorcinol > 0, "Com resorcinol", "Sem resorcinol"),
    Status_APE   = if_else(arylpolyene > 0, "Com arylpolyene", "Sem arylpolyene"),
    
    # Categoria combinada de Resorcinol + Arylpolyene
    Perfil_Res_APE = case_when(
      resorcinol > 0 & arylpolyene > 0 ~ "Ambos (Resorcinol + Arylpolyene)",
      resorcinol > 0 & arylpolyene == 0 ~ "Apenas Resorcinol",
      resorcinol == 0 & arylpolyene > 0 ~ "Apenas Arylpolyene",
      TRUE ~ "Nenhum dos dois"
    )
  )

# 2. Tabela de Porcentagens: O que acontece quando NÃO há T1PKS vs. quando HÁ T1PKS?
resumo_t1pks <- analise_t1pks %>%
  count(Status_T1PKS, Perfil_Res_APE, name = "Qtd_Genomas") %>%
  group_by(Status_T1PKS) %>%
  mutate(
    Total_Grupo = sum(Qtd_Genomas),
    Porcentagem_no_Grupo = round((Qtd_Genomas / Total_Grupo) * 100, 2)
  ) %>%
  ungroup()

print("=== DISTRIBUIÇÃO DE RESORCINOL E ARYLPOLYENE POR STATUS DE T1PKS ===")
print(resumo_t1pks)

# 3. Testes Exatos de Fisher individuais e combinados
cat("\n=== 1. TESTE DE FISHER: T1PKS vs. RESORCINOL ===\n")
tab_res <- table(analise_t1pks$Status_T1PKS, analise_t1pks$Status_Res)
print(tab_res)
cat("p-valor (T1PKS x Resorcinol):", format.pval(fisher.test(tab_res)$p.value, digits = 4), "\n")

cat("\n=== 2. TESTE DE FISHER: T1PKS vs. ARYLPOLYENE ===\n")
tab_ape <- table(analise_t1pks$Status_T1PKS, analise_t1pks$Status_APE)
print(tab_ape)
cat("p-valor (T1PKS x Arylpolyene):", format.pval(fisher.test(tab_ape)$p.value, digits = 4), "\n")

cat("\n=== 3. TESTE DE FISHER: RESORCINOL vs. ARYLPOLYENE (Co-ocorrência) ===\n")
tab_co <- table(analise_t1pks$Status_Res, analise_t1pks$Status_APE)
print(tab_co)
cat("p-valor (Resorcinol x Arylpolyene):", format.pval(fisher.test(tab_co)$p.value, digits = 4), "\n")