# =========================================================
# 05 — COMPOSIÇÃO DE EUGENIA E ENDEMISMO POR HABITAT
#
# Perguntas:
#
# 1. How compositionally distinct are Eugenia assemblages
#    among vegetation domains?
#
# 2. Are areas of high Eugenia endemism associated with
#    particular vegetation domains?
#
# =========================================================
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(vegan)
library(tibble)
library(readr)
sf::sf_use_s2(FALSE)

# =========================================================
# 1. DIRETÓRIOS
# =========================================================

dir_objetos <- file.path("outputs", "objetos")
dir_tabelas <- file.path("outputs", "tabelas")
dir_vetores <- file.path("outputs", "vetores")

dir.create(
  dir_objetos,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  dir_tabelas,
  recursive = TRUE,
  showWarnings = FALSE)

dir.create(
  dir_vetores,
  recursive = TRUE,
  showWarnings = FALSE)

# =========================================================
# 2. CARREGAR A BASE
# =========================================================

load(
  file.path(
    dir_objetos,
    "base_analises.RData"))


# Conferência dos objetos
ls()

# =========================================================
# 3. CONFERIR AS COLUNAS
# =========================================================
names(dados_hex)
names(dados_sf)
names(hex_50km)

# Conferir habitats presentes
table(dados_sf$habitat, useNA = "ifany")

# =========================================================
# 4. HABITAT DOMINANTE EM CADA HEXÁGONO
# =========================================================
#
# O habitat atualmente está associado às ocorrências.
#
# Para comparar a composição entre habitats, precisamos
# atribuir um habitat ao HEXÁGONO.
#
# Aqui usamos:
#
#     habitat que ocupa a maior proporção da área do hexágono
#
# =========================================================
# ---------------------------------------------------------
# 4.1 Interseção entre hexágonos e vegetação
# ---------------------------------------------------------
# CORREÇÃO DE GEOMETRIAS ANTES DA INTERSEÇÃO
# Garantir geometria válida nos hexágonos
hex_50km_valid <- hex_50km |> sf::st_make_valid()

# Garantir geometria válida na vegetação
vegetacao_proj_valid <- vegetacao_proj |> sf::st_make_valid()

# Remover geometrias vazias
hex_50km_valid <- hex_50km_valid |> dplyr::filter(!sf::st_is_empty(geometry))
vegetacao_proj_valid <- vegetacao_proj_valid |> dplyr::filter(!sf::st_is_empty(geometry))

# Conferir validade
table(sf::st_is_valid(hex_50km_valid))
table(sf::st_is_valid(vegetacao_proj_valid))

# INTERSEÇÃO
hex_habitat <- sf::st_intersection(hex_50km_valid, vegetacao_proj_valid)

# ---------------------------------------------------------
# 4.2 Calcular área de cada fragmento de habitat
# ---------------------------------------------------------

hex_habitat <- hex_habitat |> dplyr::mutate(area_habitat_m2 = as.numeric(sf::st_area(geometry)))

# ---------------------------------------------------------
# 4.3 Área total de cada hexágono
# ---------------------------------------------------------

area_hex <- hex_50km |>
  dplyr::mutate(
    area_hex_m2 = as.numeric(
      sf::st_area(geometry))) |>
  sf::st_drop_geometry() |>
  dplyr::select(id_hex, area_hex_m2)

# ---------------------------------------------------------
# 4.4 Calcular proporção de cada habitat no hexágono
# ---------------------------------------------------------

hex_habitat <- hex_habitat |>
  sf::st_drop_geometry() |>
  dplyr::left_join(
    area_hex,
    by = "id_hex") |>
  dplyr::mutate(proporcao_habitat = area_habitat_m2 / area_hex_m2)

# ---------------------------------------------------------
# 4.5 Selecionar o habitat dominante
# ---------------------------------------------------------

habitat_dominante <- hex_habitat |>
  dplyr::group_by(id_hex) |>
  dplyr::slice_max(
    proporcao_habitat,
    n = 1,
    with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    id_hex,
    habitat,
    proporcao_habitat)


# Conferência
habitat_dominante |> dplyr::count(habitat, sort = TRUE)

# ---------------------------------------------------------
# 4.6 Adicionar habitat dominante aos hexágonos
# ---------------------------------------------------------

hex_50km_habitat <- hex_50km |>dplyr::left_join(habitat_dominante, by = "id_hex")

# =========================================================
# 5. SALVAR HEXÁGONOS COM HABITAT DOMINANTE
# =========================================================
sf::st_write(
  hex_50km_habitat,
  dsn = file.path(
    dir_vetores,
    "hexagonos_50km_habitat.gpkg"),
  layer = "hexagonos_habitat",
  delete_dsn = TRUE,
  quiet = TRUE
)

# =========================================================
# 6. CRIAR MATRIZ ESPÉCIE × HEXÁGONO
# =========================================================
# Cada célula:
#     1 = espécie presente
#     0 = espécie ausente
# =========================================================

dados_matriz <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(
    !is.na(id_hex),
    !is.na(scientificName_det)) |>
  dplyr::select(
    id_hex,
    scientificName_det) |>
  dplyr::distinct()

# ---------------------------------------------------------
# Criar matriz presença/ausência
# ---------------------------------------------------------

matriz_pa <- dados_matriz |>
  dplyr::mutate(presenca = 1) |>
  tidyr::pivot_wider(
    names_from = scientificName_det,
    values_from = presenca,
    values_fill = 0)

# Guardar IDs dos hexágonos
ids_hex <- matriz_pa$id_hex

# Transformar em matriz numérica
matriz_pa <- matriz_pa |>
  tibble::column_to_rownames("id_hex") |>
  as.matrix()

# Conferência
dim(matriz_pa)
head(matriz_pa[, seq_len(min(10, ncol(matriz_pa)))])

# =========================================================
# 7. ADICIONAR HABITAT AOS HEXÁGONOS DA MATRIZ
# =========================================================
dados_composicao <- habitat_dominante |>
  dplyr::filter(
    id_hex %in% rownames(matriz_pa)) |>
  dplyr::arrange(
    match(id_hex, rownames(matriz_pa)))

# Padronizar IDs
dados_composicao$id_hex <- as.character(dados_composicao$id_hex)
rownames(matriz_pa) <- as.character(rownames(matriz_pa))

# Identificar hexágonos sem habitat
hex_sem_habitat <- setdiff(rownames(matriz_pa), dados_composicao$id_hex)

cat("\nHexágonos na matriz:", nrow(matriz_pa),
    "\nHexágonos com habitat:", nrow(dados_composicao),
    "\nHexágonos sem habitat:", length(hex_sem_habitat),
    "\nID(s):", paste(hex_sem_habitat, collapse = ", "), "\n")

# Manter apenas IDs presentes em ambos
ids_validos <- intersect(rownames(matriz_pa), dados_composicao$id_hex)

# Reordenar matriz e tabela
matriz_pa <- matriz_pa[ids_validos, , drop = FALSE]

dados_composicao <- dados_composicao[
  match(ids_validos, dados_composicao$id_hex),
  ,
  drop = FALSE
]

# Conferência final
cat("\nLinhas matriz:", nrow(matriz_pa),
    "\nLinhas tabela:", nrow(dados_composicao),
    "\nMesma ordem:", identical(rownames(matriz_pa), dados_composicao$id_hex),
    "\n")

stopifnot(
  all(
    rownames(matriz_pa) ==
      dados_composicao$id_hex))

# =========================================================
# 8. REMOVER HEXÁGONOS SEM HABITAT DOMINANTE
# =========================================================
idx <- !is.na(dados_composicao$habitat)

matriz_pa <- matriz_pa[
  idx,
  ,
  drop = FALSE
]

dados_composicao <- dados_composicao[
  idx,
  ,
  drop = FALSE
]

# =========================================================
# 9. RIQUEZA DE ESPÉCIES POR HEXÁGONO
# =========================================================
dados_composicao$riqueza <- rowSums(matriz_pa)

# Resumo
resumo_riqueza_habitat <- dados_composicao |>
  dplyr::group_by(habitat) |>
  dplyr::summarise(
    n_hexagonos = dplyr::n(),
    riqueza_media = mean(riqueza),
    riqueza_mediana = median(riqueza),
    riqueza_min = min(riqueza),
    riqueza_max = max(riqueza),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(riqueza_media)
  )

resumo_riqueza_habitat

# =========================================================
# 10. GRÁFICO DE RIQUEZA POR HABITAT
# =========================================================

grafico_riqueza_habitat <- ggplot(dados_composicao, aes(x = habitat, y = riqueza, fill = habitat)) +
  geom_boxplot(width = 0.65, alpha = 0.85, color = "grey20", linewidth = 0.4, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.6) +
  scale_fill_viridis_d(option = "turbo", end = 0.9) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Species richness per 50-km hexagon") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title.y = element_text(size = 13), axis.line = element_line(linewidth = 0.5), plot.margin = margin(10, 15, 10, 10))
grafico_riqueza_habitat

# =========================================================
# 11. DISTÂNCIA DE JACCARD
# =========================================================
#
# Mede a diferença na composição de espécies entre os hexágonos.
#
# binary = TRUE:
# utiliza apenas presença/ausência.
#
# =========================================================

dist_jaccard <- vegan::vegdist(
  matriz_pa,
  method = "jaccard",
  binary = TRUE)

# =========================================================
# 12. PERMANOVA
# =========================================================
#
# Testa se a composição de Eugenia difere entre habitats.
#
# =========================================================

set.seed(123)

permanova <- vegan::adonis2(
  dist_jaccard ~ habitat,
  data = dados_composicao,
  permutations = 999
)

print(permanova)

# =========================================================
# 13. TABELA DA PERMANOVA
# =========================================================

permanova_tabela <- as.data.frame(permanova) |> tibble::rownames_to_column("termo")
permanova_tabela

# =========================================================
# 14. DISPERSÃO MULTIVARIADA
# =========================================================
#
# Verifica se os habitats apresentam diferentes níveis
# de heterogeneidade composicional.
#
# Isso é importante para interpretar corretamente a
# PERMANOVA.
#
# =========================================================

dispersion <- vegan::betadisper(dist_jaccard,dados_composicao$habitat)

# ANOVA
anova_disp <- anova(dispersion)
anova_disp

# Teste por permutação
set.seed(123)

permutacao_disp <- vegan::permutest(
  dispersion,
  permutations = 999
)

permutacao_disp

# =========================================================
# 15. DISTÂNCIA DOS HEXÁGONOS AO CENTROIDE
# =========================================================
distancias_centroide <- tibble::tibble(id_hex = names(dispersion$distances), distancia_centroide = as.numeric(dispersion$distances)) |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::left_join(dados_composicao |> dplyr::select(id_hex, habitat) |> dplyr::mutate(id_hex = as.character(id_hex)), by = "id_hex")
head(distancias_centroide)
table(is.na(distancias_centroide$habitat))
write.csv(distancias_centroide, file.path(dir_tabelas, "dispersao_composicional_Eugenia.csv"), row.names = FALSE)

# =========================================================
# 16. NMDS
# =========================================================
#
# Visualização da composição das assembleias.
#
# =========================================================
set.seed(123)

nmds <- vegan::metaMDS(
  matriz_pa,
  distance = "jaccard",
  binary = TRUE,
  k = 2,
  trymax = 100,
  autotransform = FALSE
)

print(nmds)

cat("\nStress do NMDS:",round(nmds$stress, 4),"\n")

# =========================================================
# 17. COORDENADAS DO NMDS
# =========================================================

nmds_sites <- as.data.frame(vegan::scores(nmds, display = "sites"))
nmds_sites$id_hex <- rownames(nmds_sites)

nmds_sites <- nmds_sites |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::left_join(
    dados_composicao |>
      dplyr::mutate(id_hex = as.character(id_hex)) |>
      dplyr::select(id_hex, habitat, riqueza),
    by = "id_hex"
  )
head(nmds_sites)
sum(is.na(nmds_sites$habitat)) # o normal é 0

# =========================================================
# 18. GRÁFICO NMDS
# =========================================================
hulls_nmds <- nmds_sites |>
  dplyr::filter(!is.na(habitat)) |>
  dplyr::group_by(habitat) |>
  dplyr::filter(dplyr::n() >= 3) |>
  dplyr::slice(chull(NMDS1, NMDS2)) |>
  dplyr::ungroup()
grafico_nmds <- ggplot() +
  geom_polygon(data = hulls_nmds, aes(x = NMDS1, y = NMDS2, group = habitat, fill = habitat), alpha = 0.15, color = NA) +
  geom_point(data = nmds_sites, aes(x = NMDS1, y = NMDS2, color = habitat), size = 2.5, alpha = 0.7) +
  scale_color_viridis_d(option = "turbo", end = 0.9) +
  scale_fill_viridis_d(option = "turbo", end = 0.9) +
  theme_classic(base_size = 13) +
  labs(x = "NMDS1", y = "NMDS2", color = "Vegetation domain", fill = "Vegetation domain") +
  theme(legend.title = element_text(face = "bold"), legend.position = "right", axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))
grafico_nmds

# =========================================================
# 19. ANÁLISE DE ENDEMISMO × HABITAT
# =========================================================
#
# IMPORTANTE:
#
# Esta parte utiliza o objeto hex_endemismo que você já
# produziu na análise de endemismo.
#
# Espera-se que ele contenha:
#
#     id_hex
#     WE
#     CWE
#
# =========================================================

# ---------------------------------------------------------
# Conferir se o objeto existe
# ---------------------------------------------------------

if (!exists("hex_endemismo")) {
  
  stop(
    paste0(
      "\nO objeto 'hex_endemismo' não foi encontrado.\n",
      "Carregue o objeto que contém WE e CWE antes de ",
      "executar esta seção.\n"
    )
  )
  
}


# Conferir colunas
names(hex_endemismo)

# =========================================================
# 20. ASSOCIAR HABITAT AO ENDEMISMO
# =========================================================

dados_endemismo <- hex_endemismo |>
  sf::st_drop_geometry() |>
  dplyr::select(
    id_hex,
    WE,
    CWE
  ) |>
  dplyr::inner_join(
    habitat_dominante,
    by = "id_hex"
  )

# Conferência
head(dados_endemismo)

# =========================================================
# 21. RESUMO DE WE POR HABITAT
# =========================================================

resumo_WE <- dados_endemismo |>
  dplyr::filter(
    !is.na(WE)
  ) |>
  dplyr::group_by(habitat) |>
  dplyr::summarise(
    n_hexagonos = dplyr::n(),
    media = mean(WE),
    mediana = median(WE),
    minimo = min(WE),
    maximo = max(WE),
    desvio_padrao = sd(WE),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(media)
  )


resumo_WE

# =========================================================
# 22. RESUMO DE CWE POR HABITAT
# =========================================================

resumo_CWE <- dados_endemismo |>
  dplyr::filter(
    !is.na(CWE)
  ) |>
  dplyr::group_by(habitat) |>
  dplyr::summarise(
    n_hexagonos = dplyr::n(),
    media = mean(CWE),
    mediana = median(CWE),
    minimo = min(CWE),
    maximo = max(CWE),
    desvio_padrao = sd(CWE),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(media)
  )


resumo_CWE

# =========================================================
# 23. TESTE DE CWE ENTRE HABITATS
# =========================================================

modelo_CWE <- lm(
  CWE ~ habitat,
  data = dados_endemismo)

anova_CWE <- anova(modelo_CWE)
anova_CWE

# =========================================================
# 24. KRUSKAL-WALLIS PARA CWE
# =========================================================
#
# Alternativa não paramétrica.
#
# Será especialmente útil se CWE não atender aos
# pressupostos do modelo linear.
#
# =========================================================

kruskal_CWE <- kruskal.test(
  CWE ~ habitat,
  data = dados_endemismo
)


kruskal_CWE

kruskal_WE <- kruskal.test(
  WE ~ habitat,
  data = dados_endemismo
)

kruskal_WE

library(FSA)

dunn_CWE <- FSA::dunnTest(
  CWE ~ habitat,
  data = dados_endemismo,
  method = "bh"
)

dunn_CWE$res

# =========================================================
# 25. GRÁFICO CWE × HABITAT
# =========================================================

grafico_CWE <- ggplot(
  dados_endemismo,
  aes(
    x = habitat,
    y = CWE
  )
) +
  geom_boxplot(
    outlier.alpha = 0.3
  ) +
  theme_classic() +
  labs(
    x = "Vegetation domain",
    y = "Corrected weighted endemism (CWE)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


grafico_CWE

# =========================================================
# 26. RELAÇÃO RIQUEZA × CWE
# =========================================================

riqueza_CWE <- dados_composicao |>
  dplyr::select(id_hex, riqueza) |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::inner_join(
    dados_endemismo |>
      dplyr::select(id_hex, habitat, CWE) |>
      dplyr::mutate(id_hex = as.character(id_hex)),
    by = "id_hex"
  ) |>
  dplyr::filter(!is.na(CWE))

# =========================================================
# 27. CORRELAÇÃO SPEARMAN
# =========================================================

cor_riqueza_CWE <- cor.test(
  riqueza_CWE$riqueza,
  riqueza_CWE$CWE,
  method = "spearman",
  exact = FALSE)

cor_riqueza_CWE

# =========================================================
# 28. GRÁFICO RIQUEZA × CWE
# =========================================================

grafico_riqueza_CWE <- ggplot(
  riqueza_CWE,
  aes(
    x = riqueza,
    y = CWE,
    fill = habitat
  )
) +
  geom_point(
    shape = 21,
    size = 3,
    alpha = 0.7
  ) +
  theme_classic() +
  labs(
    x = "Species richness per 50-km hexagon",
    y = "Corrected weighted endemism (CWE)",
    fill = "Vegetation domain"
  )

grafico_riqueza_CWE

# =========================================================
# 29. EXPORTAR TABELAS
# =========================================================

write.csv(
  habitat_dominante,
  file.path(
    dir_tabelas,
    "habitat_dominante_hexagonos.csv"
  ),
  row.names = FALSE
)


write.csv(
  resumo_riqueza_habitat,
  file.path(
    dir_tabelas,
    "riqueza_Eugenia_por_habitat.csv"
  ),
  row.names = FALSE
)


write.csv(
  permanova_tabela,
  file.path(
    dir_tabelas,
    "PERMANOVA_composicao_Eugenia.csv"
  ),
  row.names = FALSE
)


write.csv(
  distancias_centroide,
  file.path(
    dir_tabelas,
    "dispersao_composicional_Eugenia.csv"
  ),
  row.names = FALSE
)


write.csv(
  nmds_sites,
  file.path(
    dir_tabelas,
    "NMDS_composicao_Eugenia.csv"
  ),
  row.names = FALSE
)


write.csv(
  resumo_WE,
  file.path(
    dir_tabelas,
    "WE_por_habitat.csv"
  ),
  row.names = FALSE
)


write.csv(
  resumo_CWE,
  file.path(
    dir_tabelas,
    "CWE_por_habitat.csv"
  ),
  row.names = FALSE
)


write.csv(
  riqueza_CWE,
  file.path(
    dir_tabelas,
    "riqueza_CWE_por_hexagono.csv"
  ),
  row.names = FALSE
)


# =========================================================
# 30. EXPORTAR GRÁFICOS
# =========================================================

ggsave(
  file.path(
    dir_tabelas,
    "NMDS_composicao_Eugenia.png"
  ),
  grafico_nmds,
  width = 8,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    dir_tabelas,
    "riqueza_Eugenia_por_habitat.png"
  ),
  grafico_riqueza_habitat,
  width = 8,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    dir_tabelas,
    "CWE_Eugenia_por_habitat.png"
  ),
  grafico_CWE,
  width = 8,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    dir_tabelas,
    "riqueza_CWE_Eugenia.png"
  ),
  grafico_riqueza_CWE,
  width = 8,
  height = 6,
  dpi = 300
)


# =========================================================
# 31. SALVAR OBJETOS
# =========================================================

save(
  matriz_pa,
  dados_composicao,
  habitat_dominante,
  hex_50km_habitat,
  dist_jaccard,
  permanova,
  dispersion,
  nmds,
  dados_endemismo,
  riqueza_CWE,
  file = file.path(
    dir_objetos,
    "analises_composicao_endemismo.RData"
  )
)


# =========================================================
# 32. RESUMO FINAL NO CONSOLE
# =========================================================

cat(
  "\n==========================================================",
  "\nCOMPOSIÇÃO DE EUGENIA ENTRE HABITATS",
  "\n==========================================================\n"
)

cat(
  "\nNúmero de hexágonos:",
  nrow(matriz_pa),
  "\nNúmero de espécies:",
  ncol(matriz_pa),
  "\n"
)

cat(
  "\n--- PERMANOVA ---\n"
)

print(permanova)


cat(
  "\n--- DISPERSÃO MULTIVARIADA ---\n"
)

print(anova_disp)


cat(
  "\n--- NMDS ---\n"
)

cat(
  "Stress = ",
  round(nmds$stress, 4),
  "\n",
  sep = ""
)


cat(
  "\n--- CWE × HABITAT ---\n"
)

print(anova_CWE)


cat(
  "\n--- KRUSKAL-WALLIS CWE × HABITAT ---\n"
)

print(kruskal_CWE)


cat(
  "\n--- CORRELAÇÃO RIQUEZA × CWE ---\n"
)

print(cor_riqueza_CWE)


cat(
  "\n==========================================================",
  "\nANÁLISES FINALIZADAS",
  "\n==========================================================\n"
)


sf::sf_use_s2(TRUE)