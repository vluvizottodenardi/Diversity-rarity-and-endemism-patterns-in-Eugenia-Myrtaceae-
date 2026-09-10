# =========================================================
# 06 — BETA DIVERSIDADE E NMDS
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)
library(tibble)
library(betapart)
library(vegan)
library(patchwork)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
load(arquivo_base)

riqueza_minima_beta <- 5
registros_minimos_beta <- 5
n_hex_minimo_nmds <- 10

nmds_beta <- NULL
scores_nmds <- NULL
grafico_nmds <- NULL

# Função auxiliar para média de dissimilaridade por local
media_distancia_por_local <- function(objeto_distancia, ids) {
  matriz_distancia <- as.matrix(objeto_distancia)
  if (nrow(matriz_distancia) != length(ids)) stop("Comprimento de ids não corresponde à matriz de distâncias.")
  diag(matriz_distancia) <- NA_real_
  tibble::tibble(id_hex = ids, valor = rowMeans(matriz_distancia, na.rm = TRUE))
}

# MATRIZ HEXÁGONO × ESPÉCIE
pres_abs_analises <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(
    !is.na(id_hex),
    !is.na(scientificName_det),
    scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det)

matriz_comunidade <- stats::xtabs(~ id_hex + scientificName_det, data = pres_abs_analises)
matriz_comunidade <- ifelse(matriz_comunidade > 0, 1L, 0L)
matriz_comunidade <- as.matrix(matriz_comunidade)
ids_matriz <- suppressWarnings(as.integer(rownames(matriz_comunidade)))
riqueza_matriz <- rowSums(matriz_comunidade)

esforco_por_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex)) |>
  dplyr::count(id_hex, name = "n_registros")

hex_validos_beta <- tibble::tibble(
  id_hex = ids_matriz,
  riqueza = riqueza_matriz) |>
  dplyr::left_join(esforco_por_hex,
    by = "id_hex") |>
  dplyr::mutate(
    n_registros = dplyr::coalesce(n_registros, 0L)) |>
  dplyr::filter(
    riqueza >= riqueza_minima_beta,
    n_registros >= registros_minimos_beta)

matriz_beta <- matriz_comunidade[
  rownames(matriz_comunidade) %in% as.character(hex_validos_beta$id_hex),
  ,drop = FALSE]

# Remover espécies ausentes depois do filtro de células.
matriz_beta <- matriz_beta[, colSums(matriz_beta) > 0, drop = FALSE]
matriz_beta
View(matriz_beta)

# BETA DIVERSIDADE — TURNOVER, NESTEDNESS E TOTAL
beta_core <- betapart::betapart.core(matriz_beta)
beta_multi_jaccard <- betapart::beta.multi(beta_core, index.family = "jaccard")
beta_pair_jaccard <- betapart::beta.pair(beta_core, index.family = "jaccard")

# Componentes retornados:
# beta.jtu = turnover
# beta.jne = nestedness-resultant
# beta.jac = dissimilaridade total de Jaccard

print(beta_multi_jaccard)

# -------------------------------------------------------------------------
# INTERPRETAÇÃO DA BETA DIVERSIDADE MULTISSÍTIO (ÍNDICE DE JACCARD)
#
# beta.JAC = Beta diversidade total
#   Mede a diferença total na composição de espécies entre todos os hexágonos analisados.
#
# beta.JTU = Turnover (substituição de espécies)
#   Representa a parcela da beta diversidade causada pela substituição
#   de espécies entre os hexágonos. Valores altos indicam que diferentes
#   regiões possuem conjuntos distintos de espécies.
#
# beta.JNE = Nestedness-resultant
#   Representa a parcela da beta diversidade causada por perda de espécies
#   (subconjuntos). Valores altos indicam que alguns hexágonos possuem
#   apenas subconjuntos das espécies encontradas em áreas mais ricas.
# -------------------------------------------------------------------------

resumo_beta_multissitio <- tibble::tibble(
  componente = c("Turnover", "Nestedness-resultant","Jaccard total"),
  valor = c(
    unname(beta_multi_jaccard$beta.JTU),
    unname(beta_multi_jaccard$beta.JNE),
    unname(beta_multi_jaccard$beta.JAC)))

print(resumo_beta_multissitio)

ids_beta <- suppressWarnings(as.integer(rownames(matriz_beta)))

beta_turnover_local <- media_distancia_por_local(
  beta_pair_jaccard$beta.jtu, ids_beta) |>
  dplyr::rename(beta_turnover_medio = valor)

beta_nestedness_local <- media_distancia_por_local(
  beta_pair_jaccard$beta.jne,ids_beta) |>
  dplyr::rename(beta_nestedness_medio = valor)

beta_total_local <- media_distancia_por_local(
  beta_pair_jaccard$beta.jac,ids_beta) |>
  dplyr::rename(beta_jaccard_medio = valor)

beta_local <- beta_turnover_local |>
  dplyr::left_join(beta_nestedness_local, by = "id_hex") |>
  dplyr::left_join(beta_total_local, by = "id_hex") |>
  dplyr::left_join(hex_validos_beta, by = "id_hex")
View(beta_local)

hex_beta <- hex_50km |> dplyr::left_join(beta_local, by = "id_hex")
View(hex_beta)

mapa_beta_turnover <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj,4326),
    fill = "gray95",
    color = "gray50",
    linewidth = 0.25) +
  ggplot2::geom_sf(
    data = sf::st_transform(hex_beta,4326),
    mapping = ggplot2::aes(fill = beta_turnover_medio),color = NA) +
  ggplot2::scale_fill_viridis_c(
    na.value = "transparent",
    name = "Mean\nturnover") +
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "A. Compositional turnover",
    subtitle = "Mean pairwise Jaccard turnover per 50-km cell",
    x = NULL,
    y = NULL)
mapa_beta_turnover

mapa_beta_nestedness <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326),
    fill = "gray95",
    color = "gray50",
    linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_beta,4326),
    mapping = ggplot2::aes(fill = beta_nestedness_medio),color = NA) +
  ggplot2::scale_fill_viridis_c(
    na.value = "transparent",
    name = "Mean\nnestedness") +
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "B. Nestedness-resultant component",
    subtitle = "Mean pairwise Jaccard nestedness per 50-km cell",
    x = NULL,
    y = NULL)
mapa_beta_nestedness

mapa_beta_total <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj,4326),
    fill = "gray95",
    color = "gray50",
    linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_beta, 4326),
    mapping = ggplot2::aes(fill = beta_jaccard_medio), color = NA) +
  ggplot2::scale_fill_viridis_c(
    na.value = "transparent",
    name = "Mean\nJaccard") +
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "C. Total beta diversity",
    subtitle = "Mean pairwise Jaccard dissimilarity per 50-km cell",
    x = NULL,
    y = NULL)

painel_beta <- mapa_beta_turnover | mapa_beta_nestedness | mapa_beta_total
print(painel_beta)

ggplot2::ggsave(
  filename = file.path(dir_objetos, "painel_beta_diversidade.png"),
  plot = painel_beta,
  width = 21,
  height = 7,
  dpi = 300
)

# NMDS DA COMPOSIÇÃO
if (nrow(matriz_beta) >= n_hex_minimo_nmds) {
  set.seed(123)
  nmds_beta <- vegan::metaMDS(
    matriz_beta, distance = "jaccard",
    binary = TRUE,
    k = 2, trymax = 100,
    autotransform = FALSE,
    trace = FALSE)
  
  scores_nmds <- as.data.frame(
    vegan::scores(nmds_beta, display = "sites")) |>
    tibble::rownames_to_column(var = "id_hex") |>
    dplyr::mutate(id_hex = as.integer(id_hex)) |>
    dplyr::left_join(hex_validos_beta, by = "id_hex")
  
  grafico_nmds <- ggplot2::ggplot(
    scores_nmds, ggplot2::aes(
      x = NMDS1,
      y = NMDS2,
      size = riqueza,
      alpha = n_registros)) +
    ggplot2::geom_point() +
    ggplot2::scale_alpha_continuous(
      range = c(0.25, 0.9)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Floristic differentiation among hexagonal cells",
      subtitle = paste0("NMDS based on Jaccard dissimilarity; stress = ",
        round(nmds_beta$stress, 3)),
      x = "NMDS1",
      y = "NMDS2",
      size = "Richness",
      alpha = "Records")
  
  print(grafico_nmds)
}

save(matriz_comunidade, matriz_beta, hex_validos_beta,
     beta_multi_jaccard, beta_pair_jaccard, resumo_beta_multissitio,
     beta_local, hex_beta, nmds_beta, scores_nmds, grafico_nmds,
     file = file.path("outputs", "objetos", "composicao.RData"))

