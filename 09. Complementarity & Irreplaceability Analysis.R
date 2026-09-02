# =========================================================
# 07 — COMPLEMENTARIDADE, SINGULARIDADE E INSUBSTITUIBILIDADE DESCRITIVA
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)
library(purrr)
library(tibble)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_comp <- file.path("outputs", "objetos", "composicao.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_comp)) stop("Execute primeiro 06_beta_diversidade_NMDS.R")
load(arquivo_base); load(arquivo_comp)

metas_representacao <- c(0.50, 0.75, 0.90, 0.95, 1.00)
n_repeticoes_selecao <- 100
set.seed(123)

dir_saida <- file.path("outputs", "complementaridade")
dir.create(dir_saida, recursive = TRUE, showWarnings = FALSE)

normalizar_01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  intervalo <- range(x, na.rm = TRUE)
  if (!all(is.finite(intervalo)) || diff(intervalo) == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - intervalo[1]) / diff(intervalo)
}

# COMPLEMENTARIDADE SEQUENCIAL
selecionar_complementaridade <- function( matriz_pa, desempate_aleatorio = FALSE) {
  
  matriz_pa <- as.matrix(matriz_pa)
  especies_totais <- colnames(matriz_pa)
  
  representadas <- setNames(
    rep(FALSE, length(especies_totais)), especies_totais)
  disponiveis <- seq_len(nrow(matriz_pa))
  resultados <- vector("list", nrow(matriz_pa))
  
  passo <- 0L
  
  while (
    length(disponiveis) > 0 && !all(representadas)) {
    
    ganho <- rowSums(matriz_pa[
        disponiveis, !representadas,
        drop = FALSE])
    
    maior_ganho <- max(ganho,na.rm = TRUE)
    
    if (
      !is.finite(maior_ganho) || maior_ganho <= 0) {break}
    
    candidatos <- disponiveis[ganho == maior_ganho]
    
    if (
      desempate_aleatorio &&
      length(candidatos) > 1) {
      escolhido <- sample(candidatos,1)
      
    } else {
      # Desempate determinístico:
      # 1. maior riqueza total;
      # 2. menor id de linha.
      riqueza_candidatos <- rowSums(
        matriz_pa[candidatos,,drop = FALSE])
      
      escolhido <- candidatos[
        order(-riqueza_candidatos,candidatos)][1]}
  
    passo <- passo + 1L
    
    novas_especies <- colnames(matriz_pa)[
      matriz_pa[escolhido,,
        drop = TRUE] > 0 & !representadas]
    representadas[novas_especies] <- TRUE
    
    resultados[[passo]] <- tibble::tibble(
      passo = passo,
      id_hex = as.integer(rownames(matriz_pa)[escolhido]),
      ganho_novas_especies = length(novas_especies),
      riqueza_hexagono = sum(
        matriz_pa[escolhido,,drop = TRUE]),
      n_especies_acumuladas = sum(representadas),
      proporcao_especies_acumuladas =
        sum(representadas) / length(representadas))
    disponiveis <- setdiff(disponiveis,escolhido)}
  dplyr::bind_rows(resultados[seq_len(passo)])
}

complementaridade <- selecionar_complementaridade(matriz_comunidade, desempate_aleatorio = FALSE)

print(utils::head(complementaridade, 20))

resumo_metas_complementaridade <- purrr::map_dfr(metas_representacao,
  function(meta) {
    linha_meta <- complementaridade |>
      dplyr::filter(proporcao_especies_acumuladas >= meta) |> dplyr::slice_head(
        n = 1)
    
    if (nrow(linha_meta) == 0) {
      return(
        tibble::tibble(
          meta = meta,
          n_hexagonos = NA_integer_,
          n_especies_representadas = NA_integer_))}
    
    tibble::tibble(
      meta = meta,
      n_hexagonos = linha_meta$passo,
      n_especies_representadas = linha_meta$n_especies_acumuladas)})

print(resumo_metas_complementaridade)

utils::write.csv(complementaridade,
  file.path(dir_saida, "ordem_complementaridade_hexagonos.csv"), row.names = FALSE, fileEncoding = "UTF-8")

utils::write.csv(resumo_metas_complementaridade,
  file.path(dir_saida, "hexagonos_necessarios_por_meta.csv"), row.names = FALSE, fileEncoding = "UTF-8")

grafico_complementaridade <- ggplot2::ggplot(
  complementaridade,
  ggplot2::aes(
    x = passo,
    y = proporcao_especies_acumuladas)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.3) +
  ggplot2::geom_hline(
    yintercept = metas_representacao,
    linetype = "dashed",
    linewidth = 0.3) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "A. Complementarity curve",
    subtitle = "Cumulative representation of species as cells are selected",
    x = "Number of selected hexagonal cells",
    y = "Proportion of species represented") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank())

print(grafico_complementaridade)

hex_complementaridade <- hex_50km |> dplyr::left_join(complementaridade, by = "id_hex")

mapa_complementaridade <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = sf::st_transform(brasil_proj, 4326),
    fill = "gray95",
    color = "gray50",
    linewidth = 0.25) +
  ggplot2::geom_sf(
    data = sf::st_transform(hex_complementaridade, 4326),
    ggplot2::aes(fill = passo),
    color = NA) +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    na.value = "transparent",
    name = "Selection\norder") +
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE) +
  ggplot2::geom_sf(
    data = estados_plot,
    fill = NA,
    color = "gray65",
    linewidth = 0.25) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "B. Sequential complementarity",
    subtitle = "Earlier cells add more previously unrepresented species",
    x = NULL,
    y = NULL) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank())
print(mapa_complementaridade)

painel_complementaridade <- grafico_complementaridade + mapa_complementaridade +
  patchwork::plot_layout(widths = c(1, 1.15))

print(painel_complementaridade)

ggplot2::ggsave(
  filename = file.path(dir_objetos, "painel_complementaridade.png"),
  plot = painel_complementaridade,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

# SINGULARIDADE E FREQUÊNCIA DE SELEÇÃO
n_hex_por_especie <- colSums(matriz_comunidade)
peso_restricao <- 1 / n_hex_por_especie

singularidade_hex <- tibble::tibble(
  id_hex = as.integer(rownames(matriz_comunidade)),
  riqueza = rowSums(matriz_comunidade),
  n_especies_exclusivas = rowSums(
    sweep(matriz_comunidade,2,
      n_hex_por_especie == 1,`*`)),
  indice_singularidade = as.numeric(matriz_comunidade %*% peso_restricao))

# Frequência de seleção em repetições do algoritmo de complementaridade.
frequencia_selecao <- purrr::map_dfr(seq_len(n_repeticoes_selecao),
  function(repeticao) {
    selecionar_complementaridade(matriz_comunidade,
      desempate_aleatorio = TRUE) |>
      dplyr::mutate(repeticao = repeticao)}) |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(frequencia_selecao = dplyr::n_distinct(
      repeticao) / n_repeticoes_selecao,
    passo_mediano = stats::median(
      passo, na.rm = TRUE),.groups = "drop")

indicadores_insubstituibilidade <- singularidade_hex |>
  dplyr::left_join(frequencia_selecao,
    by = "id_hex") |>
  dplyr::mutate(
    frequencia_selecao = dplyr::coalesce(frequencia_selecao,0),
    singularidade_norm = normalizar_01(indice_singularidade),
    exclusivas_norm = normalizar_01(n_especies_exclusivas),
    frequencia_norm = normalizar_01(frequencia_selecao),
    
    # Índice descritivo, não uma medida formal única de insubstituibilidade.
    indice_insubstituibilidade_descritivo =
      rowMeans(cbind(
          singularidade_norm,
          exclusivas_norm,
          frequencia_norm),na.rm = TRUE)) |>
  dplyr::arrange(dplyr::desc(indice_insubstituibilidade_descritivo))

hex_insubstituibilidade <- hex_50km |> dplyr::left_join(indicadores_insubstituibilidade, by = "id_hex")

mapa_insubstituibilidade <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj,4326),
    fill = "gray95",
    color = "gray50",
    linewidth = 0.25) +
  ggplot2::geom_sf(
    data = sf::st_transform(hex_insubstituibilidade, 4326),
    ggplot2::aes(fill = indice_insubstituibilidade_descritivo), color = NA) +
  ggplot2::scale_fill_viridis_c(
    na.value = "transparent",
    name = "Descriptive\nindex") +
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Cell singularity and selection frequency",
    subtitle = paste0(n_repeticoes_selecao,
      " complementarity runs with random tie-breaking"),
    x = NULL,
    y = NULL)

print(mapa_insubstituibilidade)

save(complementaridade, resumo_metas_complementaridade,
     singularidade_hex, frequencia_selecao,
     indicadores_insubstituibilidade, hex_insubstituibilidade,
     file = file.path("outputs", "objetos", "complementaridade_insubstituibilidade.RData"))
