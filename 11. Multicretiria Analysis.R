# =========================================================
# 09 — PRIORIZAÇÃO MULTICRITÉRIO DESCRITIVA
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)

necessarios <- c(
  file.path("outputs", "objetos", "base_analises.RData"),
  file.path("outputs", "objetos", "descritivos_hex.RData"),
  file.path("outputs", "objetos", "endemismo_hex.RData"),
  file.path("outputs", "objetos", "composicao.RData"),
  file.path("outputs", "objetos", "complementaridade_insubstituibilidade.RData"),
  file.path("outputs", "objetos", "conservacao.RData")
)
if (any(!file.exists(necessarios))) stop("Execute os scripts 00, 01, 02, 04, 06 e 07 antes deste.")
invisible(lapply(necessarios, load, envir = .GlobalEnv))

normalizar_01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  intervalo <- range(x, na.rm = TRUE)
  if (!all(is.finite(intervalo)) || diff(intervalo) == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - intervalo[1]) / diff(intervalo)
}

# PRIORIZAÇÃO MULTICRITÉRIO DESCRITIVA
# Esta etapa NÃO substitui a otimização do prioritizr. Serve para comparar espacialmente indicadores contínuos.
# Não usar pesos arbitrários sem justificar. Aqui usamos pesos iguais como cenário exploratório.

tabela_multicriterio <- hex_50km |>
  sf::st_drop_geometry() |>
  dplyr::select(id_hex) |>
  dplyr::left_join(
    hex_riqueza |>
      sf::st_drop_geometry() |>
      dplyr::select(id_hex, riqueza),
    by = "id_hex") |>
  dplyr::left_join(
    hex_endemismo |>
      sf::st_drop_geometry() |>
      dplyr::select(id_hex, WE, CWE),
    by = "id_hex") |>
  dplyr::left_join(
    beta_local |> dplyr::select(id_hex,
        beta_turnover_medio),by = "id_hex") |>
  dplyr::left_join(
    indicadores_insubstituibilidade |>
      dplyr::select(id_hex,
        indice_insubstituibilidade_descritivo),by = "id_hex")

if (
  exists("hex_protecao") && "proporcao_protegida" %in% names(hex_protecao)) {
  
  tabela_multicriterio <- tabela_multicriterio |>
    dplyr::left_join(hex_protecao |> sf::st_drop_geometry() |>
        dplyr::select(id_hex,proporcao_protegida),by = "id_hex")
  
} else {
  
  tabela_multicriterio$proporcao_protegida <- NA_real_}

tabela_multicriterio <- tabela_multicriterio |>
  dplyr::mutate(
    riqueza_norm = normalizar_01(riqueza),
    CWE_norm = normalizar_01(CWE),
    turnover_norm = normalizar_01(beta_turnover_medio),
    insubstituibilidade_norm = normalizar_01(indice_insubstituibilidade_descritivo),
    lacuna_protecao_norm =
      ifelse(is.na(proporcao_protegida),
        NA_real_,1 - normalizar_01(proporcao_protegida)),
    
    # Cenário exploratório com pesos iguais.
    prioridade_multicriterio =
      rowMeans(cbind(
          riqueza_norm,
          CWE_norm,
          turnover_norm,
          insubstituibilidade_norm,
          lacuna_protecao_norm),na.rm = TRUE),
    
    n_criterios_disponiveis =
      rowSums(!is.na(cbind(
            riqueza_norm,
            CWE_norm,
            turnover_norm,
            insubstituibilidade_norm,
            lacuna_protecao_norm)))) |>
  dplyr::filter(n_criterios_disponiveis >= 3)

hex_multicriterio <- hex_50km |>dplyr::left_join(tabela_multicriterio,by = "id_hex")

View(tabela_multicriterio)

hex_multicriterio <- hex_50km |>
  dplyr::left_join(
    tabela_multicriterio,
    by = "id_hex"
  )

# =========================================================
# TOP 20 HEXÁGONOS COM MAIOR PRIORIDADE MULTICRITÉRIO
# =========================================================

top_multicriterio <- hex_multicriterio |>
  dplyr::filter(!is.na(prioridade_multicriterio)) |>
  dplyr::slice_max(
    prioridade_multicriterio,
    n = 20,
    with_ties = FALSE
  ) |>
  dplyr::mutate(
    ranking = dplyr::row_number()
  )

# =========================================================
# COORDENADAS DOS CENTRÓIDES
# =========================================================

top_multicriterio <- top_multicriterio |>
  sf::st_transform(4326) |>
  dplyr::mutate(
    centr = sf::st_centroid(geometry),
    longitude = sf::st_coordinates(centr)[, 1],
    latitude = sf::st_coordinates(centr)[, 2]
  ) |>
  sf::st_drop_geometry() |>
  dplyr::select(
    ranking,
    id_hex,
    longitude,
    latitude,
    prioridade_multicriterio,
    riqueza,
    CWE,
    beta_turnover_medio,
    indice_insubstituibilidade_descritivo,
    proporcao_protegida
  ) |>
  dplyr::arrange(ranking)

View(top_multicriterio)

#utils::write.csv(
#  tabela_multicriterio,
#  file.path(dir_priorizacao,"priorizacao_multicriterio_exploratoria.csv"),
#  row.names = FALSE, fileEncoding = "UTF-8")

mapa_multicriterio <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = sf::st_transform(brasil_proj, 4326),
    fill = "gray95",
    color = "gray50",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf(
    data = estados_br,
    fill = NA,
    color = "gray30",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf(
    data = sf::st_transform(hex_multicriterio, 4326),
    ggplot2::aes(fill = prioridade_multicriterio),
    color = NA
  ) +
  ggplot2::scale_fill_viridis_c(
    na.value = "transparent",
    name = "Priority\nscore"
  ) +
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Multicriteria priority",
    subtitle = paste0(
      "Equal weights: richness, CWE, turnover, ",
      "singularity and protection gap"
    ),
    x = NULL,
    y = NULL
  )

print(mapa_multicriterio)

save(tabela_multicriterio, hex_multicriterio, mapa_multicriterio,
     file = file.path("outputs", "objetos", "priorizacao_multicriterio.RData"))
