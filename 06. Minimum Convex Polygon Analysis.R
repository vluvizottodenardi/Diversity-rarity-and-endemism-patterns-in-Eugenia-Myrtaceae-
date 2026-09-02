# =========================================================
# 05 — GAP ANALYSIS ADAPTADA DE RODRIGUES et al. (2004)
# POLÍGONOS DE DISTRIBUIÇÃO OBSERVADA (MCP), NÃO HEXÁGONOS
# =========================================================
# Esta versão substitui a soma de hexágonos ocupados por um polígono
# convexo mínimo (Minimum Convex Polygon; MCP) construído com os pontos
# de ocorrência de cada espécie.
#
# IMPORTANTE:
# - O MCP é uma aproximação da distribuição observada, não um mapa oficial
#   de distribuição, não é AOO e não deve ser chamado de EOO da IUCN.
# - Espécies com menos de 3 localidades espacialmente distintas, ou com
#   pontos colineares que não formem área, ficam "Não avaliadas" nesta análise.
# - A proteção por pontos continua disponível no script 04 e deve ser usada
#   em conjunto com esta análise.
# =========================================================
library(dplyr)
library(sf)
library(purrr)
library(ggplot2)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_cons <- file.path("outputs", "objetos", "conservacao.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_cons)) stop("Execute primeiro 04_conservacao_UCs_pontos_hotspots.R")
load(arquivo_base); load(arquivo_cons)

sf::sf_use_s2(FALSE)

dir_gap <- file.path("outputs", "gap_Rodrigues_poligonos")
dir_tabelas <- file.path(dir_gap, "tabelas")
dir_vetores <- file.path(dir_gap, "vetores")
dir_graficos <- file.path(dir_gap, "graficos")
dir.create(dir_tabelas, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_vetores, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_graficos, recursive = TRUE, showWarnings = FALSE)

# Limites de representação usados no esquema de Rodrigues et al. (2004).
limite_distribuicao_restrita_km2 <- 1000
limite_distribuicao_ampla_km2 <- 250000
meta_distribuicao_restrita <- 1.00
meta_distribuicao_ampla <- 0.10

# ---------------------------------------------------------
# 1. Preparar pontos e território em CRS de área
# ---------------------------------------------------------
crs_area <- sf::st_crs(brasil_proj)
pontos_proj <- dados_sf |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "", !sf::st_is_empty(geometry)) |>
  sf::st_transform(crs_area)

brasil_area <- brasil_proj |> sf::st_make_valid()
ucs_area <- ucs_union |> sf::st_transform(crs_area) |> sf::st_make_valid()

# ---------------------------------------------------------
# 2. Meta de proteção em função da área do polígono
# Interpolação log-linear entre 1.000 e 250.000 km².
# ---------------------------------------------------------
calcular_meta_rodrigues <- function(area_distribuicao_km2) {
  if (length(area_distribuicao_km2) != 1 || is.na(area_distribuicao_km2) ||
      !is.finite(area_distribuicao_km2) || area_distribuicao_km2 <= 0) return(NA_real_)
  if (area_distribuicao_km2 <= limite_distribuicao_restrita_km2) return(meta_distribuicao_restrita)
  if (area_distribuicao_km2 >= limite_distribuicao_ampla_km2) return(meta_distribuicao_ampla)

  posicao_log <- (
    log10(area_distribuicao_km2) - log10(limite_distribuicao_restrita_km2)
  ) / (
    log10(limite_distribuicao_ampla_km2) - log10(limite_distribuicao_restrita_km2))

  meta_distribuicao_restrita -
    posicao_log * (meta_distribuicao_restrita - meta_distribuicao_ampla)
}

# ---------------------------------------------------------
# 3. Construir um MCP por espécie
# ---------------------------------------------------------
criar_mcp_especie <- function(sp_name) {
  pts <- pontos_proj |> dplyr::filter(scientificName_det == sp_name)
  xy <- sf::st_coordinates(pts)
  if (nrow(xy) == 0) return(NULL)

  manter <- !duplicated(data.frame(X = xy[,1], Y = xy[,2]))
  pts <- pts[manter, ]
  n_pontos_unicos <- nrow(pts)

  if (n_pontos_unicos < 3) return(NULL)

  geom_mcp <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
  mcp <- sf::st_sf(
    scientificName_det = sp_name,
    n_pontos_unicos = n_pontos_unicos,
    geometry = geom_mcp
  ) |>
    sf::st_make_valid()

  # Restringir ao território brasileiro, pois a base analisada foi filtrada para o Brasil.
  mcp <- suppressWarnings(sf::st_intersection(mcp, brasil_area))
  if (nrow(mcp) == 0) return(NULL)

  area <- sum(as.numeric(sf::st_area(mcp)), na.rm = TRUE)
  if (!is.finite(area) || area <= 0) return(NULL)

  # Dissolver possíveis partes resultantes do recorte do Brasil.
  mcp |>
    dplyr::group_by(scientificName_det, n_pontos_unicos) |>
    dplyr::summarise(.groups = "drop")
}

especies <- sort(unique(pontos_proj$scientificName_det))
poligonos_lista <- lapply(especies, criar_mcp_especie)
poligonos_validos <- poligonos_lista[!vapply(poligonos_lista, is.null, logical(1))]

if (length(poligonos_validos) == 0) stop("Nenhuma espécie apresentou pontos suficientes para formar um polígono.")
poligonos_especies <- do.call(rbind, poligonos_validos) |> sf::st_make_valid()

# Número de localidades únicas para TODAS as espécies, inclusive as não avaliadas.
contagem_pontos <- lapply(especies, function(sp) {
  pts <- pontos_proj |> dplyr::filter(scientificName_det == sp)
  xy <- sf::st_coordinates(pts)
  n_unique <- if (nrow(xy) == 0) 0L else nrow(unique(data.frame(X = xy[,1], Y = xy[,2])))
  data.frame(scientificName_det = sp, n_pontos_unicos = n_unique)
}) |> dplyr::bind_rows()

# ---------------------------------------------------------
# 4. Área total do MCP e área sobreposta às UCs
# ---------------------------------------------------------
resultado_validos <- lapply(seq_len(nrow(poligonos_especies)), function(i) {
  poly <- poligonos_especies[i, ]
  sp_name <- poly$scientificName_det
  area_total_km2 <- as.numeric(sf::st_area(poly)) / 1e6

  tem_uc <- lengths(sf::st_intersects(poly, ucs_area)) > 0
  if (tem_uc) {
    inter <- suppressWarnings(sf::st_intersection(poly, ucs_area))
    area_protegida_km2 <- if (nrow(inter) > 0) sum(as.numeric(sf::st_area(inter)), na.rm = TRUE) / 1e6 else 0
  } else {
    area_protegida_km2 <- 0
  }

  area_protegida_km2 <- min(area_protegida_km2, area_total_km2)
  proporcao_protegida <- area_protegida_km2 / area_total_km2
  meta_protecao <- calcular_meta_rodrigues(area_total_km2)
  proporcao_meta_atendida <- proporcao_protegida / meta_protecao

  # Classificação mais próxima da lógica original de gap analysis:
  # gap = nenhuma representação; partial gap = alguma representação, mas meta não atingida;
  # protected/covered = meta atingida.
  categoria_gap <- dplyr::case_when(
    area_protegida_km2 <= sqrt(.Machine$double.eps) ~ "Não protegida (gap)",
    proporcao_protegida < meta_protecao ~ "Lacuna parcial",
    proporcao_protegida >= meta_protecao ~ "Protegida",
    TRUE ~ "Não avaliada"
  )

  data.frame(
    scientificName_det = sp_name,
    n_pontos_unicos = poly$n_pontos_unicos,
    area_poligono_km2 = area_total_km2,
    area_protegida_km2 = area_protegida_km2,
    proporcao_protegida = proporcao_protegida,
    percentual_protegido = 100 * proporcao_protegida,
    meta_protecao = meta_protecao,
    percentual_meta = 100 * meta_protecao,
    proporcao_meta_atendida = proporcao_meta_atendida,
    percentual_meta_atendida = 100 * proporcao_meta_atendida,
    categoria_gap = categoria_gap
  )
}) |> dplyr::bind_rows()

nao_avaliadas <- contagem_pontos |>
  dplyr::filter(!scientificName_det %in% resultado_validos$scientificName_det) |>
  dplyr::mutate(
    area_poligono_km2 = NA_real_,
    area_protegida_km2 = NA_real_,
    proporcao_protegida = NA_real_,
    percentual_protegido = NA_real_,
    meta_protecao = NA_real_,
    percentual_meta = NA_real_,
    proporcao_meta_atendida = NA_real_,
    percentual_meta_atendida = NA_real_,
    categoria_gap = "Não avaliada (<3 localidades distintas ou MCP sem área)"
  )

gap_rodrigues_poligonos <- dplyr::bind_rows(resultado_validos, nao_avaliadas) |>
  dplyr::arrange(categoria_gap, proporcao_meta_atendida)

resumo_gap <- gap_rodrigues_poligonos |>
  dplyr::count(categoria_gap, name = "n_especies") |>
  dplyr::mutate(percentual_especies = 100 * n_especies / sum(n_especies))

print(resumo_gap)

utils::write.csv(gap_rodrigues_poligonos, file.path(dir_tabelas, "gap_Rodrigues_poligonos_MCP.csv"), row.names = FALSE)
utils::write.csv(resumo_gap, file.path(dir_tabelas, "resumo_gap_Rodrigues_poligonos_MCP.csv"), row.names = FALSE)
sf::st_write(poligonos_especies, file.path(dir_vetores, "poligonos_distribuicao_MCP.gpkg"), layer = "MCP_especies", delete_layer = TRUE, quiet = TRUE)

# Gráfico-resumo
grafico_gap <- ggplot2::ggplot(resumo_gap, ggplot2::aes(x = reorder(categoria_gap, n_especies), y = n_especies)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Species representation in protected areas", subtitle = "Gap analysis using observed MCP range polygons", x = NULL, y = "Number of species")
print(grafico_gap)
ggplot2::ggsave(file.path(dir_graficos, "categorias_gap_Rodrigues_MCP.png"), grafico_gap, width = 8, height = 5, dpi = 300)

save(poligonos_especies, gap_rodrigues_poligonos, resumo_gap,
     file = file.path("outputs", "objetos", "gap_Rodrigues_poligonos.RData"))

sf::sf_use_s2(TRUE)

# =========================================================
# COMPARAÇÃO ENTRE PROTEÇÃO POR PONTOS E GAP POR MCP
# =========================================================

comparacao_protecao <- protecao_especies |>
  dplyr::select(
    scientificName_det,
    status_protecao_pontos
  ) |>
  dplyr::left_join(
    gap_rodrigues_poligonos |>
      dplyr::select(
        scientificName_det,
        categoria_gap
      ),
    by = "scientificName_det"
  )

# Conferir a tabela cruzada
tabela_comparacao <- comparacao_protecao |>
  dplyr::count(
    status_protecao_pontos,
    categoria_gap,
    name = "n_especies"
  )

print(tabela_comparacao)

# =========================================================
# FIGURA — PROPORÇÃO: PROTEÇÃO POR PONTOS × GAP POR MCP
# =========================================================

grafico_comparacao_protecao_prop <- ggplot2::ggplot(
  tabela_comparacao,
  ggplot2::aes(
    x = status_protecao_pontos,
    y = n_especies,
    fill = categoria_gap
  )
) +
  
  # Barras proporcionais
  ggplot2::geom_col(
    position = "fill",
    width = 0.8
  ) +
  
  # Número absoluto de espécies dentro de cada segmento
  ggplot2::geom_text(
    ggplot2::aes(
      label = n_especies
    ),
    position = ggplot2::position_fill(
      vjust = 0.5
    ),
    size = 3.5
  ) +
  
  # Eixo Y em porcentagem
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = c(0, 0)
  ) +
  
  # -------------------------------------------------------
# EDITAR NOMES E CORES DA LEGENDA
# -------------------------------------------------------

ggplot2::scale_fill_manual(
  values = c(
    "Protegida" = "#2E8B57",
    "Lacuna parcial" = "#E6AB02",
    "Não protegida (gap)" = "#D73027",
    "Não avaliada (<3 localidades distintas ou MCP sem área)" = "#BDBDBD"
  ),
  labels = c(
    "Protegida" = "Protected",
    "Lacuna parcial" = "Partial gap",
    "Não protegida (gap)" = "Unprotected gap",
    "Não avaliada (<3 localidades distintas ou MCP sem área)" = "Not evaluated"
  ),
  name = "MCP-based gap status"
) +
  
  ggplot2::theme_minimal() +
  
  ggplot2::labs(
    title = "Agreement between species protection assessments",
    subtitle = paste0(
      "Occurrence-based protection status versus ",
      "MCP-based gap analysis"
    ),
    x = "Occurrence-based protection status",
    y = "Proportion of species"
  ) +
  
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    legend.position = "right",
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 13
    ),
    axis.title = ggplot2::element_text(
      size = 11
    ),
    axis.text = ggplot2::element_text(
      size = 10
    ),
    legend.title = ggplot2::element_text(
      face = "bold"
    )
  )

print(grafico_comparacao_protecao_prop)

# =========================================================
# SALVAR FIGURA EM ALTA RESOLUÇÃO
# =========================================================
ggplot2::ggsave(
  filename = "protecao_pontos_vs_gap_MCP.png",
  plot = grafico_comparacao_protecao_prop,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

