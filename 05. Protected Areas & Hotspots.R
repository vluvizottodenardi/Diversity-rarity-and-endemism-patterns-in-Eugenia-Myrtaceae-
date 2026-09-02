# =========================================================
# 04 — CONSERVAÇÃO: UCs, HOTSPOTS E PROTEÇÃO POR PONTOS
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_desc <- file.path("outputs", "objetos", "descritivos_hex.RData")
arquivo_end <- file.path("outputs", "objetos", "endemismo_hex.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_desc)) stop("Execute primeiro 01_descritivos_esforco_riqueza.R")
if (!file.exists(arquivo_end)) stop("Execute primeiro 02_endemismo_WE_CWE.R")
load(arquivo_base); load(arquivo_desc); load(arquivo_end)

sf::sf_use_s2(FALSE)

arquivo_ucs <- "C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/Uc_Brazil/UC_Brazil.shp"
percentil_hotspot <- 0.90
# 100% = célula integralmente coberta. Isso NÃO é usado para dizer que uma espécie está protegida.
limiar_protecao_hex <- 1
# Mantém células com ao menos uma espécie; interpretar CWE de células muito pobres com cautela.
riqueza_minima_cwe <- 1

dir_conservacao <- file.path("outputs", "conservacao")
dir_tabelas <- file.path(dir_conservacao, "tabelas")
dir_graficos <- file.path(dir_conservacao, "graficos")
dir_vetores <- file.path(dir_conservacao, "vetores")
dir.create(dir_tabelas, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_graficos, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_vetores, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# UCs
# ---------------------------------------------------------
ucs <- sf::st_read(arquivo_ucs, quiet = TRUE) |>
  dplyr::filter(
    ISO3 == "BRA",
    REALM != "Marine",
    STATUS %in% c("Designated", "Established", "Inscribed")
  ) |>
  sf::st_make_valid() |>
  dplyr::filter(!sf::st_is_empty(geometry))

ucs_hex <- ucs |> sf::st_transform(sf::st_crs(hex_50km))
ucs_union <- ucs_hex |> sf::st_union() |> sf::st_make_valid()
ucs_union <- sf::st_sf(id_uc_union = 1, geometry = ucs_union)

# ---------------------------------------------------------
# Proteção espacial dos hexágonos: riqueza, WE e CWE
# ---------------------------------------------------------
hex_protecao <- hex_50km |>
  sf::st_make_valid() |>
  dplyr::mutate(area_hex_km2 = as.numeric(sf::st_area(geometry)) / 1e6)

indices_com_uc <- lengths(sf::st_intersects(hex_protecao, ucs_union)) > 0
hex_com_uc <- hex_protecao[indices_com_uc, ]

if (nrow(hex_com_uc) > 0) {
  intersecao_hex_uc <- suppressWarnings(sf::st_intersection(hex_com_uc |> dplyr::select(id_hex), ucs_union))
  area_protegida_hex <- intersecao_hex_uc |>
    dplyr::mutate(area_protegida_km2 = as.numeric(sf::st_area(geometry)) / 1e6) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(id_hex) |>
    dplyr::summarise(area_protegida_km2 = sum(area_protegida_km2, na.rm = TRUE), .groups = "drop")
} else {
  area_protegida_hex <- data.frame(id_hex = integer(0), area_protegida_km2 = numeric(0))
}

hex_protecao <- hex_protecao |>
  dplyr::left_join(area_protegida_hex, by = "id_hex") |>
  dplyr::mutate(
    area_protegida_km2 = dplyr::coalesce(area_protegida_km2, 0),
    area_protegida_km2 = pmin(area_protegida_km2, area_hex_km2),
    proporcao_protegida = area_protegida_km2 / area_hex_km2,
    percentual_protegido = 100 * proporcao_protegida,
    protegido = proporcao_protegida >= limiar_protecao_hex,
    completamente_fora_uc = proporcao_protegida == 0
  )

gap_riqueza <- hex_riqueza |>
  dplyr::select(id_hex, riqueza) |>
  dplyr::left_join(hex_protecao |> sf::st_drop_geometry() |> dplyr::select(id_hex, proporcao_protegida, percentual_protegido, protegido, completamente_fora_uc), by = "id_hex")

gap_endemismo <- hex_endemismo |>
  dplyr::select(id_hex, riqueza, WE, CWE) |>
  dplyr::left_join(hex_protecao |> sf::st_drop_geometry() |> dplyr::select(id_hex, proporcao_protegida, percentual_protegido, protegido, completamente_fora_uc), by = "id_hex")

limiar_riqueza <- stats::quantile(gap_riqueza$riqueza, percentil_hotspot, na.rm = TRUE, names = FALSE)
hotspots_riqueza <- gap_riqueza |>
  dplyr::filter(!is.na(riqueza), riqueza >= limiar_riqueza, completamente_fora_uc)

limiar_cwe <- stats::quantile(
  gap_endemismo$CWE[!is.na(gap_endemismo$CWE) & !is.na(gap_endemismo$riqueza) & gap_endemismo$riqueza >= riqueza_minima_cwe],
  percentil_hotspot, na.rm = TRUE, names = FALSE)
hotspots_cwe <- gap_endemismo |>
  dplyr::filter(!is.na(CWE), riqueza >= riqueza_minima_cwe, CWE >= limiar_cwe, completamente_fora_uc)

# OBS.: hotspots acima são "100% fora de UCs". A variável protegido (100% coberto)
# pode ser usada separadamente para comparações conservadoras de células.

# ---------------------------------------------------------
# Proteção das espécies baseada NOS PONTOS DE OCORRÊNCIA
# ---------------------------------------------------------
dados_uc <- dados_sf |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "", !sf::st_is_empty(geometry)) |>
  sf::st_transform(sf::st_crs(ucs_hex)) |>
  dplyr::mutate(registro_em_uc = lengths(sf::st_intersects(geometry, ucs_union)) > 0)

protecao_especies <- dados_uc |>
  sf::st_drop_geometry() |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(
    n_registros = dplyr::n(),
    n_registros_em_uc = sum(registro_em_uc, na.rm = TRUE),
    proporcao_registros_em_uc = n_registros_em_uc / n_registros,
    percentual_registros_em_uc = 100 * proporcao_registros_em_uc,
    .groups = "drop") |>
  dplyr::mutate(
    status_protecao_pontos = dplyr::case_when(
      n_registros_em_uc == 0 ~ "UPS - Unprotected",
      n_registros_em_uc == n_registros ~ "PS - Protected",
      TRUE ~ "PPS - Partially protected"
    ),
    pelo_menos_um_registro_em_uc = n_registros_em_uc > 0,
    todos_registros_em_uc = n_registros_em_uc == n_registros
  )

utils::write.csv(protecao_especies, file.path(dir_tabelas, "protecao_especies_por_pontos.csv"), row.names = FALSE)
utils::write.csv(sf::st_drop_geometry(hotspots_riqueza), file.path(dir_tabelas, "hotspots_riqueza_100pct_fora_UCs.csv"), row.names = FALSE)
utils::write.csv(sf::st_drop_geometry(hotspots_cwe), file.path(dir_tabelas, "hotspots_CWE_100pct_fora_UCs.csv"), row.names = FALSE)

sf::st_write(hotspots_riqueza, file.path(dir_vetores, "hotspots_conservacao.gpkg"), layer = "hotspots_riqueza_fora", delete_layer = TRUE, quiet = TRUE)
sf::st_write(hotspots_cwe, file.path(dir_vetores, "hotspots_conservacao.gpkg"), layer = "hotspots_CWE_fora", delete_layer = TRUE, quiet = TRUE)

save(ucs, ucs_hex, ucs_union, hex_protecao, gap_riqueza, gap_endemismo,
     hotspots_riqueza, hotspots_cwe, protecao_especies,
     file = file.path("outputs", "objetos", "conservacao.RData"))

sf::sf_use_s2(TRUE)

# ---------------------------------------------------------
# ESPÉCIES PRESENTES NOS HOTSPOTS
# ---------------------------------------------------------
# IDs dos hotspots
ids_hotspots_riqueza <- hotspots_riqueza$id_hex
ids_hotspots_cwe <- hotspots_cwe$id_hex

# ---------------------------------------------------------
# 1. Espécies nos hotspots de riqueza
# ---------------------------------------------------------
especies_hotspots_riqueza <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(
    id_hex %in% ids_hotspots_riqueza,
    !is.na(scientificName_det),
    scientificName_det != "") |>
  dplyr::distinct(
    id_hex,
    scientificName_det) |>
  dplyr::left_join(
    protecao_especies |>
      dplyr::select(
        scientificName_det,
        status_protecao_pontos,
        n_registros,
        n_registros_em_uc,
        percentual_registros_em_uc),
    by = "scientificName_det") |>
  dplyr::arrange(
    id_hex,
    scientificName_det
  )

# ---------------------------------------------------------
# 2. Espécies nos hotspots de CWE
# ---------------------------------------------------------
especies_hotspots_cwe <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(
    id_hex %in% ids_hotspots_cwe,
    !is.na(scientificName_det),
    scientificName_det != "") |>
  dplyr::distinct(
    id_hex,
    scientificName_det) |>
  dplyr::left_join(
    protecao_especies |>
      dplyr::select(
        scientificName_det,
        status_protecao_pontos,
        n_registros,
        n_registros_em_uc,
        percentual_registros_em_uc),
    by = "scientificName_det") |>
  dplyr::arrange(
    id_hex,
    scientificName_det
  )

# Ver resultados
print(especies_hotspots_riqueza)
print(especies_hotspots_cwe)

utils::write.csv(
  especies_hotspots_riqueza,
  file.path(
    dir_tabelas,
    "especies_hotspots_riqueza.csv"),
  row.names = FALSE
)

utils::write.csv(
  especies_hotspots_cwe,
  file.path(
    dir_tabelas,
    "especies_hotspots_CWE.csv"),
  row.names = FALSE
)
# ---------------------------------------------------------
# MAPA — HOTSPOTS DE CWE NÃO PROTEGIDOS
# ---------------------------------------------------------
# Brasil em WGS84 para plotagem
brasil_plot <- brasil_proj |> sf::st_transform(4326)

# Identificar quais hexágonos intersectam UCs
ucs_intersects <- sf::st_intersects(
  hex_50km,
  ucs_hex,
  sparse = TRUE
)

# Classificar hexágonos como protegidos/não protegidos
hex_ucs <- hex_50km |>
  dplyr::mutate(
    tem_uc = lengths(ucs_intersects) > 0
  )

# Apenas hexágonos que possuem UC
ucs_plot <- hex_ucs |>
  dplyr::filter(tem_uc) |>
  sf::st_transform(4326)

# Hotspots de CWE
hotspots_cwe_plot <- hotspots_cwe |>
  sf::st_transform(4326)

# ---------------------------------------------------------
# IDENTIFICAR HOTSPOTS DE CWE NÃO PROTEGIDOS
# ---------------------------------------------------------

hotspots_cwe_gap <- hotspots_cwe |>
  dplyr::left_join(
    hex_ucs |>
      sf::st_drop_geometry() |>
      dplyr::select(
        id_hex,
        tem_uc
      ),
    by = "id_hex"
  ) |>
  dplyr::filter(
    is.na(tem_uc) | !tem_uc
  ) |>
  sf::st_transform(4326)

# ---------------------------------------------------------
# CAMADAS PARA OS MAPAS
# ---------------------------------------------------------
estados_plot <- estados_br |>
  sf::st_transform(4326)

ucs_plot <- ucs_hex |>
  sf::st_transform(4326)

hotspots_cwe_plot <- hotspots_cwe_gap |>
  sf::st_transform(4326)

# ---------------------------------------------------------
# MAPA — HOTSPOTS DE CWE NÃO PROTEGIDOS
# ---------------------------------------------------------
mapa_gap_cwe <- ggplot2::ggplot() +
  
  # Limite do Brasil
  ggplot2::geom_sf(
    data = brasil_plot,
    fill = "gray97",
    color = "gray45",
    linewidth = 0.35
  ) +
  
  # Estados
  ggplot2::geom_sf(
    data = estados_plot,
    fill = NA,
    color = "gray65",
    linewidth = 0.25
  ) +
  
  # Unidades de Conservação
  ggplot2::geom_sf(
    data = ucs_plot,
    fill = NA,
    color = "darkgreen",
    linewidth = 0.25
  ) +
  
  # Hotspots não protegidos
  ggplot2::geom_sf(
    data = hotspots_cwe_plot,
    fill = "red3",
    color = "black",
    linewidth = 0.2
  ) +
  
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::labs(
    title = "Priority gaps in protection",
    subtitle = paste0(
      "Unprotected hotspots in the upper ",
      round((1 - percentil_hotspot) * 100, 1),
      "% of CWE; minimum richness = ",
      riqueza_minima_cwe
    )
  ) +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10
    ),
    plot.margin = ggplot2::margin(
      10, 10, 10, 10
    )
  )

print(mapa_gap_cwe)

# ---------------------------------------------------------
# SALVAR MAPA CWE
# ---------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    dir_graficos,
    "mapa_gap_CWE_hotspots_nao_protegidos.png"
  ),
  plot = mapa_gap_cwe,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)


# =========================================================
# MAPA — HOTSPOTS DE RIQUEZA NÃO PROTEGIDOS
# =========================================================

# Hotspots de riqueza
hotspots_riqueza_plot <- hotspots_riqueza |>
  sf::st_transform(4326)

# ---------------------------------------------------------
# IDENTIFICAR HOTSPOTS DE RIQUEZA NÃO PROTEGIDOS
# ---------------------------------------------------------

hotspots_riqueza_gap <- hotspots_riqueza |>
  dplyr::left_join(
    hex_ucs |>
      sf::st_drop_geometry() |>
      dplyr::select(
        id_hex,
        tem_uc
      ),
    by = "id_hex"
  ) |>
  dplyr::filter(
    is.na(tem_uc) | !tem_uc
  ) |>
  sf::st_transform(4326)

hotspots_riqueza_plot <- hotspots_riqueza_gap |>
  sf::st_transform(4326)

# ---------------------------------------------------------
# MAPA
# ---------------------------------------------------------
mapa_gap_riqueza <- ggplot2::ggplot() +
  
  # Limite do Brasil
  ggplot2::geom_sf(
    data = brasil_plot,
    fill = "gray97",
    color = "gray45",
    linewidth = 0.35
  ) +
  
  # Estados
  ggplot2::geom_sf(
    data = estados_plot,
    fill = NA,
    color = "gray65",
    linewidth = 0.25
  ) +
  
  # Unidades de Conservação
  ggplot2::geom_sf(
    data = ucs_plot,
    fill = NA,
    color = "darkgreen",
    linewidth = 0.25
  ) +
  
  # Hotspots não protegidos
  ggplot2::geom_sf(
    data = hotspots_riqueza_plot,
    fill = "red3",
    color = "black",
    linewidth = 0.2
  ) +
  
  ggplot2::coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::labs(
    title = "Priority gaps in protection",
    subtitle = paste0(
      "Unprotected hotspots in the upper ",
      round((1 - percentil_hotspot) * 100, 1),
      "% of observed richness"
    )
  ) +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10
    ),
    plot.margin = ggplot2::margin(
      10, 10, 10, 10
    )
  )

print(mapa_gap_riqueza)

# ---------------------------------------------------------
# SALVAR MAPA RIQUEZA
# ---------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    dir_graficos,
    "mapa_gap_riqueza_hotspots_nao_protegidos.png"
  ),
  plot = mapa_gap_riqueza,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggplot2::ggsave(
  filename = file.path(
    dir_graficos,
    "mapa_gap_riqueza_hotspots_nao_protegidos.png"),
  plot = mapa_gap_riqueza,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

# FIGURA CONJUNTA — HOTSPOTS DE RIQUEZA E CWE
library(patchwork)
mapas_hotspots <- mapa_gap_riqueza + mapa_gap_cwe + patchwork::plot_annotation(tag_levels = "A")

print(mapas_hotspots)

ggplot2::ggsave(
  filename = file.path(
    dir_graficos,"mapas_hotspots_riqueza_CWE.png"),
  plot = mapas_hotspots,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

# =========================================================
# GRÁFICO — STATUS DE PROTEÇÃO POR CATEGORIA DE RARIDADE
# =========================================================
# Exemplo:
# raridade_adaptada = tabela com categoria de raridade por espécie
# protecao_especies = tabela com PS / PPS / UPS por espécie

dados_raridade_protecao <- raridade_adaptada |>
  dplyr::select(
    scientificName_det,
    categoria_raridade) |>
  dplyr::left_join(
    protecao_especies |>
      dplyr::select(
        scientificName_det,
        status_protecao_pontos
      ),
    by = "scientificName_det"
  ) |>
  dplyr::filter(
    !is.na(categoria_raridade),
    !is.na(status_protecao_pontos),
    scientificName_det != "Eugenia sprengelli"
  )

# ---------------------------------------------------------
# Contagem por categoria
# ---------------------------------------------------------

resumo_raridade_protecao <- dados_raridade_protecao |>
  dplyr::count(
    categoria_raridade,
    status_protecao_pontos,
    name = "n"
  ) |>
  dplyr::group_by(categoria_raridade) |>
  dplyr::mutate(
    total_categoria = sum(n),
    proporcao = n / total_categoria
  ) |>
  dplyr::ungroup()

# ---------------------------------------------------------
# Gráfico
# ---------------------------------------------------------
grafico_raridade_protecao <- ggplot2::ggplot(
  resumo_raridade_protecao,
  ggplot2::aes(
    x = categoria_raridade,
    y = proporcao,
    fill = status_protecao_pontos)) +
  ggplot2::geom_col(
    position = "fill",
    width = 0.9) +
  
  # número de espécies dentro de cada segmento
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_fill(vjust = 0.5),
    size = 3) +
  
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)) +
  
  ggplot2::theme_minimal() +
  
  ggplot2::labs(
    title = "Protection status across multidimensional rarity categories",
    x = "Rarity category",
    y = "Proportion of species",
    fill = "Protection status") +
  
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1),
    panel.grid.major.x = ggplot2::element_blank())

print(grafico_raridade_protecao)

ggplot2::ggsave(
  filename = "grafico_raridade_status_protecao.png",
  plot = grafico_raridade_protecao,
  width = 10,
  height = 6,
  dpi = 300
)
