# =========================================================
# 03 — ENDEMISMO TAXONÔMICO: WE E CWE
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
load(arquivo_base)

dir_saida <- file.path("outputs", "endemismo")
dir.create(dir_saida, recursive = TRUE, showWarnings = FALSE)

pres_abs_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det)

range_sp_hex <- pres_abs_hex |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(n_hexagonos = dplyr::n_distinct(id_hex), .groups = "drop") |>
  dplyr::mutate(peso_endemismo = 1 / n_hexagonos)

pres_abs_peso_hex <- pres_abs_hex |> dplyr::left_join(range_sp_hex, by = "scientificName_det")

riqueza_minima <- 5  # <- ajuste esse valor conforme sua decisão metodológica, influencia o CWE

endemismo_hex <- pres_abs_peso_hex |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(
    riqueza = dplyr::n_distinct(scientificName_det),
    WE = sum(peso_endemismo),
    CWE = dplyr::if_else(riqueza >= riqueza_minima, WE / riqueza, NA_real_),
    .groups = "drop")

hex_endemismo <- hex_50km |> dplyr::left_join(endemismo_hex, by = "id_hex")

mapa_we <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray60", linewidth = 0.2) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_endemismo, 4326), ggplot2::aes(fill = WE), color = NA) +
  ggplot2::scale_fill_viridis_c(na.value = "transparent", name = "WE") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Weighted endemism", subtitle = "WE per 50 km hexagonal cell", x = NULL, y = NULL)

sf::st_write(hex_endemismo |> dplyr::select(id_hex, WE, geometry),
  file.path(dir_saida, "we_hexagonos_50km.gpkg"),delete_dsn = TRUE, quiet = TRUE
)

mapa_cwe <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray60", linewidth = 0.2) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_endemismo, 4326), ggplot2::aes(fill = CWE), color = NA) +
  ggplot2::scale_fill_viridis_c(na.value = "transparent", name = "CWE") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Corrected weighted endemism", subtitle = "CWE per 50 km hexagonal cell", x = NULL, y = NULL)

sf::st_write(hex_endemismo |> dplyr::select(id_hex, CWE, geometry),
  file.path(dir_saida, "cwe_hexagonos_50km.gpkg"), delete_dsn = TRUE, quiet = TRUE
)

print(mapa_we | mapa_cwe)
ggplot2::ggsave(file.path(dir_saida, "WE_CWE.png"), mapa_we | mapa_cwe, width = 14, height = 7, dpi = 300)

# Top 20 apenas como produto descritivo; riqueza mínima pode ser ajustada.
top_endemismo <- hex_endemismo |>
  dplyr::filter(!is.na(CWE), riqueza >= 1) |>
  dplyr::slice_max(CWE, n = 20, with_ties = FALSE)

top_endemismo_tabela <- top_endemismo |>
  sf::st_transform(4326) |>
  dplyr::mutate(
    longitude = sf::st_coordinates(sf::st_centroid(geometry))[, 1],
    latitude = sf::st_coordinates(sf::st_centroid(geometry))[, 2]) |>
  sf::st_drop_geometry() |>
  dplyr::select(dplyr::everything(), longitude, latitude)

utils::write.csv(
  top_endemismo_tabela,
  file.path(dir_tabelas, "top_20_hexagonos_endemismo.csv"),
  row.names = FALSE)

save(pres_abs_hex, range_sp_hex, endemismo_hex, hex_endemismo, mapa_we, mapa_cwe,
     file = file.path("outputs", "objetos", "endemismo_hex.RData"))

# =========================================================
# COORDENADAS DOS CENTROIDES DOS HEXÁGONOS
# =========================================================
hex_centroides <- hex_50km |>
  sf::st_centroid() |>
  sf::st_transform(4326)

coords_hex <- sf::st_coordinates(hex_centroides)

tabela_hex_coords <- hex_centroides |>
  dplyr::mutate(
    longitude = coords_hex[, 1],
    latitude = coords_hex[, 2]) |>
  sf::st_drop_geometry() |>
  dplyr::select(
    id_hex,
    longitude,
    latitude
  )

print(tabela_hex_coords)

tabela_hex_coords |> dplyr::filter(id_hex == 7391)

# =========================================================
# EXPORTAR CWE PARA O QGIS — APENAS SÃO PAULO
# As análises continuam sendo realizadas para todo o Brasil
# =========================================================

# Selecionar o estado de São Paulo
sp <- estados_br |>
  dplyr::filter(abbrev_state == "SP") |>
  sf::st_transform(sf::st_crs(hex_endemismo))

# Selecionar apenas os hexágonos que intersectam São Paulo
hex_cwe_sp <- hex_endemismo |>
  sf::st_filter(sp, .predicate = sf::st_intersects) |>
  dplyr::select(
    id_hex,
    riqueza,
    WE,
    CWE,
    geometry
  )

# Exportar para GeoPackage
sf::st_write(
  hex_cwe_sp,
  file.path(dir_saida, "cwe_hexagonos_50km_SP.gpkg"),
  layer = "cwe_hexagonos_50km_SP",
  delete_dsn = TRUE,
  quiet = TRUE
)
