# =========================================================
# 01 — ANÁLISES DESCRITIVAS: DISTRIBUIÇÃO, ESFORÇO E RIQUEZA
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
load(arquivo_base)

dir_saida <- file.path("outputs", "descritivos")
dir.create(dir_saida, recursive = TRUE, showWarnings = FALSE)

# Mapa geral
mapa_geral <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray30", linewidth = 0.5) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray55", linewidth = 0.25) +
  ggplot2::geom_sf(data = dados_sf, color = "red3", size = 0.5, alpha = 0.5) +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "General distribution of *Eugenia* records", x = NULL, y = NULL)
print(mapa_geral)

# Mapas automáticos por espécie
mapa_distribuicao_especie <- function(sp_name) {
  sp <- dados |> dplyr::filter(scientificName_det == sp_name)
  if (nrow(sp) == 0) return(NULL)
  mapa <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray60", linewidth = 0.2) +
    ggplot2::geom_point(data = sp, ggplot2::aes(x = longitude.gazetteer, y = latitude.gazetteer), color = "red3", size = 1, alpha = 0.7) +
    ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sp_name, subtitle = paste(nrow(sp), "records"), x = NULL, y = NULL)
  nome_sp <- gsub("[^A-Za-z0-9_]+", "_", sp_name)
  ggplot2::ggsave(file.path(dir_saida, paste0(nome_sp, ".jpg")), mapa, width = 7, height = 6, dpi = 300)
  mapa
}

species_all <- dados |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(scientificName_det) |>
  dplyr::pull(scientificName_det)

invisible(lapply(species_all, mapa_distribuicao_especie))

# Esforço amostral
esforco_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex)) |>
  dplyr::count(id_hex, name = "n_registros")
hex_esforco <- hex_50km |> dplyr::left_join(esforco_hex, by = "id_hex")

# Riqueza observada
riqueza_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(riqueza = dplyr::n_distinct(scientificName_det), .groups = "drop")
hex_riqueza <- hex_50km |> dplyr::left_join(riqueza_hex, by = "id_hex")

sf::st_write(
  hex_riqueza,
  file.path(dir_saida, "riqueza_hexagonos_50km.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

mapa_registros <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray50", linewidth = 0.25) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_esforco, 4326), ggplot2::aes(fill = n_registros), color = NA) +
  ggplot2::scale_fill_viridis_c(na.value = "transparent", name = "Records") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "A. Sampling effort", subtitle = "Number of records per 50 km hexagonal cell", x = NULL, y = NULL)

mapa_riqueza <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray50", linewidth = 0.25) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_riqueza, 4326), ggplot2::aes(fill = riqueza), color = NA) +
  ggplot2::scale_fill_viridis_c(na.value = "transparent", name = "Richness") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "B. Observed richness", subtitle = "Number of species per 50 km hexagonal cell", x = NULL, y = NULL)

mapa_hex_comparativo <- mapa_registros | mapa_riqueza
print(mapa_hex_comparativo)
ggplot2::ggsave(file.path(dir_saida, "mapa_registros_riqueza.png"), mapa_hex_comparativo, width = 14, height = 7, dpi = 300)

resumo_hexagonos <- hex_esforco |>
  dplyr::summarise(
    total_hexagonos = dplyr::n(),
    hexagonos_com_registro = sum(n_registros > 0, na.rm = TRUE),
    percentual = 100 * hexagonos_com_registro / total_hexagonos)
resumo_hexagonos

resumo_esforco <- hex_esforco |>
  dplyr::filter(n_registros > 0) |>
  dplyr::summarise(
    celulas_amostradas = dplyr::n(),
    celulas_1_5 = sum(n_registros >= 1 & n_registros <= 5),
    celulas_mais_100 = sum(n_registros > 100))
resumo_esforco

# Esforço x riqueza
comparacao_hex <- esforco_hex |> dplyr::left_join(riqueza_hex, by = "id_hex")
print(stats::cor.test(comparacao_hex$n_registros, comparacao_hex$riqueza, method = "spearman", exact = FALSE))

# Top 20 hexágonos mais ricos
centroides <- hex_riqueza |> sf::st_centroid() |> sf::st_transform(4326) #parei no 589
coords <- sf::st_coordinates(centroides)
hex_riqueza_coords <- hex_riqueza |> dplyr::mutate(longitude = coords[,1], latitude = coords[,2])
top_riqueza <- hex_riqueza_coords |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(riqueza)) |>
  dplyr::slice_max(riqueza, n = 20, with_ties = FALSE) |>
  dplyr::arrange(dplyr::desc(riqueza)) |>
  dplyr::mutate(posicao = dplyr::row_number())
utils::write.csv(top_riqueza, file.path(dir_saida, "top20_riqueza.csv"), row.names = FALSE)

save(esforco_hex, hex_esforco, riqueza_hex, hex_riqueza, comparacao_hex,
     mapa_registros, mapa_riqueza, top_riqueza,
     file = file.path("outputs", "objetos", "descritivos_hex.RData"))

especies_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(
    !is.na(id_hex),
    !is.na(scientificName_det),
    scientificName_det != ""
  ) |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::arrange(id_hex, scientificName_det)

lista_especies_hex <- especies_hex |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(
    riqueza = dplyr::n_distinct(scientificName_det),
    especies = paste(
      sort(unique(scientificName_det)),
      collapse = "; "
    ),
    .groups = "drop"
  )

# Juntar a geometria dos hexágonos
hex_especies <- hex_50km |>  dplyr::left_join(lista_especies_hex, by = "id_hex")

#se quiser salvar a planilha com os centróides dos hexágonos
centroides <- hex_riqueza |> 
  sf::st_centroid() |> 
  sf::st_transform(4326)

centroides_df <- centroides |>  sf::st_drop_geometry()

centroides_df$x <- sf::st_coordinates(centroides)[, 1]
centroides_df$y <- sf::st_coordinates(centroides)[, 2]

write.csv(
  centroides_df,
  "centroides_hexagonos.csv",
  row.names = FALSE
)

