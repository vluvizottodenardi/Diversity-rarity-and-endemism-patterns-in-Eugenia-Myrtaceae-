# =========================================================
# 02 — ANÁLISES DESCRITIVAS: DISTRIBUIÇÃO, ESFORÇO E RIQUEZA
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)
library(writexl)

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
    ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray55", linewidth = 0.25) +
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

sf::st_write(
  hex_esforco,file.path(dir_saida, "riqueza_esforco_50km.gpkg"),
  delete_dsn = TRUE, quiet = TRUE
)

# Riqueza observada
riqueza_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(riqueza = dplyr::n_distinct(scientificName_det), .groups = "drop")
hex_riqueza <- hex_50km |> dplyr::left_join(riqueza_hex, by = "id_hex")

sf::st_write(
  hex_riqueza,file.path(dir_saida, "riqueza_hexagonos_50km.gpkg"),
  delete_dsn = TRUE, quiet = TRUE
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

# ---------------------------------------------------------
# FUNÇÃO GENÉRICA: HEXÁGONOS × CAMADA POLIGONAL
# ---------------------------------------------------------
resumo_hex_categoria <- function(hex_sf, camada_sf, coluna_categoria) {
  
  # Garantir mesma projeção antes do join espacial
  hex_sf <- sf::st_transform(hex_sf, sf::st_crs(camada_sf))
  
  # Garantir que n_registros não tenha NA (hexágonos sem registro = 0)
  hex_sf <- hex_sf |>
    dplyr::mutate(n_registros = dplyr::coalesce(n_registros, 0L))
  
  # Join espacial (um hexágono pode "cruzar" mais de uma categoria, se estiver na borda)
  join_sf <- sf::st_join(
    hex_sf,
    camada_sf[, coluna_categoria],
    join = sf::st_intersects,
    left = TRUE
  )
  
  join_df <- sf::st_drop_geometry(join_sf) |>
    dplyr::filter(!is.na(.data[[coluna_categoria]]))
  
  # Resumo por categoria (ex: por estado, por bioma, por domínio)
  resumo_por_categoria <- join_df |>
    dplyr::group_by(.data[[coluna_categoria]]) |>
    dplyr::summarise(
      n_hexagonos = dplyr::n_distinct(id_hex),
      n_hexagonos_com_registro = dplyr::n_distinct(id_hex[n_registros > 0]),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_hexagonos))
  
  # Totais gerais (hexágonos únicos, já que um mesmo hexágono pode aparecer em mais de uma categoria)
  total_hexagonos <- dplyr::n_distinct(join_df$id_hex)
  total_com_registro <- dplyr::n_distinct(join_df$id_hex[join_df$n_registros > 0])
  
  list(
    por_categoria = resumo_por_categoria,
    total_hexagonos = total_hexagonos,
    total_hexagonos_com_registro = total_com_registro
  )
}

# ---------------------------------------------------------
# Sampled effort
# ---------------------------------------------------------
# 1. HEXÁGONOS × ESTADOS
res_estados <- resumo_hex_categoria(hex_esforco, estados_br, "name_state")
cat("Total de hexágonos que cruzam algum estado:", res_estados$total_hexagonos, "\n")
cat("Desses, com pelo menos 1 registro:", res_estados$total_hexagonos_com_registro, "\n")
View(res_estados$por_categoria)

# HEXÁGONOS × REGIÕES DO BRASIL
res_regioes <- resumo_hex_categoria(hex_esforco, estados_br, "name_region")
cat("Total de hexágonos que cruzam alguma região:", res_regioes$total_hexagonos, "\n")
cat("Desses, com pelo menos 1 registro:", res_regioes$total_hexagonos_com_registro, "\n")
print(res_regioes$por_categoria)

# 2. HEXÁGONOS × BIOMAS
res_biomas <- resumo_hex_categoria(hex_esforco, biomas_proj, "bioma")
cat("Total de hexágonos que cruzam algum bioma:", res_biomas$total_hexagonos, "\n")
cat("Desses, com pelo menos 1 registro:", res_biomas$total_hexagonos_com_registro, "\n")
print(res_biomas$por_categoria)

# 3. HEXÁGONOS × DOMÍNIOS VEGETACIONAIS (habitat)
res_dominios <- resumo_hex_categoria(hex_esforco, vegetacao_proj, "habitat")
cat("Total de hexágonos que cruzam algum domínio vegetacional:", res_dominios$total_hexagonos, "\n")
cat("Desses, com pelo menos 1 registro:", res_dominios$total_hexagonos_com_registro, "\n")
print(res_dominios$por_categoria)

# ---------------------------------------------------------
# Riqueza
# ---------------------------------------------------------
riqueza_por_categoria <- function(dados_sf_pontos, camada_sf, coluna_categoria) {
  
  pontos_proj <- sf::st_transform(dados_sf_pontos, sf::st_crs(camada_sf)) |>
    dplyr::select(-dplyr::any_of(coluna_categoria))  # remove coluna conflitante, se já existir
  
  join_sf <- sf::st_join(pontos_proj, camada_sf[, coluna_categoria], join = sf::st_intersects, left = TRUE)
  join_df <- sf::st_drop_geometry(join_sf) |>
    dplyr::filter(!is.na(.data[[coluna_categoria]]), !is.na(scientificName_det), scientificName_det != "")
  
  join_df |>
    dplyr::group_by(.data[[coluna_categoria]]) |>
    dplyr::summarise(
      riqueza_total = dplyr::n_distinct(scientificName_det),
      n_registros = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(riqueza_total))
}

riqueza_estados <- riqueza_por_categoria(dados_sf, estados_br, "name_state")
riqueza_regioes <- riqueza_por_categoria(dados_sf, estados_br, "name_region")
riqueza_biomas  <- riqueza_por_categoria(dados_sf, biomas_proj, "bioma")
riqueza_dominios <- riqueza_por_categoria(dados_sf, vegetacao_proj, "habitat")

print(riqueza_regioes)
print(riqueza_biomas)
print(riqueza_dominios)

especies_biomas_lista <- especies_por_bioma |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(
    n_biomas = dplyr::n_distinct(bioma),
    biomas = paste(sort(unique(bioma)), collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::arrange(scientificName_det)

print(especies_biomas_lista)
write_xlsx(especies_biomas_lista,path = file.path(dir_tabelas, "especies_por_bioma_wide.xlsx"))

# ----------------------------------------------------------
# Top 20 hexágonos mais ricos
centroides <- hex_riqueza |> sf::st_centroid() |> sf::st_transform(4326) 
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

#Se quiser selecionar um estao
# Selecionar São Paulo
sp <- estados_br |>
  dplyr::filter(name_state == "São Paulo") |>
  sf::st_make_valid() |> sf::st_transform(5880)

# Selecionar hexágonos que intersectam São Paulo
hex_sp <- sf::st_filter(hex_riqueza, sp, .predicate = sf::st_intersects)

# Calcular centroides
centroides_sp <- sf::st_centroid(hex_sp) |> sf::st_transform(4326)

# Extrair coordenadas
coords_sp <- sf::st_coordinates(centroides_sp)

# Criar tabela com ID, riqueza e coordenadas
ids_centroides_sp <- hex_sp |>
  sf::st_drop_geometry() |>
  dplyr::select(id_hex, riqueza) |>
  dplyr::mutate(
    longitude = coords_sp[, 1],
    latitude = coords_sp[, 2])
ids_centroides_sp

# =========================================================
# EXPORTAR RESULTADOS PARA O QGIS — APENAS Algum estado
# As análises continuam sendo realizadas para todo o Brasil
# =========================================================

# Selecionar São Paulo
sp <- estados_br |>
  dplyr::filter(abbrev_state == "SP") |>
  sf::st_transform(sf::st_crs(hex_50km))

# Associar os resultados das análises aos hexágonos do Brasil
hex_resultados <- hex_50km |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::left_join(
    dados_completo |>
      dplyr::select(
        id_hex,
        n_registros,
        riqueza,
        WE,
        CWE,
        habitat,
        bioma,
        regiao
      ),
    by = "id_hex"
  )

# Selecionar apenas os hexágonos que intersectam São Paulo
hex_resultados_sp <- hex_resultados |>
  sf::st_filter(sp, .predicate = sf::st_intersects)

# Exportar para GeoPackage
sf::st_write(
  hex_resultados_sp,
  file.path(dir_vetores, "hex_resultados_SP.gpkg"),
  layer = "hex_resultados_SP",
  delete_dsn = TRUE,
  quiet = TRUE
)
