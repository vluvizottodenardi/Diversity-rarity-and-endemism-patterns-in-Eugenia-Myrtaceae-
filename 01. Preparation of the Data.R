# =========================================================
# 00 — PREPARAÇÃO DA BASE
# Limpeza, filtro espacial, habitat e grade hexagonal de 50 km
# =========================================================
library(readxl)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(geobr)

sf::sf_use_s2(FALSE)

# ---------------------------------------------------------
# CONFIGURAÇÕES
# ---------------------------------------------------------
arquivo_dados <- "Teste.xlsx"
arquivo_vegetacao <- "C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/vege_area/vege_area.shp"
dir_objetos <- file.path("outputs", "objetos")
dir_tabelas <- file.path("outputs", "tabelas")
dir.create(dir_objetos, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabelas, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 1. IMPORTAÇÃO E LIMPEZA
# ---------------------------------------------------------
dados <- readxl::read_xlsx(arquivo_dados)

dados <- dados |>
  dplyr::mutate(
    longitude_original = longitude.gazetteer,
    latitude_original = latitude.gazetteer,
    recordedBy = trimws(tolower(recordedBy)),
    recordNumber = trimws(tolower(recordNumber)),
    longitude.gazetteer = suppressWarnings(as.numeric(gsub(",", ".", longitude.gazetteer))),
    latitude.gazetteer = suppressWarnings(as.numeric(gsub(",", ".", latitude.gazetteer)))
  ) |>
  dplyr::filter(
    !is.na(longitude.gazetteer),
    !is.na(latitude.gazetteer),
    longitude.gazetteer >= -120,
    longitude.gazetteer <= -25,
    latitude.gazetteer >= -60,
    latitude.gazetteer <= 35
  ) |>
  dplyr::group_by(
    recordedBy,
    ifelse(recordNumber == "s.n.", dplyr::row_number(), recordNumber)
  ) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

# Coordenadas suspeitas próximas a centroides do Brasil/estados.
centroides_suspeitos <- data.frame(
  localidade = c(
    "Brasil", "Acre", "Alagoas", "Amapa", "Amazonas", "Bahia","Ceara", 
    "Distrito_Federal", "Espirito_Santo", "Goias", "Maranhao", "Mato_Grosso","Mato_Grosso_do_Sul", "Minas_Gerais",
    "Para", "Paraiba", "Parana","Pernambuco", "Piaui", "Rio_de_Janeiro", "Rio_Grande_do_Norte",
    "Rio_Grande_do_Sul", "Rondonia", "Roraima", "Santa_Catarina", "Sao_Paulo", "Sergipe", "Tocantins"),
  longitude = c(
    -53.08032013, -70.44666618, -36.62320845, -51.96677136, -64.69850145, -41.73166786, -39.61790011, 
    -47.79722597, -40.66423643, -49.62332081, -45.29157809, -55.91920097, -54.84547026, -44.65770348, 
    -53.06377009, -36.82516143, -51.62693244, -38.00094632, -42.96912486, -42.65292206, -36.68393113,
    -53.31706492, -62.85654951, -61.39801946, -50.49738044, -48.7420254, -37.44607871, -48.3300012),
  latitude = c(
    -10.74257447, -9.309916736, -9.515547672, 1.437893618, -4.182240778, -12.47297852, -5.091472019, 
    -15.78111952, -19.6424025, -16.04442885, -5.079714699, -12.94959988, -20.32589839, -18.45798563, 
    -3.974215002, -7.125595572, -24.63603905, -8.326396202, -7.384484436, -22.18913072, -5.847932943,
    -29.73744154, -10.89856477, 2.084692212, -27.24916594, -22.2584562, -10.59086452, -10.15018372)
)

dados$coord_centroide <- FALSE
dados$centroide_detectado <- NA_character_

for (i in seq_len(nrow(centroides_suspeitos))) {
  teste_centroide <- abs(dados$longitude.gazetteer - centroides_suspeitos$longitude[i]) <= 0.05 &
    abs(dados$latitude.gazetteer - centroides_suspeitos$latitude[i]) <= 0.05 &
    tolower(trimws(dados$resolution.gazetteer)) != "locality"
  
  dados$coord_centroide <- dados$coord_centroide | teste_centroide
  dados$centroide_detectado[teste_centroide & is.na(dados$centroide_detectado)] <- centroides_suspeitos$localidade[i]
}

coordenadas_excluidas <- dados |> dplyr::filter(coord_centroide)

utils::write.csv(
  coordenadas_excluidas,
  file.path(dir_tabelas, "coordenadas_excluidas_centroides.csv"),
  row.names = FALSE)

dados <- dados |>
  dplyr::filter(!coord_centroide) |>
  dplyr::select(-coord_centroide, -centroide_detectado)

dados_sf <- sf::st_as_sf(
  dados,
  coords = c("longitude.gazetteer", "latitude.gazetteer"),
  crs = 4326,
  remove = FALSE)

# ---------------------------------------------------------
# 2. MAPA-BASE E FILTRO DO BRASIL
# ---------------------------------------------------------
mundo <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  dplyr::select(sovereignt, admin, continent, geometry) |>
  sf::st_make_valid()

america <- sf::st_crop(mundo, xmin = -120, xmax = -30, ymin = -60, ymax = 35)
crs_mollweide <- "+proj=moll +lon_0=-60 +datum=WGS84 +units=m +no_defs"
america_proj <- sf::st_transform(america, crs_mollweide)

estados_br <- geobr::read_state(year = 2020, simplified = TRUE) |> sf::st_make_valid()

brasil <- estados_br |>
  dplyr::summarise() |>
  sf::st_transform(sf::st_crs(dados_sf))

dados_sf <- sf::st_filter(dados_sf, brasil, .predicate = sf::st_intersects)

# Conta o número de registros por estado
estados_br <- geobr::read_state(year = 2020, simplified = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(sf::st_crs(dados_sf))

registros_estado <- sf::st_join(
  dados_sf,
  estados_br |> dplyr::select(abbrev_state, name_state),
  join = sf::st_intersects) |>
  sf::st_drop_geometry() |>
  dplyr::count(name_state, abbrev_state, name = "n_registros") |>
  dplyr::arrange(dplyr::desc(n_registros))
View(registros_estado)

#dados de frequência de espécies
distribuicao_especies <- dados_sf |>
  sf::st_drop_geometry() |>
  dplyr::count(scientificName_det, name = "n_registros")

resumo_registros <- distribuicao_especies |>
  dplyr::summarise(
    n_especies = dplyr::n(),
    minimo = min(n_registros),
    maximo = max(n_registros),
    media = mean(n_registros),
    mediana = median(n_registros),
    especies_ate_5 = sum(n_registros <= 5),
    especies_ate_10 = sum(n_registros <= 10),
    especies_ate_100 = sum(n_registros > 10),
    especies_mais_100 = sum(n_registros > 100)
  )

resumo_registros

# ---------------------------------------------------------
# 3. HABITAT/VEGETAÇÃO POR COORDENADA
# A classe usada é a coluna legenda_1 do BDiA.
# ---------------------------------------------------------
vegetacao <- sf::st_read(arquivo_vegetacao, quiet = TRUE) |> sf::st_make_valid()
dados_sf <- dados_sf |> dplyr::select(-dplyr::any_of(c("nm_pretet", "nm_pretet.x", "nm_pretet.y", "habitat")))

vegetacao_proj <- vegetacao |>
  sf::st_transform(5880) |>
  dplyr::select(habitat = legenda_1, geometry)

dados_sf_proj <- dados_sf |> sf::st_transform(5880)
dados_sf_proj <- sf::st_join(
  dados_sf_proj,
  vegetacao_proj,
  join = sf::st_intersects,
  left = TRUE
)

dados_sf <- sf::st_transform(dados_sf_proj, 4326)
dados <- sf::st_drop_geometry(dados_sf)

# ---------------------------------------------------------
# 3.1 SALVAR OCORRÊNCIAS FILTRADAS PARA USO NO QGIS
# ---------------------------------------------------------
dir_vetores <- file.path("outputs", "vetores")
dir.create(
  dir_vetores,
  recursive = TRUE,
  showWarnings = FALSE
)

# Shapefile das ocorrências finais
sf::st_write(
  dados_sf,
  dsn = file.path(
    dir_vetores,
    "ocorrencias_Eugenia_filtradas.shp"
  ),
  delete_dsn = TRUE,
  quiet = TRUE
)

# Opcional: GeoPackage, mais robusto que shapefile
sf::st_write(
  dados_sf,
  dsn = file.path(
    dir_vetores,
    "ocorrencias_Eugenia_filtradas.gpkg"),
  layer = "ocorrencias",
  delete_dsn = TRUE,
  quiet = TRUE
)

# ---------------------------------------------------------
# 4. GRADE HEXAGONAL DE 50 KM
# ---------------------------------------------------------
dados_proj <- sf::st_transform(dados_sf, 5880)
brasil_proj <- estados_br |>
  sf::st_make_valid() |>
  sf::st_transform(5880) |>
  dplyr::summarise()

hex_50km <- sf::st_make_grid(brasil_proj, cellsize = 50000, square = FALSE)
hex_50km <- sf::st_sf(
  id_hex = seq_along(hex_50km),
  geometry = hex_50km,
  crs = sf::st_crs(brasil_proj))
hex_50km <- sf::st_filter(hex_50km, brasil_proj, .predicate = sf::st_intersects)
hex_50km <- sf::st_transform(hex_50km, sf::st_crs(dados_proj))

dados_hex <- sf::st_join(
  dados_proj,
  hex_50km,
  join = sf::st_intersects,
  left = TRUE)

# ---------------------------------------------------------
# SALVAR OBJETOS-BASE
# ---------------------------------------------------------
save(
  dados, dados_sf, dados_proj, dados_hex,
  mundo, america, america_proj,
  estados_br, brasil, brasil_proj,
  hex_50km, crs_mollweide,
  file = file.path(dir_objetos, "base_analises.RData")
)

cat("Base preparada e salva em outputs/objetos/base_analises.RData\n")
sf::sf_use_s2(TRUE)

