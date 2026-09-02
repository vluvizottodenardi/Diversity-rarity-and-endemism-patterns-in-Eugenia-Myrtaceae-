# =========================================================
# 03 — RARIDADE MULTIDIMENSIONAL — ADAPTAÇÃO DE RABINOWITZ
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
load(arquivo_base)

dir_saida <- file.path("outputs", "raridade")
dir.create(dir_saida, recursive = TRUE, showWarnings = FALSE)

pres_abs_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det)

# 1. Amplitude geográfica observada
amplitude_geo <- pres_abs_hex |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(n_hexagonos = dplyr::n_distinct(id_hex), .groups = "drop")

# 2. Amplitude de habitat observada
amplitude_habitat <- dados |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "", !is.na(habitat), habitat != "") |>
  dplyr::distinct(scientificName_det, habitat) |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(
    n_habitats = dplyr::n(),
    habitats = paste(sort(habitat), collapse = "; "),
    .groups = "drop"
  )

# 3. Frequência no banco de ocorrências
freq_registros <- dados |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "") |>
  dplyr::count(scientificName_det, name = "n_registros")

# 4. Reunir dimensões e definir cortes pela mediana
raridade_adaptada <- amplitude_geo |>
  dplyr::left_join(amplitude_habitat, by = "scientificName_det") |>
  dplyr::left_join(freq_registros, by = "scientificName_det")

corte_geo <- stats::median(raridade_adaptada$n_hexagonos, na.rm = TRUE)
corte_habitat <- stats::median(raridade_adaptada$n_habitats, na.rm = TRUE)
corte_registros <- stats::median(raridade_adaptada$n_registros, na.rm = TRUE)

raridade_adaptada <- raridade_adaptada |>
  dplyr::mutate(
    amplitude_geografica = dplyr::case_when(
      is.na(n_hexagonos) ~ NA_character_,
      n_hexagonos <= corte_geo ~ "restrita",
      TRUE ~ "ampla"),
    especificidade_habitat = dplyr::case_when(
      is.na(n_habitats) ~ NA_character_,
      n_habitats <= corte_habitat ~ "restrita",
      TRUE ~ "ampla"),
    frequencia_no_banco = dplyr::case_when(
      is.na(n_registros) ~ NA_character_,
      n_registros <= corte_registros ~ "baixa",
      TRUE ~ "alta"),
    avaliacao_completa = !is.na(amplitude_geografica) & !is.na(especificidade_habitat) & !is.na(frequencia_no_banco),
    categoria_raridade = dplyr::if_else(
      avaliacao_completa,
      paste(amplitude_geografica, especificidade_habitat, frequencia_no_banco, sep = " | "),
      "não avaliada")
  )

cortes_raridade <- data.frame(
  dimensao = c("Número de hexágonos", "Número de habitats", "Número de registros"),
  mediana = c(corte_geo, corte_habitat, corte_registros)
)

utils::write.csv(raridade_adaptada, file.path(dir_saida, "raridade_multidimensional_adaptada.csv"), row.names = FALSE)
utils::write.csv(cortes_raridade, file.path(dir_saida, "cortes_raridade_mediana.csv"), row.names = FALSE)

raridade_grafico <- raridade_adaptada |>
  dplyr::filter(categoria_raridade != "não avaliada") |>
  dplyr::mutate(
    codigo_categoria = paste0(
      ifelse(amplitude_geografica == "restrita", "G-", "G+"), " / ",
      ifelse(especificidade_habitat == "restrita", "H-", "H+"), " / ",
      ifelse(frequencia_no_banco == "baixa", "F-", "F+")
    )
  ) |>
  dplyr::count(codigo_categoria, name = "n")

grafico_raridade <- ggplot2::ggplot(raridade_grafico, ggplot2::aes(x = reorder(codigo_categoria, n), y = n)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Multidimensional rarity categories", subtitle = "G = geographic range; H = habitat breadth; F = database frequency", x = NULL, y = "Number of species")
print(grafico_raridade)

especies_mais_restritas <- raridade_adaptada |>
  dplyr::filter(amplitude_geografica == "restrita", especificidade_habitat == "restrita", frequencia_no_banco == "baixa") |>
  dplyr::arrange(n_hexagonos, n_habitats, n_registros)

# Frequência de registros por espécie
freq <- dados |> dplyr::count(scientificName_det, sort = TRUE)
grafico_freq <- ggplot2::ggplot(freq, ggplot2::aes(n)) +
  ggplot2::geom_histogram(bins = 30) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Number of records", y = "Number of species", title = "Frequency of records per species")
print(grafico_freq)

save(raridade_adaptada, cortes_raridade, especies_mais_restritas,
     file = file.path("outputs", "objetos", "raridade.RData"))


