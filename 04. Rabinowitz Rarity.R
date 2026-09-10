# =========================================================
# ÍNDICE CONTÍNUO DE RARIDADE + CLASSIFICAÇÃO + MAPAS + FIGURAS
# =========================================================
library(dplyr)
library(sf)
library(classInt)
library(ggplot2)
library(patchwork)
library(writexl)
library(tidyr)

# =========================================================
# PARTE 1 — ÍNDICE CONTÍNUO DE RARIDADE
# =========================================================
# 1. Presença por hexágono
pres_abs_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  distinct(id_hex, scientificName_det)

# 2. Amplitude geográfica
amplitude_geo <- pres_abs_hex |>
  group_by(scientificName_det) |>
  summarise(n_hexagonos = n_distinct(id_hex), .groups = "drop")

# 3. Distribuição dos registros entre os hexágonos
registros_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  count(scientificName_det, id_hex, name = "n_registros_hex")

# 4. Concentração espacial dos registros
concentracao_local <- registros_hex |>
  group_by(scientificName_det) |>
  mutate(
    total_registros = sum(n_registros_hex, na.rm = TRUE),
    proporcao_hex = n_registros_hex / total_registros
  ) |>
  summarise(
    concentracao_local = sum(proporcao_hex^2, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Especificidade de habitat (usando habitat dominante por hexágono, não por ponto)
amplitude_habitat <- pres_abs_hex |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::left_join(
    habitat_dominante |> dplyr::mutate(id_hex = as.character(id_hex)) |> dplyr::select(id_hex, habitat),
    by = "id_hex"
  ) |>
  dplyr::filter(!is.na(habitat)) |>
  dplyr::distinct(scientificName_det, habitat) |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(n_habitats = dplyr::n_distinct(habitat), .groups = "drop")

# 5b. Lista de vegetações onde cada espécie ocorre (por hexágono dominante, não por ponto)
habitats_por_especie <- pres_abs_hex |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::left_join(
    habitat_dominante |> dplyr::mutate(id_hex = as.character(id_hex)) |> dplyr::select(id_hex, habitat),
    by = "id_hex"
  ) |>
  dplyr::filter(!is.na(habitat)) |>
  dplyr::distinct(scientificName_det, habitat) |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(vegetacoes = paste(sort(unique(habitat)), collapse = "; "), .groups = "drop")

# 6. Lista de todas as espécies
especies <- dados |>
  filter(!is.na(scientificName_det), scientificName_det != "") |>
  distinct(scientificName_det)

# 7. Juntar as três dimensões
raridade <- especies |>
  left_join(amplitude_geo, by = "scientificName_det") |>
  left_join(concentracao_local, by = "scientificName_det") |>
  left_join(amplitude_habitat, by = "scientificName_det")

# 8. Transformar as dimensões em escores contínuos de raridade
raridade <- raridade |>
  mutate(
    log_geo = log1p(n_hexagonos),
    log_habitat = log1p(n_habitats),
    
    R_geo = 1 - (log_geo - min(log_geo, na.rm = TRUE)) /
      (max(log_geo, na.rm = TRUE) - min(log_geo, na.rm = TRUE)),
    
    R_local = (concentracao_local - min(concentracao_local, na.rm = TRUE)) /
      (max(concentracao_local, na.rm = TRUE) - min(concentracao_local, na.rm = TRUE)),
    
    R_habitat = 1 - (log_habitat - min(log_habitat, na.rm = TRUE)) /
      (max(log_habitat, na.rm = TRUE) - min(log_habitat, na.rm = TRUE))
  )

# 9. Índice multidimensional de raridade (média geométrica das 3 dimensões)
raridade <- raridade |>
  mutate(
    indice_raridade = (R_geo * R_local * R_habitat)^(1/3)
  ) |>
  arrange(desc(indice_raridade))

# Espécies não classificáveis por falta de dado de hexágono/habitat
n_nao_classificaveis <- sum(is.na(raridade$indice_raridade))
if (n_nao_classificaveis > 0) {
  message(n_nao_classificaveis, " espécie(s) ficaram sem índice de raridade calculado ",
          "(faltam dados de hexágono e/ou habitat) e não entrarão na classificação Jenks.")
}

# =========================================================
# PARTE 2 — CLASSIFICAÇÃO EM CLASSES (JENKS NATURAL BREAKS)
# =========================================================
# Remove NAs antes de calcular os cortes de Jenks (evita erro/viés nos breaks)
indice_valido <- raridade$indice_raridade[!is.na(raridade$indice_raridade)]

jenks_breaks <- classIntervals(indice_valido, n = 5, style = "jenks")

# Garante breaks únicos (Jenks pode gerar cortes duplicados se houver muitos empates)
breaks_unicos <- unique(jenks_breaks$brks)

print(jenks_breaks$brks)  # conferir os cortes

raridade <- raridade |>
  mutate(
    classe_raridade = cut(
      indice_raridade,
      breaks = breaks_unicos,
      labels = c("Common", "Uncommon", "Rare", "Very Rare", "Extremely Rare")[1:(length(breaks_unicos) - 1)],
      include.lowest = TRUE,
      ordered_result = TRUE
    )
  )

# Conferir a distribuição de espécies por classe
raridade |> count(classe_raridade)

# GVF (Goodness of Variance Fit): 1 - (variância dentro das classes / variância total)
# Quanto mais perto de 1, melhor os breaks capturam a estrutura dos dados
sdam <- sum((indice_valido - mean(indice_valido))^2)  # variância total (SDAM)
sdcm <- sum(sapply(1:(length(breaks_unicos) - 1), function(i) {
  grupo <- if (i == 1) {
    indice_valido[indice_valido >= breaks_unicos[i] & indice_valido <= breaks_unicos[i + 1]]
  } else {
    indice_valido[indice_valido > breaks_unicos[i] & indice_valido <= breaks_unicos[i + 1]]
  }
  if (length(grupo) > 1) sum((grupo - mean(grupo))^2) else 0
}))
gvf <- 1 - (sdcm / sdam)
message("GVF (Goodness of Variance Fit) = ", round(gvf, 3))

# =========================================================
# PARTE 3 — AGREGAÇÃO POR HEXÁGONO
# =========================================================
raridade_hex <- pres_abs_hex |>
  left_join(raridade |> select(scientificName_det, indice_raridade), 
            by = "scientificName_det") |>
  filter(!is.na(indice_raridade)) |>
  group_by(id_hex) |>
  summarise(
    riqueza = n_distinct(scientificName_det),
    soma_raridade = sum(indice_raridade, na.rm = TRUE),
    media_raridade = mean(indice_raridade, na.rm = TRUE),
    .groups = "drop"
  )

hex_raridade_sf <- hex_endemismo |>  # ou o objeto sf original da malha hexagonal
  left_join(raridade_hex, by = "id_hex")

# Classe de raridade dominante por hexágono (para o mapa categórico)
classe_dominante_hex <- pres_abs_hex |>
  left_join(raridade |> select(scientificName_det, classe_raridade), by = "scientificName_det") |>
  filter(!is.na(classe_raridade)) |>
  group_by(id_hex, classe_raridade) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(id_hex) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup()

hex_classe_sf <- hex_endemismo |> left_join(classe_dominante_hex, by = "id_hex")

# =========================================================
# PARTE 4 — EXPORTAR PLANILHA FINAL
# =========================================================
raridade_export <- raridade |>
  left_join(habitats_por_especie, by = "scientificName_det") |>
  select(scientificName_det, n_hexagonos, concentracao_local, n_habitats, vegetacoes,
         R_geo, R_local, R_habitat, indice_raridade, classe_raridade)

write_xlsx(
  list("Raridade_por_especie" = raridade_export),
  path = "indice_raridade_resultado.xlsx"
)

# =========================================================
# PARTE 5 — MAPAS (SOMA E MÉDIA LADO A LADO)
# =========================================================
mapa_soma <- ggplot() +
  geom_sf(data = st_transform(brasil_proj, 4326), fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  geom_sf(data = st_transform(hex_raridade_sf, 4326), aes(fill = soma_raridade), color = NA) +
  scale_fill_viridis_c(na.value = "transparent", name = "Soma") +
  coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  theme_minimal() +
  labs(title = "Soma da raridade", subtitle = "Soma do índice por célula (50 km)", x = NULL, y = NULL)

mapa_media <- ggplot() +
  geom_sf(data = st_transform(brasil_proj, 4326), fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  geom_sf(data = st_transform(hex_raridade_sf, 4326), aes(fill = media_raridade), color = NA) +
  scale_fill_viridis_c(na.value = "transparent", name = "Média") +
  coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  theme_minimal() +
  labs(title = "Média da raridade", subtitle = "Média do índice por célula (50 km)", x = NULL, y = NULL)
mapa_media

mapa_soma + mapa_media

# Top 10 hexágonos por SOMA de raridade
top10_soma <- hex_raridade_sf |>
  sf::st_drop_geometry() |>
  dplyr::slice_max(soma_raridade, n = 10, with_ties = FALSE) |>
  dplyr::select(id_hex, riqueza = riqueza.y, soma_raridade, media_raridade) |>
  dplyr::arrange(dplyr::desc(soma_raridade))
print(top10_soma)

# Top 10 hexágonos por MÉDIA de raridade
top10_media <- hex_raridade_sf |>
  sf::st_drop_geometry() |>
  dplyr::slice_max(media_raridade, n = 10, with_ties = FALSE) |>
  dplyr::select(id_hex, riqueza = riqueza.y, soma_raridade, media_raridade) |>
  dplyr::arrange(dplyr::desc(media_raridade))
print(top10_media)

# Soma
sf::st_write(
  hex_raridade_sf |> 
    sf::st_transform(4326) |> 
    dplyr::select(id_hex, soma_raridade, geometry),
  "outputs/vetores/raridade.gpkg",
  layer = "soma_raridade",
  delete_layer = TRUE
)

# Média
sf::st_write(
  hex_raridade_sf |> 
    sf::st_transform(4326) |> 
    dplyr::select(id_hex, media_raridade, geometry),
  "outputs/vetores/raridade.gpkg",
  layer = "media_raridade",
  delete_layer = TRUE
)

# Mapa categórico — classe de raridade dominante por hexágono
mapa_classe <- ggplot() +
  geom_sf(data = st_transform(brasil_proj, 4326), fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_sf(data = estados_br, fill = NA, color = "gray30", linewidth = 0.25) +
  geom_sf(data = st_transform(hex_classe_sf, 4326), aes(fill = classe_raridade), color = NA) +
  scale_fill_viridis_d(na.value = "transparent", name = "Classe\ndominante") +
  coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  theme_minimal() +
  labs(title = "Classe de raridade dominante por hexágono")

mapa_classe

# =========================================================
# PARTE 6 — FIGURAS DE APOIO PARA A DISSERTAÇÃO
# =========================================================

# 6.1 Histograma do índice com os cortes de Jenks sobrepostos
fig_histograma <- ggplot(raridade, aes(x = indice_raridade)) +
  geom_histogram(bins = 40, fill = "gray70", color = "white") +
  geom_vline(xintercept = breaks_unicos, linetype = "dashed", color = "firebrick") +
  theme_minimal() +
  labs(title = "Distribuição do índice de raridade",
       subtitle = "Linhas tracejadas = cortes de Jenks (n = 5 classes)",
       x = "Índice de raridade", y = "Número de espécies")

fig_histograma

# 6.2 Boxplot/violin das três dimensões por classe
fig_dimensoes <- raridade |>
  filter(!is.na(classe_raridade)) |>
  pivot_longer(cols = c(R_geo, R_local, R_habitat), names_to = "dimensao", values_to = "valor") |>
  ggplot(aes(x = classe_raridade, y = valor, fill = classe_raridade)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ dimensao, scales = "free_y") +
  scale_fill_viridis_d(guide = "none") +
  theme_minimal() +
  labs(title = "Contribuição de cada dimensão por classe de raridade", x = NULL, y = "Valor normalizado")

fig_dimensoes

# 6.3 Número de espécies por classe de raridade
fig_barplot_classes <- raridade |>
  filter(!is.na(classe_raridade)) |>
  count(classe_raridade) |>
  ggplot(aes(x = classe_raridade, y = n, fill = classe_raridade)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_viridis_d(guide = "none") +
  theme_minimal() +
  labs(title = "Número de espécies por classe de raridade", x = NULL, y = "Nº de espécies")

fig_barplot_classes

# 6.4 Ranking das 20 espécies com maior índice de raridade
fig_ranking <- raridade |>
  filter(!is.na(indice_raridade)) |>
  slice_max(indice_raridade, n = 20) |>
  ggplot(aes(x = reorder(scientificName_det, indice_raridade), y = indice_raridade, fill = classe_raridade)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_d(name = "Classe") +
  theme_minimal() +
  labs(title = "20 espécies com maior índice de raridade", x = NULL, y = "Índice de raridade")

fig_ranking

# 6.5 Correlação entre riqueza e soma da raridade por hexágono
fig_correlacao <- ggplot(hex_raridade_sf, aes(x = riqueza, y = soma_raridade)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "firebrick") +
  theme_minimal() +
  labs(title = "Soma da raridade vs. riqueza de espécies por hexágono",
       x = "Riqueza (nº de espécies)", y = "Soma do índice de raridade")

fig_correlacao

# =========================================================
# PARTE 7 — PAINEL COMBINADO (FIGURA ÚNICA PARA A DISSERTAÇÃO)
# =========================================================

painel_metodologico <- (fig_histograma + fig_barplot_classes) / fig_dimensoes +
  plot_annotation(tag_levels = "A")

painel_metodologico
