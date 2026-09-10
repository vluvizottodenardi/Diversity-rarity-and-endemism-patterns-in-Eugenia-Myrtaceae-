# =========================================================
# SÍNTESE DE PRIORIDADE DE CONSERVAÇÃO
# Convergência entre:
# 1. Raridade multidimensional de Rabinowitz
# 2. Proteção baseada nos pontos (PS/PPS/UPS)
# 3. Gap analysis por MCP - Rodrigues et al. (2004)
# 4. Categoria de ameaça
# =========================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# Para o UpSet plot:
#install.packages("ComplexUpset")
library(ComplexUpset)
avaliacao_ameaca <- readxl::read_xlsx("avaliacao_ameaca.xlsx")

# =========================================================
# 1. CONFIGURAÇÕES
# =========================================================
dir_sintese <- file.path(
  "outputs",
  "sintese_conservacao"
)

dir.create(
  dir_sintese,
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================================================
# 2. PREPARAR RARIDADE DE RABINOWITZ
# =========================================================

# Consideramos "mais rara" a combinação:
#
# distribuição geográfica restrita
# habitat restrito
# baixa frequência no banco

raridade_sintese <- raridade_adaptada |>
  dplyr::transmute(
    scientificName_det,
    
    raridade_maxima =
      amplitude_geografica == "restrita" &
      especificidade_habitat == "restrita" &
      frequencia_no_banco == "baixa",
    
    categoria_raridade
  )

# =========================================================
# 3. PREPARAR PROTEÇÃO POR PONTOS
# =========================================================
protecao_pontos_sintese <- protecao_especies |>
  dplyr::transmute(
    scientificName_det,
    
    status_protecao_pontos,
    
    # TRUE = nenhum registro conhecido dentro de UC
    sem_protecao_pontos =
      status_protecao_pontos ==
      "UPS - Unprotected"
  )

# =========================================================
# 4. PREPARAR GAP ANALYSIS DE RODRIGUES / MCP
# =========================================================

gap_mcp_sintese <- gap_rodrigues_poligonos |>
  dplyr::transmute(
    scientificName_det,
    
    categoria_gap,
    
    # TRUE quando a meta de representação
    # espacial NÃO foi atingida.
    #
    # Inclui:
    # - Lacuna parcial
    # - Não protegida (gap)
    #
    # Espécies "Não avaliadas" não entram
    # como gap automaticamente.
    meta_mcp_nao_atingida =
      categoria_gap %in% c(
        "Lacuna parcial",
        "Não protegida (gap)"
      ),
    
    mcp_nao_avaliado =
      grepl(
        "^Não avaliada",
        categoria_gap
      )
  )

# =========================================================
# 5. PREPARAR CATEGORIA DE AMEAÇA
# =========================================================

ameaca_sintese <- avaliacao_ameaca |>
  dplyr::transmute(
    scientificName_det,
    
    categoria_ameaca,
    
    # Espécies ameaçadas:
    # CR + EN + VU
    ameacada =
      categoria_ameaca %in%
      c("CR", "EN", "VU"),
    
    nao_avaliada_ameaca =
      categoria_ameaca %in%
      c("NE", NA_character_)
  )


# =========================================================
# 6. JUNTAR TODAS AS ANÁLISES
# =========================================================

sintese_conservacao <- raridade_sintese |>
  
  dplyr::full_join(
    protecao_pontos_sintese,
    by = "scientificName_det"
  ) |>
  
  dplyr::full_join(
    gap_mcp_sintese,
    by = "scientificName_det"
  ) |>
  
  dplyr::full_join(
    ameaca_sintese,
    by = "scientificName_det"
  )


# =========================================================
# 7. GARANTIR TRUE/FALSE
# =========================================================

sintese_conservacao <- sintese_conservacao |>
  dplyr::mutate(
    
    raridade_maxima =
      dplyr::coalesce(
        raridade_maxima,
        FALSE
      ),
    
    sem_protecao_pontos =
      dplyr::coalesce(
        sem_protecao_pontos,
        FALSE
      ),
    
    meta_mcp_nao_atingida =
      dplyr::coalesce(
        meta_mcp_nao_atingida,
        FALSE
      ),
    
    ameacada =
      dplyr::coalesce(
        ameacada,
        FALSE
      )
  )


# =========================================================
# 8. NÚMERO DE CRITÉRIOS DE PREOCUPAÇÃO
# =========================================================
#
# IMPORTANTE:
#
# Isto NÃO é um novo índice de raridade.
#
# É apenas o número de sinais independentes
# de preocupação/conservação que convergem
# para cada espécie.
# =========================================================

sintese_conservacao <- sintese_conservacao |>
  dplyr::mutate(
    
    n_criterios_conservacao =
      as.integer(raridade_maxima) +
      as.integer(sem_protecao_pontos) +
      as.integer(meta_mcp_nao_atingida) +
      as.integer(ameacada)
  )


# =========================================================
# 9. ORDENAR ESPÉCIES
# =========================================================

sintese_conservacao <- sintese_conservacao |>
  dplyr::arrange(
    dplyr::desc(
      n_criterios_conservacao
    ),
    scientificName_det
  )


# Ver as espécies com maior convergência
print(
  sintese_conservacao |>
    dplyr::select(
      scientificName_det,
      categoria_raridade,
      status_protecao_pontos,
      categoria_gap,
      categoria_ameaca,
      raridade_maxima,
      sem_protecao_pontos,
      meta_mcp_nao_atingida,
      ameacada,
      n_criterios_conservacao
    ),
  n = 100
)


# =========================================================
# 10. ESPÉCIES EM QUE OS 4 CRITÉRIOS CONVERGEM
# =========================================================

prioridade_convergente <- sintese_conservacao |>
  dplyr::filter(
    raridade_maxima,
    sem_protecao_pontos,
    meta_mcp_nao_atingida,
    ameacada
  ) |>
  dplyr::arrange(
    scientificName_det
  )

print(prioridade_convergente)


# =========================================================
# 11. ESPÉCIES COM PELO MENOS 3 CRITÉRIOS
# =========================================================

alta_convergencia <- sintese_conservacao |>
  dplyr::filter(
    n_criterios_conservacao >= 3
  ) |>
  dplyr::arrange(
    dplyr::desc(
      n_criterios_conservacao
    ),
    scientificName_det
  )

print(alta_convergencia)


# =========================================================
# 12. RESUMO DO NÚMERO DE CRITÉRIOS
# =========================================================

resumo_criterios <- sintese_conservacao |>
  dplyr::count(
    n_criterios_conservacao,
    name = "n_especies"
  ) |>
  dplyr::mutate(
    percentual =
      100 *
      n_especies /
      sum(n_especies)
  )

print(resumo_criterios)


# =========================================================
# 13. UPSET PLOT
# =========================================================
#
# Mostra quantas espécies pertencem simultaneamente
# aos diferentes conjuntos.
# =========================================================

dados_upset <- sintese_conservacao |>
  dplyr::select(
    scientificName_det,
    raridade_maxima,
    sem_protecao_pontos,
    meta_mcp_nao_atingida,
    ameacada
  ) |>
  dplyr::rename(
    `Highest Rabinowitz rarity` =
      raridade_maxima,
    
    `Unprotected occurrences (UPS)` =
      sem_protecao_pontos,
    
    `MCP target not achieved` =
      meta_mcp_nao_atingida,
    
    `Threatened (CR/EN/VU)` =
      ameacada
  )


grafico_upset <- ComplexUpset::upset(
  dados_upset,
  
  intersect = c(
    "Highest Rabinowitz rarity",
    "Unprotected occurrences (UPS)",
    "MCP target not achieved",
    "Threatened (CR/EN/VU)"
  ),
  
  min_size = 1,
  
  width_ratio = 0.20,
  
  base_annotations = list(
    "Intersection size" =
      ComplexUpset::intersection_size(
        text = list(
          size = 3.5
        )
      )
  )
) +
  
  ggplot2::labs(
    title =
      "Convergence of conservation concern criteria in Eugenia"
  )

print(grafico_upset)


# =========================================================
# 14. SALVAR UPSET
# =========================================================

ggplot2::ggsave(
  filename = file.path(
    dir_sintese,
    "UpSet_convergencia_conservacao_Eugenia.png"
  ),
  plot = grafico_upset,
  width = 11,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)


# =========================================================
# 15. GRÁFICO — QUANTOS CRITÉRIOS CADA ESPÉCIE ACUMULA
# =========================================================

grafico_n_criterios <- ggplot2::ggplot(
  resumo_criterios,
  ggplot2::aes(
    x = factor(
      n_criterios_conservacao
    ),
    y = n_especies
  )
) +
  
  ggplot2::geom_col(
    width = 0.75
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = n_especies
    ),
    vjust = -0.3,
    size = 3.5
  ) +
  
  ggplot2::theme_minimal() +
  
  ggplot2::labs(
    title =
      "Convergence of conservation concern criteria",
    x =
      "Number of conservation concern criteria",
    y =
      "Number of species"
  )

print(grafico_n_criterios)


ggplot2::ggsave(
  filename = file.path(
    dir_sintese,
    "numero_criterios_conservacao.png"
  ),
  plot = grafico_n_criterios,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# =========================================================
# 16. SALVAR TABELAS
# =========================================================

utils::write.csv(
  sintese_conservacao,
  file.path(
    dir_sintese,
    "sintese_conservacao_todas_especies.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  prioridade_convergente,
  file.path(
    dir_sintese,
    "especies_4_criterios_convergentes.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  alta_convergencia,
  file.path(
    dir_sintese,
    "especies_3_ou_4_criterios.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  resumo_criterios,
  file.path(
    dir_sintese,
    "resumo_numero_criterios.csv"
  ),
  row.names = FALSE
)


# =========================================================
# ESPÉCIES COM PERFIL MAIS EXTREMO DE RARIDADE OBSERVADA
# Critério:
# 1 hexágono + 1 habitat + 1 registro
# =========================================================

library(dplyr)
library(sf)

# ---------------------------------------------------------
# 1. IDENTIFICAR AS ESPÉCIES MAIS EXTREMAMENTE RARAS
# ---------------------------------------------------------

especies_raridade_extrema <- raridade_adaptada |>
  dplyr::filter(
    n_hexagonos == 1,
    n_habitats == 1,
    n_registros == 1
  ) |>
  dplyr::arrange(scientificName_det)

print(especies_raridade_extrema)

cat(
  "\nNúmero de espécies com 1 hexágono + 1 habitat + 1 registro:",
  nrow(especies_raridade_extrema),
  "\n"
)


# =========================================================
# 2. JUNTAR INFORMAÇÕES DE CONSERVAÇÃO
# =========================================================

# Junta o status PS / PPS / UPS, caso o objeto exista.
if (exists("protecao_especies")) {
  
  especies_raridade_extrema <- especies_raridade_extrema |>
    dplyr::left_join(
      protecao_especies |>
        dplyr::select(
          scientificName_det,
          status_protecao_pontos,
          n_registros_em_uc,
          percentual_registros_em_uc
        ),
      by = "scientificName_det"
    )
}


# Junta o resultado de Rodrigues/MCP, caso exista.
# Muitas dessas espécies provavelmente serão "Não avaliadas",
# pois possuem apenas uma localidade.
if (exists("gap_rodrigues_poligonos")) {
  
  especies_raridade_extrema <- especies_raridade_extrema |>
    dplyr::left_join(
      gap_rodrigues_poligonos |>
        dplyr::select(
          scientificName_det,
          categoria_gap,
          area_poligono_km2,
          percentual_protegido,
          percentual_meta
        ),
      by = "scientificName_det"
    )
}


# ---------------------------------------------------------
# OPCIONAL — CATEGORIA DE AMEAÇA
# ---------------------------------------------------------
#
# Se sua tabela tiver outro nome, troque "avaliacao_ameaca".
# Também ajuste "categoria_ameaca" se necessário.

if (exists("avaliacao_ameaca")) {
  
  especies_raridade_extrema <- especies_raridade_extrema |>
    dplyr::left_join(
      avaliacao_ameaca |>
        dplyr::select(
          scientificName_det,
          categoria_ameaca
        ),
      by = "scientificName_det"
    )
}


# =========================================================
# 3. SALVAR TABELA
# =========================================================

dir_raridade_extrema <- file.path(
  "outputs",
  "raridade_extrema"
)

dir.create(
  dir_raridade_extrema,
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.csv(
  especies_raridade_extrema,
  file.path(
    dir_raridade_extrema,
    "especies_raridade_extrema.csv"
  ),
  row.names = FALSE
)


# =========================================================
# 4. EXTRAIR OS PONTOS DESSAS ESPÉCIES
# =========================================================

nomes_raridade_extrema <-
  especies_raridade_extrema$scientificName_det

pontos_raridade_extrema <- dados_sf |>
  dplyr::filter(
    scientificName_det %in% nomes_raridade_extrema
  )


# =========================================================
# 5. JUNTAR INFORMAÇÕES NA CAMADA ESPACIAL
# =========================================================

pontos_raridade_extrema <- pontos_raridade_extrema |>
  dplyr::left_join(
    especies_raridade_extrema,
    by = "scientificName_det"
  )


# =========================================================
# 6. EXPORTAR PARA QGIS
# =========================================================

dir_vetores <- file.path(
  dir_raridade_extrema,
  "vetores"
)

dir.create(
  dir_vetores,
  recursive = TRUE,
  showWarnings = FALSE
)

# GeoPackage — recomendado
sf::st_write(
  pontos_raridade_extrema,
  dsn = file.path(
    dir_vetores,
    "especies_raridade_extrema.gpkg"
  ),
  layer = "ocorrencias_raridade_extrema",
  delete_layer = TRUE,
  quiet = TRUE
)

# Shapefile — opcional
sf::st_write(
  pontos_raridade_extrema,
  dsn = file.path(
    dir_vetores,
    "especies_raridade_extrema.shp"
  ),
  delete_dsn = TRUE,
  quiet = TRUE
)


# =========================================================
# 7. RESUMO
# =========================================================

cat(
  "\nEspécies identificadas:",
  nrow(especies_raridade_extrema),
  "\nPontos exportados:",
  nrow(pontos_raridade_extrema),
  "\n"
)

# =========================================================
# 8. MAPA DAS ESPÉCIES COM RARIDADE EXTREMA
# =========================================================

library(ggplot2)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)

# ---------------------------------------------------------
# 8.1. Limites do Brasil
# ---------------------------------------------------------

brasil <- rnaturalearth::ne_countries(
  country = "Brazil",
  scale = "medium",
  returnclass = "sf"
)

# Garantir mesmo CRS
brasil <- sf::st_transform(
  brasil,
  sf::st_crs(pontos_raridade_extrema)
)


# ---------------------------------------------------------
# 8.2. Preparar nomes para os rótulos
# ---------------------------------------------------------

pontos_mapa <- pontos_raridade_extrema |>
  dplyr::mutate(
    label = scientificName_det
  )


# ---------------------------------------------------------
# 8.3. Mapa
# ---------------------------------------------------------
estados <- rnaturalearth::ne_states(
  country = "Brazil",
  returnclass = "sf"
) |>
  sf::st_transform(sf::st_crs(pontos_raridade_extrema))

mapa_raridade_extrema <- ggplot() +
  geom_sf(data = brasil, fill = "grey95", color = "grey30", linewidth = 0.3) +
  geom_sf(data = estados, fill = NA, color = "grey60", linewidth = 0.2) +
  geom_sf(data = pontos_mapa, aes(color = scientificName_det), size = 3, alpha = 0.9) +
  geom_text_repel(
    data = sf::st_drop_geometry(pontos_mapa) |>
      dplyr::mutate(
        x = sf::st_coordinates(pontos_mapa)[, 1],
        y = sf::st_coordinates(pontos_mapa)[, 2]
      ),
    aes(x = x, y = y, label = label, color = scientificName_det),
    size = 3,
    fontface = "bold.italic",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_discrete(name = "Species") +
  coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  labs(
    title = "Extremely rare species",
    subtitle = "Species represented by a single record in one hexagon and one habitat",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10)
  )

print(mapa_raridade_extrema)


# ---------------------------------------------------------
# 8.5. Salvar
# ---------------------------------------------------------

ggsave(
  filename = file.path(
    dir_raridade_extrema,
    "mapa_especies_raridade_extrema.png"
  ),
  plot = mapa_raridade_extrema,
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)

