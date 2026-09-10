# =========================================================
# 5.2 — RARIDADE × CONSERVAÇÃO
# REPRESENTAÇÃO DAS ESPÉCIES RARAS EM ÁREAS PROTEGIDAS
# =========================================================
library(dplyr); library(tidyr); library(sf); library(ggplot2); library(viridis); library(patchwork); library(scales)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_raridade <- file.path("outputs", "objetos", "raridade.RData")
arquivo_conservacao <- file.path("outputs", "objetos", "conservacao.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_raridade)) stop("Execute primeiro o script de raridade 03.")
if (!file.exists(arquivo_conservacao)) stop("O arquivo conservacao.RData não foi encontrado.")

dir_saida <- file.path("outputs", "raridade_conservacao")
dir_tabelas <- file.path(dir_saida, "tabelas")
dir_mapas <- file.path(dir_saida, "mapas")
invisible(lapply(list(dir_saida, dir_tabelas, dir_mapas), dir.create, recursive = TRUE, showWarnings = FALSE))

# =========================================================
# 1. CARREGAR OBJETOS
# =========================================================
load(arquivo_base); load(arquivo_raridade); load(arquivo_conservacao)

# Critérios de proteção testados (mesma lógica de sensibilidade usada nos hotspots gerais)
criterios_protecao <- list(0.20, "totalmente_fora")
raio_buffer_m <- 25000  # raio do buffer (m) para espécies com <3 pontos únicos (não é possível calcular MCP)

# =========================================================
# 2. IDENTIFICAR ÁREAS PROTEGIDAS
# =========================================================
nomes_preferenciais <- c("ucs_union", "areas_protegidas", "protected_areas", "ucs", "UCs", "unidades_conservacao", "unidades_conservacao_union")
objetos_encontrados <- nomes_preferenciais[nomes_preferenciais %in% ls()]
if (length(objetos_encontrados) == 0) stop("Nenhum objeto de áreas protegidas foi encontrado. Defina manualmente o objeto na variável 'objeto_ucs'.")
objeto_ucs <- objetos_encontrados[1]
ucs_rede <- get(objeto_ucs)
if (!inherits(ucs_rede, "sf")) stop(paste0("O objeto '", objeto_ucs, "' não é um objeto sf."))
cat("\nCamada utilizada para áreas protegidas:", objeto_ucs, "\n")

# =========================================================
# 3. PREPARAR GEOMETRIAS E CRS
# =========================================================
crs_proj <- sf::st_crs(brasil_proj)
sf::sf_use_s2(FALSE)
ucs_rede <- sf::st_transform(ucs_rede, crs_proj) |> sf::st_make_valid()
ucs_rede <- ucs_rede[!sf::st_is_empty(ucs_rede), ]
ucs_rede <- ucs_rede[sf::st_is_valid(ucs_rede), ]
if (nrow(ucs_rede) == 0) stop("Nenhuma geometria válida permaneceu na camada de áreas protegidas.")
cat("\nGeometrias válidas de UCs:", nrow(ucs_rede), "\n")

# =========================================================
# 4. UNIR REDE DE ÁREAS PROTEGIDAS
# =========================================================
ucs_union_geom <- sf::st_union(ucs_rede)
ucs_union_proj <- sf::st_as_sf(data.frame(id = 1), geometry = ucs_union_geom) |> sf::st_make_valid()
ucs_union_map <- sf::st_transform(ucs_union_proj, 4326)

# =========================================================
# 5. PREPARAR OCORRÊNCIAS
# =========================================================
dados_sf_proj <- sf::st_transform(dados_sf, crs_proj) |> sf::st_make_valid()
dados_sf_proj <- dados_sf_proj[!sf::st_is_empty(dados_sf_proj), ]

# =========================================================
# 6. DEFINIR CATEGORIA MAIS RARA
# =========================================================
especies_mais_raras <- raridade |>
  dplyr::filter(classe_raridade == "Extremely Rare") |>
  dplyr::pull(scientificName_det)

# =========================================================
# 7. CLASSIFICAÇÃO DE PROTEÇÃO POR ESPÉCIE (pontos de ocorrência)
# ---------------------------------------------------------
# Reaproveita "protecao_especies" já calculado no script de hotspots
# (conservacao.RData), evitando recomputar a mesma coisa duas vezes
# com possível divergência entre scripts.
# =========================================================
if (exists("protecao_especies") && "status_protecao_pontos" %in% names(protecao_especies)) {
  protecao_especies <- protecao_especies |>
    dplyr::mutate(
      categoria_protecao = dplyr::case_when(
        status_protecao_pontos == "PS - Protected" ~ "PS",
        status_protecao_pontos == "PPS - Partially protected" ~ "PPS",
        status_protecao_pontos == "UPS - Unprotected" ~ "UPS"
      )
    )
  cat("\nUsando 'protecao_especies' já calculado em conservacao.RData.\n")
} else {
  dados_protecao <- dados_sf_proj |>
    dplyr::filter(!is.na(scientificName_det), scientificName_det != "")
  intersecoes_ucs <- sf::st_intersects(dados_protecao, ucs_union_proj)
  dados_protecao$em_uc <- lengths(intersecoes_ucs) > 0
  protecao_especies <- dados_protecao |>
    sf::st_drop_geometry() |>
    dplyr::group_by(scientificName_det) |>
    dplyr::summarise(
      n_registros_protecao = dplyr::n(), n_registros_uc = sum(em_uc),
      proporcao_registros_uc = n_registros_uc / n_registros_protecao,
      categoria_protecao = dplyr::case_when(
        n_registros_uc == n_registros_protecao ~ "PS",
        n_registros_uc == 0 ~ "UPS",
        TRUE ~ "PPS"
      ), .groups = "drop"
    )
  cat("\n'protecao_especies' não encontrado; recalculado neste script.\n")
}

# =========================================================
# 8. ASSOCIAR RARIDADE + PROTEÇÃO
# =========================================================
raridade_protecao <- raridade |>
  dplyr::left_join(protecao_especies |> dplyr::select(scientificName_det, categoria_protecao), by = "scientificName_det") |>
  dplyr::mutate(
    grupo_raridade = dplyr::case_when(
      is.na(classe_raridade) ~ "Não avaliada",
      classe_raridade == "Extremely Rare" ~ "Mais rara",
      TRUE ~ "Outras categorias de raridade"
    )
  )

# =========================================================
# 9. PROPORÇÃO DE PROTEÇÃO POR CATEGORIA
# =========================================================
resumo_raridade_protecao <- raridade_protecao |>
  dplyr::filter(!is.na(classe_raridade), !is.na(categoria_protecao)) |>
  dplyr::group_by(classe_raridade, categoria_protecao) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(classe_raridade) |>
  dplyr::mutate(proporcao = n / sum(n)) |>
  dplyr::ungroup()

write.csv(resumo_raridade_protecao, file.path(dir_tabelas, "protecao_por_categoria_raridade.csv"), row.names = FALSE)

# =========================================================
# 10. TABELA RARIDADE × PS/PPS/UPS
# =========================================================
tabela_raridade_protecao <- raridade_protecao |>
  dplyr::filter(!is.na(classe_raridade), !is.na(categoria_protecao)) |>
  dplyr::count(classe_raridade, categoria_protecao) |>
  tidyr::pivot_wider(names_from = categoria_protecao, values_from = n, values_fill = 0)

write.csv(tabela_raridade_protecao, file.path(dir_tabelas, "tabela_raridade_PS_PPS_UPS.csv"), row.names = FALSE)
print(tabela_raridade_protecao)

# =========================================================
# 11. TESTE DE ASSOCIAÇÃO
# =========================================================
tabela_chi <- table(raridade_protecao$classe_raridade, raridade_protecao$categoria_protecao)
teste_chi <- suppressWarnings(stats::chisq.test(tabela_chi))
print(teste_chi)
if (any(teste_chi$expected < 5)) {
  set.seed(123)
  teste_fisher <- stats::fisher.test(tabela_chi, simulate.p.value = TRUE, B = 10000)
  print(teste_fisher)
} else {
  teste_fisher <- NULL
}
teste_chi$stdres
tabela_chi

# =========================================================
# 12. FOCO NA CATEGORIA MAIS RARA
# =========================================================
mais_raras <- raridade_protecao |> dplyr::filter(grupo_raridade == "Mais rara")
resumo_mais_raras <- mais_raras |>
  dplyr::count(categoria_protecao, name = "n_especies") |>
  dplyr::mutate(proporcao = n_especies / sum(n_especies))
print(resumo_mais_raras)
write.csv(resumo_mais_raras, file.path(dir_tabelas, "protecao_especies_mais_raras.csv"), row.names = FALSE)

# =========================================================
# 13. FUNÇÃO PARA MCP (mínimo de 3 pontos únicos)
# =========================================================
calcular_mcp_protecao <- function(sp_name) {
  sp <- dados_sf_proj |> dplyr::filter(scientificName_det == sp_name)
  if (nrow(sp) < 3) return(NULL)
  coords <- sf::st_coordinates(sp)
  sp <- sp[!duplicated(data.frame(x = coords[, 1], y = coords[, 2])), ]
  if (nrow(sp) < 3) return(NULL)
  pontos_union <- sf::st_union(sp)
  mcp_geom <- sf::st_convex_hull(pontos_union)
  mcp <- sf::st_as_sf(data.frame(species = sp_name), geometry = mcp_geom) |> sf::st_make_valid()
  area_total <- as.numeric(sf::st_area(mcp))
  if (length(area_total) == 0 || is.na(area_total) || area_total <= 0) return(NULL)
  area_protegida_geom <- suppressWarnings(sf::st_intersection(mcp, ucs_union_proj))
  area_uc <- if (nrow(area_protegida_geom) == 0) 0 else sum(as.numeric(sf::st_area(area_protegida_geom)), na.rm = TRUE)
  tibble::tibble(scientificName_det = sp_name, metodo = "MCP", area_km2 = area_total / 1e6,
                 area_protegida_km2 = area_uc / 1e6, proporcao_protegida = area_uc / area_total)
}

# =========================================================
# 13b. FUNÇÃO ALTERNATIVA: BUFFER (para espécies com <3 pontos únicos)
# ---------------------------------------------------------
# O MCP exige >=3 pontos únicos para formar um polígono; espécies com
# 1-2 registros ficam de fora dessa análise (ver discussão anterior).
# Como alternativa, aplicamos um buffer fixo ao redor dos pontos —
# abordagem análoga à métrica de Area of Occupancy (AOO) da IUCN,
# que também usa uma célula/buffer fixo quando há poucos registros.
# =========================================================
calcular_buffer_protecao <- function(sp_name, raio_m = raio_buffer_m) {
  sp <- dados_sf_proj |> dplyr::filter(scientificName_det == sp_name)
  if (nrow(sp) == 0) return(NULL)
  buffer_geom <- sf::st_buffer(sp, dist = raio_m) |> sf::st_union() |> sf::st_make_valid()
  buffer_sf <- sf::st_as_sf(data.frame(species = sp_name), geometry = buffer_geom)
  area_total <- as.numeric(sf::st_area(buffer_sf))
  if (length(area_total) == 0 || is.na(area_total) || area_total <= 0) return(NULL)
  area_protegida_geom <- suppressWarnings(sf::st_intersection(buffer_sf, ucs_union_proj))
  area_uc <- if (nrow(area_protegida_geom) == 0) 0 else sum(as.numeric(sf::st_area(area_protegida_geom)), na.rm = TRUE)
  tibble::tibble(scientificName_det = sp_name, metodo = paste0("Buffer_", raio_m / 1000, "km"),
                 area_km2 = area_total / 1e6, area_protegida_km2 = area_uc / 1e6, proporcao_protegida = area_uc / area_total)
}

# =========================================================
# 14. COBERTURA DA DISTRIBUIÇÃO — ESPÉCIES MAIS RARAS (MCP + BUFFER)
# ---------------------------------------------------------
# Para cada espécie "Extremely Rare": tenta MCP primeiro; se não for
# possível (< 3 pontos únicos), usa o buffer como alternativa. Assim,
# TODAS as espécies mais raras entram na análise, não só as poucas
# com registros suficientes para formar um polígono.
# =========================================================
resultado_mcp <- dplyr::bind_rows(lapply(especies_mais_raras, calcular_mcp_protecao))
especies_sem_mcp <- setdiff(especies_mais_raras, resultado_mcp$scientificName_det)
cat("\nEspécies mais raras com MCP calculável:", nrow(resultado_mcp),
    "| sem pontos suficientes (usarão buffer):", length(especies_sem_mcp), "\n")

resultado_buffer <- dplyr::bind_rows(lapply(especies_sem_mcp, calcular_buffer_protecao))

representacao_distribuicao <- dplyr::bind_rows(resultado_mcp, resultado_buffer) |>
  dplyr::left_join(raridade |> dplyr::select(scientificName_det, n_hexagonos, n_habitats, concentracao_local, classe_raridade), by = "scientificName_det") |>
  dplyr::mutate(cobertura_classe = dplyr::case_when(
    proporcao_protegida == 0 ~ "0%", proporcao_protegida < 0.25 ~ "<25%",
    proporcao_protegida < 0.50 ~ "25–50%", proporcao_protegida < 0.75 ~ "50–75%", TRUE ~ "≥75%"
  ))

write.csv(representacao_distribuicao, file.path(dir_tabelas, "representacao_especies_mais_raras.csv"), row.names = FALSE)

# =========================================================
# 15. RESUMO DA COBERTURA DAS MAIS RARAS (por método)
# =========================================================
resumo_cobertura <- representacao_distribuicao |>
  dplyr::group_by(metodo) |>
  dplyr::summarise(n_especies = dplyr::n(), cobertura_media = mean(proporcao_protegida, na.rm = TRUE),
                   cobertura_mediana = median(proporcao_protegida, na.rm = TRUE),
                   cobertura_min = min(proporcao_protegida, na.rm = TRUE),
                   cobertura_max = max(proporcao_protegida, na.rm = TRUE), .groups = "drop")
print(resumo_cobertura)
write.csv(resumo_cobertura, file.path(dir_tabelas, "cobertura_por_metodo.csv"), row.names = FALSE)

# =========================================================
# 16. MCP + BUFFER PARA TODAS AS ESPÉCIES CLASSIFICADAS
# =========================================================
todas_especies <- raridade |> dplyr::filter(!is.na(classe_raridade)) |> dplyr::pull(scientificName_det)
resultado_mcp_todas <- dplyr::bind_rows(lapply(todas_especies, calcular_mcp_protecao))
especies_sem_mcp_todas <- setdiff(todas_especies, resultado_mcp_todas$scientificName_det)
resultado_buffer_todas <- dplyr::bind_rows(lapply(especies_sem_mcp_todas, calcular_buffer_protecao))

representacao_todas <- dplyr::bind_rows(resultado_mcp_todas, resultado_buffer_todas) |>
  dplyr::left_join(raridade |> dplyr::select(scientificName_det, classe_raridade), by = "scientificName_det")

write.csv(representacao_todas, file.path(dir_tabelas, "representacao_todas_especies.csv"), row.names = FALSE)

# =========================================================
# 17. COBERTURA POR CATEGORIA DE RARIDADE (todas as espécies)
# =========================================================
resumo_cobertura_raridade <- representacao_todas |>
  dplyr::group_by(classe_raridade) |>
  dplyr::summarise(n_especies = dplyr::n(), cobertura_media = mean(proporcao_protegida, na.rm = TRUE),
                   cobertura_mediana = median(proporcao_protegida, na.rm = TRUE),
                   cobertura_min = min(proporcao_protegida, na.rm = TRUE),
                   cobertura_max = max(proporcao_protegida, na.rm = TRUE), .groups = "drop")
print(resumo_cobertura_raridade)
write.csv(resumo_cobertura_raridade, file.path(dir_tabelas, "cobertura_por_raridade.csv"), row.names = FALSE)

# =========================================================
# 18. GRÁFICO DE COBERTURA POR CATEGORIA
# =========================================================
grafico_cobertura_raridade <- ggplot2::ggplot(representacao_todas, ggplot2::aes(x = classe_raridade, y = proporcao_protegida)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.35) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Rarity category", y = "Proportion of distribution within protected areas (MCP or buffer)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
print(grafico_cobertura_raridade)
ggplot2::ggsave(file.path(dir_saida, "cobertura_distribuicao_por_raridade.png"), grafico_cobertura_raridade, width = 11, height = 7, dpi = 300)

# =========================================================
# 19. ESPÉCIES MAIS RARAS E NÃO PROTEGIDAS
# =========================================================
especies_mais_raras_UPS <- mais_raras |>
  dplyr::filter(categoria_protecao == "UPS") |>
  dplyr::select(scientificName_det, n_hexagonos, n_habitats, concentracao_local, classe_raridade) |>
  dplyr::arrange(n_hexagonos, n_habitats, concentracao_local)
print(especies_mais_raras_UPS)
write.csv(especies_mais_raras_UPS, file.path(dir_tabelas, "especies_mais_raras_UPS.csv"), row.names = FALSE)

# =========================================================
# 20. ESPÉCIES MAIS RARAS POR HEXÁGONO + PROPORÇÃO PROTEGIDA
# =========================================================
dados_hex_raridade <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::inner_join(raridade |> dplyr::filter(classe_raridade == "Extremely Rare") |> dplyr::select(scientificName_det), by = "scientificName_det") |>
  dplyr::count(id_hex, name = "n_especies_mais_raras")

estados_br <- geobr::read_state(year = 2020, simplified = TRUE) |> sf::st_transform(4326)

hex_protecao_raridade <- hex_50km |>
  sf::st_transform(crs_proj) |>
  sf::st_make_valid() |>
  dplyr::mutate(area_hex_m2 = as.numeric(sf::st_area(geometry)))

intersecao_hex_uc <- suppressWarnings(sf::st_intersection(hex_protecao_raridade |> dplyr::select(id_hex), ucs_union_proj))

area_protegida_hex <- if (nrow(intersecao_hex_uc) > 0) {
  intersecao_hex_uc |>
    dplyr::mutate(area_uc_m2 = as.numeric(sf::st_area(geometry))) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(id_hex) |>
    dplyr::summarise(area_uc_m2 = sum(area_uc_m2, na.rm = TRUE), .groups = "drop")
} else {
  tibble::tibble(id_hex = integer(), area_uc_m2 = numeric())
}

hex_protecao_raridade <- hex_protecao_raridade |>
  dplyr::left_join(area_protegida_hex, by = "id_hex") |>
  dplyr::mutate(area_uc_m2 = dplyr::coalesce(area_uc_m2, 0), proporcao_protegida = area_uc_m2 / area_hex_m2) |>
  dplyr::left_join(dados_hex_raridade, by = "id_hex") |>
  dplyr::mutate(n_especies_mais_raras = dplyr::coalesce(n_especies_mais_raras, 0L))

# =========================================================
# 21. FUNÇÃO GENÉRICA — HEXÁGONOS PRIORITÁRIOS, DADO UM CRITÉRIO DE PROTEÇÃO
# (mesma lógica de sensibilidade usada no script de hotspots gerais: 04)
# =========================================================
brasil_plot <- sf::st_transform(brasil_proj, 4326)

gerar_mapa_raras <- function(hex_dados, criterio_protecao, titulo_sufixo) {
  if (identical(criterio_protecao, "totalmente_fora")) {
    gap <- hex_dados |> dplyr::mutate(eh_gap = proporcao_protegida == 0)
    rotulo <- "0% protected coverage"
  } else {
    gap <- hex_dados |> dplyr::mutate(eh_gap = proporcao_protegida < criterio_protecao)
    rotulo <- paste0("<", round(criterio_protecao * 100), "% protected coverage")
  }
  
  hex_prioritarios <- gap |> dplyr::filter(n_especies_mais_raras > 0, eh_gap)
  
  mapa <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = brasil_plot, fill = "gray95", color = "gray40", linewidth = 0.3) +
    ggplot2::geom_sf(data = ucs_union_map, fill = NA, color = "darkgreen", linewidth = 0.25) +
    ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray20", linewidth = 0.25) +
    ggplot2::geom_sf(data = sf::st_transform(hex_prioritarios, 4326), fill = "red3", color = "black", linewidth = 0.25) +
    ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Potential conservation gaps for the rarest species (", titulo_sufixo, ")"),
                  subtitle = paste0("Hexagons containing Extremely Rare species, with ", rotulo), x = NULL, y = NULL) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray88", linewidth = 0.3),
                   panel.grid.minor = ggplot2::element_line(color = "gray94", linewidth = 0.2),
                   panel.background = ggplot2::element_rect(fill = "white", color = NA))
  
  list(mapa = mapa, hex = hex_prioritarios, n = nrow(hex_prioritarios))
}

res_raras_02 <- gerar_mapa_raras(hex_protecao_raridade, criterios_protecao[[1]], "A")
res_raras_10 <- gerar_mapa_raras(hex_protecao_raridade, criterios_protecao[[2]], "B")

cat("\nHexágonos prioritários (espécies mais raras) — <20% protegido:", res_raras_02$n,
    "| totalmente fora (0%):", res_raras_10$n, "\n")

# =========================================================
# 22. MAPA GERAL DE DISTRIBUIÇÃO DAS ESPÉCIES MAIS RARAS (contexto, sem filtro de proteção)
# =========================================================
max_raridade <- max(hex_protecao_raridade$n_especies_mais_raras, na.rm = TRUE)
if (!is.finite(max_raridade) || max_raridade < 1) max_raridade <- 1

mapa_raridade <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = brasil_plot, fill = "gray95", color = "gray40", linewidth = 0.3) +
  ggplot2::geom_sf(data = ucs_union_map, fill = NA, color = "darkgreen", linewidth = 0.25) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray25", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_protecao_raridade, 4326), ggplot2::aes(fill = n_especies_mais_raras), color = NA) +
  viridis::scale_fill_viridis(option = "viridis", direction = 1, na.value = "transparent",
                              limits = c(1, max_raridade), breaks = seq(1, max_raridade, by = 3), name = "Rarest\nspecies") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "A. Spatial distribution of the rarest species", subtitle = "Species classified as Extremely Rare by the multidimensional rarity index", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "right", panel.grid.major = ggplot2::element_line(color = "gray85", linewidth = 0.3),
                 panel.grid.minor = ggplot2::element_line(color = "gray92", linewidth = 0.2), panel.background = ggplot2::element_rect(fill = "white", color = NA))
print(mapa_raridade)
ggplot2::ggsave(file.path(dir_mapas, "mapa_especies_mais_raras.png"), mapa_raridade, width = 10, height = 7, dpi = 300)

print(res_raras_02$mapa); print(res_raras_10$mapa)
ggplot2::ggsave(file.path(dir_mapas, "mapa_prioridades_menos20pct.png"), res_raras_02$mapa, width = 10, height = 7, dpi = 300)
ggplot2::ggsave(file.path(dir_mapas, "mapa_prioridades_totalmente_fora.png"), res_raras_10$mapa, width = 10, height = 7, dpi = 300)
sf::st_write(res_raras_02$hex,file.path(dir_mapas, "hex_prioridades_menos20pct.gpkg"),delete_dsn = TRUE, quiet = TRUE)
sf::st_write(res_raras_10$hex,file.path(dir_mapas, "hex_prioridades_totalmente_fora.gpkg"),delete_dsn = TRUE, quiet = TRUE)

# =========================================================
# 23. ESPÉCIES EM CADA HEXÁGONO PRIORITÁRIO (para os dois critérios)
# =========================================================
listar_especies_hex_prioritario <- function(hex_prioritarios) {
  dados_hex |>
    sf::st_drop_geometry() |>
    dplyr::filter(id_hex %in% hex_prioritarios$id_hex, !is.na(scientificName_det), scientificName_det != "") |>
    dplyr::distinct(id_hex, scientificName_det) |>
    dplyr::filter(scientificName_det %in% especies_mais_raras) |>
    dplyr::group_by(id_hex) |>
    dplyr::summarise(n_especies = dplyr::n(), especies = paste(sort(scientificName_det), collapse = "; "), .groups = "drop")
}

especies_por_hex_02 <- listar_especies_hex_prioritario(res_raras_02$hex)
especies_por_hex_10 <- listar_especies_hex_prioritario(res_raras_10$hex)

tabela_prioridades_02 <- res_raras_02$hex |> sf::st_drop_geometry() |> dplyr::select(id_hex, n_especies_mais_raras, proporcao_protegida) |>
  dplyr::left_join(especies_por_hex_02, by = "id_hex") |> dplyr::arrange(dplyr::desc(n_especies_mais_raras))
tabela_prioridades_10 <- res_raras_10$hex |> sf::st_drop_geometry() |> dplyr::select(id_hex, n_especies_mais_raras, proporcao_protegida) |>
  dplyr::left_join(especies_por_hex_10, by = "id_hex") |> dplyr::arrange(dplyr::desc(n_especies_mais_raras))

write.csv(tabela_prioridades_02, file.path(dir_tabelas, "hexagonos_prioritarios_menos20pct.csv"), row.names = FALSE)
write.csv(tabela_prioridades_10, file.path(dir_tabelas, "hexagonos_prioritarios_totalmente_fora.csv"), row.names = FALSE)

# =========================================================
# 24. MAPA DAS ESPÉCIES MAIS RARAS E NÃO PROTEGIDAS (por pontos, UPS)
# =========================================================
dados_hex_raras_UPS <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::inner_join(especies_mais_raras_UPS |> dplyr::select(scientificName_det), by = "scientificName_det") |>
  dplyr::count(id_hex, name = "n_raras_UPS")

hex_raras_UPS <- hex_50km |>
  dplyr::left_join(dados_hex_raras_UPS, by = "id_hex") |>
  dplyr::mutate(n_raras_UPS = dplyr::coalesce(n_raras_UPS, 0L))

mapa_raras_UPS <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = brasil_plot, fill = "gray95", color = "gray55", linewidth = 0.25) +
  ggplot2::geom_sf(data = ucs_union_map, fill = NA, color = "darkgreen", linewidth = 0.25) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray25", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_raras_UPS, 4326), ggplot2::aes(fill = n_raras_UPS), color = NA) +
  viridis::scale_fill_viridis(option = "viridis", direction = 1, na.value = "transparent", name = "Rarest UPS\nspecies") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Unprotected concentrations of the rarest *Eugenia* species", subtitle = "Hexagons containing species classified as unprotected (UPS)", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom", panel.grid = ggplot2::element_blank())
print(mapa_raras_UPS)
ggplot2::ggsave(file.path(dir_mapas, "mapa_especies_mais_raras_UPS.png"), mapa_raras_UPS, width = 10, height = 7, dpi = 300)

# =========================================================
# 25. PAINÉIS FINAIS
# =========================================================
painel_raridade_conservacao <- mapa_raridade | res_raras_02$mapa
print(painel_raridade_conservacao)
ggplot2::ggsave(file.path(dir_saida, "painel_raridade_conservacao.png"), painel_raridade_conservacao, width = 16, height = 7, dpi = 300)

painel_sensibilidade_raras <- res_raras_02$mapa + res_raras_10$mapa + patchwork::plot_annotation(tag_levels = "A", title = "Sensitivity to protection threshold — Rarest species")
print(painel_sensibilidade_raras)
ggplot2::ggsave(file.path(dir_saida, "sensibilidade_limiar_especies_raras.png"), painel_sensibilidade_raras, width = 14, height = 7, dpi = 300)

# =========================================================
# 26. SALVAR OBJETOS
# =========================================================
save(raridade_protecao, resumo_raridade_protecao, tabela_raridade_protecao, teste_chi, teste_fisher,
     mais_raras, resumo_mais_raras, representacao_distribuicao, resumo_cobertura,
     representacao_todas, resumo_cobertura_raridade, especies_mais_raras_UPS,
     dados_hex_raridade, hex_protecao_raridade, res_raras_02, res_raras_10,
     tabela_prioridades_02, tabela_prioridades_10, hex_raras_UPS,
     mapa_raridade, mapa_raras_UPS, painel_raridade_conservacao, painel_sensibilidade_raras,
     file = file.path("outputs", "objetos", "raridade_conservacao.RData"))

# =========================================================
# 27. RESUMO FINAL
# =========================================================
cat("\n=========================================================\nRARIDADE × CONSERVAÇÃO\n=========================================================\n")
cat("\nEspécies na categoria mais rara:", length(especies_mais_raras), "\n")
cat("Cobertura calculada via MCP:", nrow(resultado_mcp), "| via buffer (<3 pontos):", nrow(resultado_buffer), "\n")
cat("\nHexágonos prioritários — <20% protegido:", res_raras_02$n, "| totalmente fora (0%):", res_raras_10$n, "\n")
cat("\n--- Proteção das espécies mais raras ---\n"); print(resumo_mais_raras)
cat("\n--- Cobertura da distribuição das espécies mais raras (por método) ---\n"); print(resumo_cobertura)
cat("\n--- Teste de associação entre raridade e proteção ---\n"); print(teste_chi)
if (!is.null(teste_fisher)) print(teste_fisher)
cat("\n=========================================================\nANÁLISE FINALIZADA\n=========================================================\n")

sf::sf_use_s2(TRUE)