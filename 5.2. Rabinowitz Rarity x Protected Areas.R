# =========================================================
# 04 — RARIDADE × CONSERVAÇÃO
# REPRESENTAÇÃO DAS ESPÉCIES RARAS EM ÁREAS PROTEGIDAS
# =========================================================
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_raridade <- file.path("outputs", "objetos", "raridade.RData")
arquivo_conservacao <- file.path("outputs", "objetos", "conservacao.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_raridade)) stop("Execute primeiro o script de raridade 03.")
if (!file.exists(arquivo_conservacao)) stop("O arquivo conservacao.RData não foi encontrado.")

dir_saida <- file.path("outputs", "raridade_conservacao")
dir_tabelas <- file.path(dir_saida, "tabelas")
dir_mapas <- file.path(dir_saida, "mapas")
dir.create(dir_saida, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabelas, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_mapas, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1. CARREGAR OBJETOS
# =========================================================
load(arquivo_base)
load(arquivo_raridade)
load(arquivo_conservacao)

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
ucs_rede <- sf::st_transform(ucs_rede, crs_proj)
ucs_rede <- sf::st_make_valid(ucs_rede)
ucs_rede <- ucs_rede[!sf::st_is_empty(ucs_rede), ]
ucs_rede <- ucs_rede[sf::st_is_valid(ucs_rede), ]
if (nrow(ucs_rede) == 0) stop("Nenhuma geometria válida permaneceu na camada de áreas protegidas.")
cat("\nGeometrias válidas de UCs:", nrow(ucs_rede), "\n")

# =========================================================
# 4. UNIR REDE DE ÁREAS PROTEGIDAS
# =========================================================
ucs_union_geom <- sf::st_union(ucs_rede)
ucs_union_proj <- sf::st_as_sf(data.frame(id = 1), geometry = ucs_union_geom)
ucs_union_proj <- sf::st_make_valid(ucs_union_proj)
ucs_union_map <- sf::st_transform(ucs_union_proj, 4326)

# =========================================================
# 5. PREPARAR OCORRÊNCIAS
# =========================================================
dados_sf_proj <- sf::st_transform(dados_sf, crs_proj)
dados_sf_proj <- sf::st_make_valid(dados_sf_proj)
dados_sf_proj <- dados_sf_proj[!sf::st_is_empty(dados_sf_proj), ]

# =========================================================
# 6. DEFINIR CATEGORIA MAIS RARA
# =========================================================
especies_mais_raras <- raridade_adaptada |>
  dplyr::filter(amplitude_geografica == "restrita", especificidade_habitat == "restrita", frequencia_no_banco == "baixa") |>
  dplyr::pull(scientificName_det)
cat("\nNúmero de espécies na categoria mais rara:", length(especies_mais_raras), "\n")

# =========================================================
# 7. CLASSIFICAÇÃO DE PROTEÇÃO POR ESPÉCIE
# =========================================================
dados_protecao <- dados_sf_proj |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "")
intersecoes_ucs <- sf::st_intersects(dados_protecao, ucs_union_proj)
dados_protecao$em_uc <- lengths(intersecoes_ucs) > 0
protecao_especies <- dados_protecao |>
  sf::st_drop_geometry() |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(n_registros_protecao = dplyr::n(), n_registros_uc = sum(em_uc), proporcao_registros_uc = n_registros_uc / n_registros_protecao, categoria_protecao = dplyr::case_when(n_registros_uc == n_registros_protecao ~ "PS", n_registros_uc == 0 ~ "UPS", TRUE ~ "PPS"), .groups = "drop")

# =========================================================
# 8. ASSOCIAR RARIDADE + PROTEÇÃO
# =========================================================
raridade_protecao <- raridade_adaptada |>
  dplyr::left_join(protecao_especies, by = "scientificName_det") |>
  dplyr::mutate(grupo_raridade = dplyr::case_when(categoria_raridade == "não avaliada" ~ "Não avaliada", amplitude_geografica == "restrita" & especificidade_habitat == "restrita" & frequencia_no_banco == "baixa" ~ "Mais rara", TRUE ~ "Outras categorias de raridade"))

# =========================================================
# 9. PROPORÇÃO DE PROTEÇÃO POR CATEGORIA
# =========================================================
resumo_raridade_protecao <- raridade_protecao |>
  dplyr::filter(categoria_raridade != "não avaliada", !is.na(categoria_protecao)) |>
  dplyr::group_by(categoria_raridade, categoria_protecao) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(categoria_raridade) |>
  dplyr::mutate(proporcao = n / sum(n)) |>
  dplyr::ungroup()
write.csv(resumo_raridade_protecao, file.path(dir_tabelas, "protecao_por_categoria_raridade.csv"), row.names = FALSE)

# =========================================================
# 10. TABELA RARIDADE × PS/PPS/UPS
# =========================================================
tabela_raridade_protecao <- raridade_protecao |>
  dplyr::filter(categoria_raridade != "não avaliada", !is.na(categoria_protecao)) |>
  dplyr::count(categoria_raridade, categoria_protecao) |>
  tidyr::pivot_wider(names_from = categoria_protecao, values_from = n, values_fill = 0)
write.csv(tabela_raridade_protecao, file.path(dir_tabelas, "tabela_raridade_PS_PPS_UPS.csv"), row.names = FALSE)
print(tabela_raridade_protecao)

# =========================================================
# 11. TESTE DE ASSOCIAÇÃO
# =========================================================
tabela_chi <- table(raridade_protecao$categoria_raridade, raridade_protecao$categoria_protecao)
teste_chi <- suppressWarnings(stats::chisq.test(tabela_chi))
print(teste_chi)
if (any(teste_chi$expected < 5)) { set.seed(123); teste_fisher <- stats::fisher.test(tabela_chi, simulate.p.value = TRUE, B = 10000); print(teste_fisher) } else { teste_fisher <- NULL }

# =========================================================
# 12. FOCO NA CATEGORIA MAIS RARA
# =========================================================
mais_raras <- raridade_protecao |>
  dplyr::filter(grupo_raridade == "Mais rara")
resumo_mais_raras <- mais_raras |>
  dplyr::count(categoria_protecao, name = "n_especies") |>
  dplyr::mutate(proporcao = n_especies / sum(n_especies))
print(resumo_mais_raras)
write.csv(resumo_mais_raras, file.path(dir_tabelas, "protecao_especies_mais_raras.csv"), row.names = FALSE)

# =========================================================
# 13. FUNÇÃO PARA MCP
# =========================================================
calcular_mcp_protecao <- function(sp_name) {
  sp <- dados_sf_proj |>
    dplyr::filter(scientificName_det == sp_name)
  if (nrow(sp) < 3) return(NULL)
  coords <- sf::st_coordinates(sp)
  sp <- sp[!duplicated(data.frame(x = coords[, 1], y = coords[, 2])), ]
  if (nrow(sp) < 3) return(NULL)
  pontos_union <- sf::st_union(sp)
  mcp_geom <- sf::st_convex_hull(pontos_union)
  mcp <- sf::st_as_sf(data.frame(species = sp_name), geometry = mcp_geom)
  mcp <- sf::st_make_valid(mcp)
  area_total <- as.numeric(sf::st_area(mcp))
  if (length(area_total) == 0 || is.na(area_total) || area_total <= 0) return(NULL)
  area_protegida_geom <- suppressWarnings(sf::st_intersection(mcp, ucs_union_proj))
  if (nrow(area_protegida_geom) == 0) { area_uc <- 0 } else { area_uc <- sum(as.numeric(sf::st_area(area_protegida_geom)), na.rm = TRUE) }
  tibble::tibble(scientificName_det = sp_name, area_mcp = area_total, area_protegida = area_uc, proporcao_mcp_protegida = area_uc / area_total)
}

# =========================================================
# 14. MCP DAS ESPÉCIES MAIS RARAS
# =========================================================
representacao_distribuicao <- dplyr::bind_rows(lapply(especies_mais_raras, calcular_mcp_protecao)) |>
  dplyr::left_join(raridade_adaptada |>
                     dplyr::select(scientificName_det, n_hexagonos, n_habitats, n_registros, categoria_raridade), by = "scientificName_det") |>
  dplyr::mutate(cobertura_mcp_classe = dplyr::case_when(proporcao_mcp_protegida == 0 ~ "0%", proporcao_mcp_protegida < 0.25 ~ "<25%", proporcao_mcp_protegida < 0.50 ~ "25–50%", proporcao_mcp_protegida < 0.75 ~ "50–75%", TRUE ~ "≥75%"))
write.csv(representacao_distribuicao, file.path(dir_tabelas, "representacao_MCP_especies_mais_raras.csv"), row.names = FALSE)

# =========================================================
# 15. RESUMO DA COBERTURA DAS MAIS RARAS
# =========================================================
resumo_MCP <- representacao_distribuicao |>
  dplyr::summarise(n_especies = dplyr::n(), cobertura_media = mean(proporcao_mcp_protegida, na.rm = TRUE), cobertura_mediana = median(proporcao_mcp_protegida, na.rm = TRUE), cobertura_min = min(proporcao_mcp_protegida, na.rm = TRUE), cobertura_max = max(proporcao_mcp_protegida, na.rm = TRUE))
print(resumo_MCP)

# =========================================================
# 16. MCP PARA TODAS AS CATEGORIAS
# =========================================================
todas_especies_mcp <- raridade_adaptada |>
  dplyr::filter(categoria_raridade != "não avaliada") |>
  dplyr::pull(scientificName_det)
representacao_todas <- dplyr::bind_rows(lapply(todas_especies_mcp, calcular_mcp_protecao)) |>
  dplyr::left_join(raridade_adaptada |>
                     dplyr::select(scientificName_det, categoria_raridade), by = "scientificName_det")
write.csv(representacao_todas, file.path(dir_tabelas, "representacao_MCP_todas_especies.csv"), row.names = FALSE)

# =========================================================
# 17. COBERTURA POR CATEGORIA DE RARIDADE
# =========================================================
resumo_cobertura_raridade <- representacao_todas |>
  dplyr::group_by(categoria_raridade) |>
  dplyr::summarise(n_especies = dplyr::n(), cobertura_media = mean(proporcao_mcp_protegida, na.rm = TRUE), cobertura_mediana = median(proporcao_mcp_protegida, na.rm = TRUE), cobertura_min = min(proporcao_mcp_protegida, na.rm = TRUE), cobertura_max = max(proporcao_mcp_protegida, na.rm = TRUE), .groups = "drop")
print(resumo_cobertura_raridade)
write.csv(resumo_cobertura_raridade, file.path(dir_tabelas, "cobertura_MCP_por_raridade.csv"), row.names = FALSE)

# =========================================================
# 18. GRÁFICO DE COBERTURA POR CATEGORIA
# =========================================================
grafico_cobertura_raridade <- ggplot2::ggplot(representacao_todas, ggplot2::aes(x = categoria_raridade, y = proporcao_mcp_protegida)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.35) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Rarity category", y = "Proportion of observed distribution within protected areas") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
print(grafico_cobertura_raridade)
ggplot2::ggsave(file.path(dir_saida, "cobertura_distribuicao_por_raridade.png"), grafico_cobertura_raridade, width = 11, height = 7, dpi = 300)

# =========================================================
# 19. ESPÉCIES MAIS RARAS E NÃO PROTEGIDAS
# =========================================================
especies_mais_raras_UPS <- mais_raras |>
  dplyr::filter(categoria_protecao == "UPS") |>
  dplyr::select(scientificName_det, n_hexagonos, n_habitats, n_registros, categoria_raridade) |>
  dplyr::arrange(n_hexagonos, n_habitats, n_registros)
print(especies_mais_raras_UPS)
write.csv(especies_mais_raras_UPS, file.path(dir_tabelas, "especies_mais_raras_UPS.csv"), row.names = FALSE)

# =========================================================
# 20. ESPÉCIES MAIS RARAS POR HEXÁGONO
# =========================================================
dados_hex_raridade <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::inner_join(raridade_adaptada |>
                      dplyr::filter(amplitude_geografica == "restrita", especificidade_habitat == "restrita", frequencia_no_banco == "baixa") |>
                      dplyr::select(scientificName_det), by = "scientificName_det") |>
  dplyr::count(id_hex, name = "n_especies_mais_raras")

# =========================================================
# 21. MAPA DA DISTRIBUIÇÃO DAS ESPÉCIES MAIS RARAS
# =========================================================
hex_raridade <- hex_50km |>
  dplyr::left_join(dados_hex_raridade, by = "id_hex") |>
  dplyr::mutate(n_especies_mais_raras = dplyr::coalesce(n_especies_mais_raras, 0L))
max_raridade <- max(hex_raridade$n_especies_mais_raras, na.rm = TRUE)
if (!is.finite(max_raridade) || max_raridade < 1) max_raridade <- 1

estados_br <- geobr::read_state(year = 2020, simplified = TRUE) |>
  sf::st_transform(4326)

mapa_raridade <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray40", linewidth = 0.3) +
  ggplot2::geom_sf(data = ucs_union_map, fill = NA, color = "darkgreen", linewidth = 0.25) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray25", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_raridade, 4326), ggplot2::aes(fill = n_especies_mais_raras), color = NA) +
  viridis::scale_fill_viridis(option = "viridis", direction = 1, na.value = "transparent", limits = c(1, max_raridade), breaks = seq(1, max_raridade, by = 3), name = "Rarest\nspecies") +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "A. Spatial distribution of the rarest species", subtitle = "Species simultaneously restricted in geographic range, habitat breadth and occurrence frequency", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "right", panel.grid.major = ggplot2::element_line(color = "gray85", linewidth = 0.3), panel.grid.minor = ggplot2::element_line(color = "gray92", linewidth = 0.2), panel.background = ggplot2::element_rect(fill = "white", color = NA))
print(mapa_raridade)
ggplot2::ggsave(file.path(dir_mapas, "mapa_especies_mais_raras.png"), mapa_raridade, width = 10, height = 7, dpi = 300)

# =========================================================
# 22. PROPORÇÃO DE CADA HEXÁGONO DENTRO DE UCS
# =========================================================
hex_protecao_raridade <- hex_50km |>
  sf::st_transform(crs_proj) |>
  sf::st_make_valid()

area_total_hex <- hex_protecao_raridade |>
  dplyr::mutate(area_hex_m2 = as.numeric(sf::st_area(geometry))) |>
  sf::st_drop_geometry() |>
  dplyr::select(id_hex, area_hex_m2)

intersecao_hex_uc <- suppressWarnings(sf::st_intersection(hex_protecao_raridade, ucs_union_proj))

if (nrow(intersecao_hex_uc) > 0) {
  area_protegida_hex <- intersecao_hex_uc |>
    dplyr::mutate(area_uc_m2 = as.numeric(sf::st_area(geometry))) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(id_hex) |>
    dplyr::summarise(area_uc_m2 = sum(area_uc_m2, na.rm = TRUE), .groups = "drop")
} else {
  area_protegida_hex <- tibble::tibble(id_hex = integer(), area_uc_m2 = numeric())
}

hex_protecao_raridade <- hex_protecao_raridade |>
  dplyr::left_join(area_total_hex, by = "id_hex") |>
  dplyr::left_join(area_protegida_hex, by = "id_hex") |>
  dplyr::mutate(area_uc_m2 = dplyr::coalesce(area_uc_m2, 0), proporcao_protegida = area_uc_m2 / area_hex_m2) |>
  dplyr::left_join(dados_hex_raridade, by = "id_hex") |>
  dplyr::mutate(n_especies_mais_raras = dplyr::coalesce(n_especies_mais_raras, 0L))

# =========================================================
# 23. HEXÁGONOS PRIORITÁRIOS
# =========================================================
hex_prioritarios_raros <- hex_protecao_raridade |>
  dplyr::filter(n_especies_mais_raras > 0, proporcao_protegida <= 1e-10)

cat("\nNúmero de hexágonos potenciais de conservação:", nrow(hex_prioritarios_raros), "\n")

# =========================================================
# 24. ESPÉCIES EM CADA HEXÁGONO PRIORITÁRIO
# =========================================================
especies_por_hex_prioritario <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(id_hex %in% hex_prioritarios_raros$id_hex, !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::filter(scientificName_det %in% especies_mais_raras) |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(n_especies = dplyr::n(), especies = paste(sort(scientificName_det), collapse = "; "), .groups = "drop")

tabela_prioridades_raridade <- hex_prioritarios_raros |>
  sf::st_drop_geometry() |>
  dplyr::select(id_hex, n_especies_mais_raras, proporcao_protegida) |>
  dplyr::left_join(especies_por_hex_prioritario, by = "id_hex") |>
  dplyr::arrange(dplyr::desc(n_especies_mais_raras))

print(tabela_prioridades_raridade)
write.csv(tabela_prioridades_raridade, file.path(dir_tabelas, "hexagonos_prioritarios_especies_raras.csv"), row.names = FALSE)

# Tabela espécie × hexágono
tabela_especie_hex_prioritario <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(id_hex %in% hex_prioritarios_raros$id_hex, scientificName_det %in% especies_mais_raras) |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::left_join(raridade_adaptada |>
                     dplyr::select(scientificName_det, n_hexagonos, n_habitats, n_registros, categoria_raridade), by = "scientificName_det") |>
  dplyr::arrange(id_hex, scientificName_det)

write.csv(tabela_especie_hex_prioritario, file.path(dir_tabelas, "especie_hexagono_prioritario_raridade.csv"), row.names = FALSE)

# =========================================================
# 25. MAPA DOS HEXÁGONOS PRIORITÁRIOS
# =========================================================
mapa_prioridade_raridade <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray40", linewidth = 0.3) +
  ggplot2::geom_sf(data = ucs_union_map, fill = NA, color = "darkgreen", linewidth = 0.25) +
  ggplot2::geom_sf(data = estados_br, fill = NA, color = "gray20", linewidth = 0.25) +
  ggplot2::geom_sf(data = sf::st_transform(hex_prioritarios_raros, 4326), fill = "red3", color = "black", linewidth = 0.25) +
  ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "B. Potential conservation gaps for the rarest species", subtitle = "Hexagons containing rare species and entirely outside protected areas", x = NULL, y = NULL) +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray88", linewidth = 0.3), panel.grid.minor = ggplot2::element_line(color = "gray94", linewidth = 0.2), panel.background = ggplot2::element_rect(fill = "white", color = NA))
print(mapa_prioridade_raridade)
ggplot2::ggsave(file.path(dir_mapas, "mapa_prioridades_especies_raras.png"), mapa_prioridade_raridade, width = 10, height = 7, dpi = 300)

# =========================================================
# 26. MAPA DAS ESPÉCIES MAIS RARAS E NÃO PROTEGIDAS
# =========================================================
dados_hex_raras_UPS <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::inner_join(especies_mais_raras_UPS |>
                      dplyr::select(scientificName_det), by = "scientificName_det") |>
  dplyr::count(id_hex, name = "n_raras_UPS")

hex_raras_UPS <- hex_50km |>
  dplyr::left_join(dados_hex_raras_UPS, by = "id_hex") |>
  dplyr::mutate(n_raras_UPS = dplyr::coalesce(n_raras_UPS, 0L))

mapa_raras_UPS <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = sf::st_transform(brasil_proj, 4326), fill = "gray95", color = "gray55", linewidth = 0.25) +
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
# 27. PAINEL FINAL
# =========================================================
painel_raridade_conservacao <- mapa_raridade | mapa_prioridade_raridade
print(painel_raridade_conservacao)
ggplot2::ggsave(file.path(dir_saida, "painel_raridade_conservacao.png"), painel_raridade_conservacao, width = 16, height = 7, dpi = 300)

# =========================================================
# 28. SALVAR OBJETOS
# =========================================================
save(raridade_protecao, resumo_raridade_protecao, tabela_raridade_protecao, teste_chi, teste_fisher, mais_raras, resumo_mais_raras, representacao_distribuicao, representacao_todas, resumo_MCP, resumo_cobertura_raridade, especies_mais_raras_UPS, dados_hex_raridade, hex_raridade, hex_protecao_raridade, hex_prioritarios_raros, especies_por_hex_prioritario, tabela_prioridades_raridade, tabela_especie_hex_prioritario, hex_raras_UPS, mapa_raridade, mapa_prioridade_raridade, mapa_raras_UPS, painel_raridade_conservacao, file = file.path("outputs", "objetos", "raridade_conservacao.RData"))

# =========================================================
# 29. RESUMO FINAL
# =========================================================
cat("\n=========================================================\nRARIDADE × CONSERVAÇÃO\n=========================================================\n")
cat("\nEspécies na categoria mais rara:", length(especies_mais_raras), "\n")
cat("\nHexágonos contendo espécies mais raras e 100% fora de UCs:", nrow(hex_prioritarios_raros), "\n")
cat("\n--- Proteção das espécies mais raras ---\n")
print(resumo_mais_raras)
cat("\n--- Cobertura da distribuição das espécies mais raras ---\n")
print(resumo_MCP)
cat("\n--- Teste de associação entre raridade e proteção ---\n")
print(teste_chi)
if (!is.null(teste_fisher)) print(teste_fisher)
cat("\n=========================================================\nANÁLISE FINALIZADA\n=========================================================\n")
sf::sf_use_s2(TRUE)
