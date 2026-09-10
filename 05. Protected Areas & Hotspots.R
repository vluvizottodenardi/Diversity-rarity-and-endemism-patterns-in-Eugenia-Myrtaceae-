# =========================================================
# 05 — CONSERVAÇÃO: UCs, HOTSPOTS E PROTEÇÃO POR PONTOS
# =========================================================
library(dplyr); library(sf); library(ggplot2); library(patchwork)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_desc <- file.path("outputs", "objetos", "descritivos_hex.RData")
arquivo_end  <- file.path("outputs", "objetos", "endemismo_hex.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_desc)) stop("Execute primeiro 01_descritivos_esforco_riqueza.R")
if (!file.exists(arquivo_end))  stop("Execute primeiro 02_endemismo_WE_CWE.R")
load(arquivo_base); load(arquivo_desc); load(arquivo_end)

sf::sf_use_s2(FALSE)

# =========================================================
# 1. PARÂMETROS
# =========================================================
arquivo_ucs <- "C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/Uc_Brazil/UC_Brazil.shp"
percentil_hotspot <- 0.90
# Critérios de proteção testados na análise de sensibilidade (Araújo, 2004):
# 0.20            = hotspot = hexágono com MENOS DE 20% de área coberta por UC (critério permissivo)
# "totalmente_fora" = hotspot = hexágono com 0% de área coberta por UC (critério mais conservador)
criterios_protecao <- list(0.20, "totalmente_fora")
# Mantém células com ao menos 1 espécie; interpretar CWE de células muito pobres com cautela.
riqueza_minima_cwe <- 1

dir_conservacao <- file.path("outputs", "conservacao")
dir_tabelas  <- file.path(dir_conservacao, "tabelas")
dir_graficos <- file.path(dir_conservacao, "graficos")
dir_vetores  <- file.path(dir_conservacao, "vetores")
invisible(lapply(list(dir_tabelas, dir_graficos, dir_vetores), dir.create, recursive = TRUE, showWarnings = FALSE))

# =========================================================
# 2. UNIDADES DE CONSERVAÇÃO
# =========================================================
ucs <- sf::st_read(arquivo_ucs, quiet = TRUE) |>
  dplyr::filter(ISO3 == "BRA", REALM != "Marine", STATUS %in% c("Designated", "Established", "Inscribed")) |>
  sf::st_make_valid() |>
  dplyr::filter(!sf::st_is_empty(geometry))

ucs_hex <- ucs |> sf::st_transform(sf::st_crs(hex_50km))
ucs_union <- ucs_hex |> sf::st_union() |> sf::st_make_valid()
ucs_union <- sf::st_sf(id_uc_union = 1, geometry = ucs_union)

# =========================================================
# 3. PROPORÇÃO DE ÁREA PROTEGIDA POR HEXÁGONO
# =========================================================
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
    completamente_fora_uc = proporcao_protegida == 0
  )
# NOTA: "protegido" (acima de um limiar) não é mais um valor fixo — é calculado
# dinamicamente na Seção 6, para cada limiar testado na análise de sensibilidade.

# =========================================================
# 4. JUNTAR PROTEÇÃO À RIQUEZA E AO ENDEMISMO
# =========================================================
gap_riqueza <- hex_riqueza |>
  dplyr::select(id_hex, riqueza) |>
  dplyr::left_join(hex_protecao |> sf::st_drop_geometry() |> dplyr::select(id_hex, proporcao_protegida, percentual_protegido, completamente_fora_uc), by = "id_hex")

gap_endemismo <- hex_endemismo |>
  dplyr::select(id_hex, riqueza, WE, CWE) |>
  dplyr::left_join(hex_protecao |> sf::st_drop_geometry() |> dplyr::select(id_hex, proporcao_protegida, percentual_protegido, completamente_fora_uc), by = "id_hex")

# Limiares de percentil das métricas (não dependem do limiar de proteção)
limiar_riqueza <- stats::quantile(gap_riqueza$riqueza, percentil_hotspot, na.rm = TRUE, names = FALSE)
limiar_cwe <- stats::quantile(
  gap_endemismo$CWE[!is.na(gap_endemismo$CWE) & !is.na(gap_endemismo$riqueza) & gap_endemismo$riqueza >= riqueza_minima_cwe],
  percentil_hotspot, na.rm = TRUE, names = FALSE
)

# =========================================================
# 5. PROTEÇÃO DAS ESPÉCIES, BASEADA NOS PONTOS DE OCORRÊNCIA
# =========================================================
dados_uc <- dados_sf |>
  dplyr::filter(!is.na(scientificName_det), scientificName_det != "", !sf::st_is_empty(geometry)) |>
  sf::st_transform(sf::st_crs(ucs_hex)) |>
  dplyr::mutate(registro_em_uc = lengths(sf::st_intersects(geometry, ucs_union)) > 0)

protecao_especies <- dados_uc |>
  sf::st_drop_geometry() |>
  dplyr::group_by(scientificName_det) |>
  dplyr::summarise(
    n_registros = dplyr::n(), n_registros_em_uc = sum(registro_em_uc, na.rm = TRUE),
    proporcao_registros_em_uc = n_registros_em_uc / n_registros,
    percentual_registros_em_uc = 100 * proporcao_registros_em_uc, .groups = "drop"
  ) |>
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

# =========================================================
# 6. FUNÇÃO GENÉRICA: MAPA DE HOTSPOT, DADO UM LIMIAR DE PROTEÇÃO
# =========================================================
brasil_plot  <- brasil_proj |> sf::st_transform(4326)
estados_plot <- estados_br  |> sf::st_transform(4326)
ucs_plot     <- ucs_hex     |> sf::st_transform(4326)

gerar_mapa_hotspot <- function(dados, coluna_metrica, limiar_metrica, criterio_protecao, riqueza_min = NULL, titulo_metrica) {
  
  # criterio_protecao: um número (ex: 0.20 -> "menos de 20% protegido") OU
  # a string "totalmente_fora" (-> hexágono com 0% de área protegida, completamente_fora_uc)
  if (identical(criterio_protecao, "totalmente_fora")) {
    gap <- dados |> dplyr::mutate(eh_hotspot_protecao = completamente_fora_uc)
    rotulo_protecao <- "cells with 0% protected coverage"
  } else {
    gap <- dados |> dplyr::mutate(eh_hotspot_protecao = proporcao_protegida < criterio_protecao)
    rotulo_protecao <- paste0("cells with <", round(criterio_protecao * 100), "% protected coverage")
  }
  
  hotspots <- if (!is.null(riqueza_min)) {
    gap |> dplyr::filter(!is.na(.data[[coluna_metrica]]), riqueza >= riqueza_min, .data[[coluna_metrica]] >= limiar_metrica, eh_hotspot_protecao)
  } else {
    gap |> dplyr::filter(!is.na(.data[[coluna_metrica]]), .data[[coluna_metrica]] >= limiar_metrica, eh_hotspot_protecao)
  }
  
  hotspots_plot <- hotspots |> sf::st_transform(4326)
  
  mapa <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = brasil_plot, fill = "gray97", color = "gray45", linewidth = 0.35) +
    ggplot2::geom_sf(data = estados_plot, fill = NA, color = "gray65", linewidth = 0.25) +
    ggplot2::geom_sf(data = ucs_plot, fill = "darkgreen", color = "darkgreen", linewidth = 0.25) +
    ggplot2::geom_sf(data = hotspots_plot, fill = "red3", color = "black", linewidth = 0.2) +
    ggplot2::coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Priority gaps in protection",
      subtitle = paste0(rotulo_protecao, ", upper ", round((1 - percentil_hotspot) * 100, 1), "% of ", titulo_metrica)
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 10),
                   plot.margin = ggplot2::margin(10, 10, 10, 10))
  
  list(mapa = mapa, hotspots = hotspots, n_hotspots = nrow(hotspots))
}

# =========================================================
# 7. GERAR OS 4 MAPAS (riqueza x 2 limiares, CWE x 2 limiares)
# =========================================================
res_riqueza_02 <- gerar_mapa_hotspot(gap_riqueza, "riqueza", limiar_riqueza, criterios_protecao[[1]], titulo_metrica = "observed richness")
res_riqueza_10 <- gerar_mapa_hotspot(gap_riqueza, "riqueza", limiar_riqueza, criterios_protecao[[2]], titulo_metrica = "observed richness")

res_cwe_02 <- gerar_mapa_hotspot(gap_endemismo, "CWE", limiar_cwe, criterios_protecao[[1]], riqueza_min = riqueza_minima_cwe, titulo_metrica = "CWE")
res_cwe_10 <- gerar_mapa_hotspot(gap_endemismo, "CWE", limiar_cwe, criterios_protecao[[2]], riqueza_min = riqueza_minima_cwe, titulo_metrica = "CWE")

cat("Hotspots de riqueza — <20% protegido:", res_riqueza_02$n_hotspots, "| totalmente fora (0%):", res_riqueza_10$n_hotspots, "\n")
cat("Hotspots de CWE — <20% protegido:", res_cwe_02$n_hotspots, "| totalmente fora (0%):", res_cwe_10$n_hotspots, "\n")

# =========================================================
# 8. PAINÉIS COMPARATIVOS (SENSIBILIDADE AO LIMIAR)
# =========================================================
painel_sensibilidade_riqueza <- res_riqueza_02$mapa + res_riqueza_10$mapa +
  patchwork::plot_annotation(tag_levels = "A", title = "Sensitivity to protection threshold — Species richness")

painel_sensibilidade_cwe <- res_cwe_02$mapa + res_cwe_10$mapa +
  patchwork::plot_annotation(tag_levels = "A", title = "Sensitivity to protection threshold — CWE")

print(painel_sensibilidade_riqueza)
print(painel_sensibilidade_cwe)

ggplot2::ggsave(file.path(dir_graficos, "sensibilidade_limiar_riqueza.png"), painel_sensibilidade_riqueza, width = 14, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(dir_graficos, "sensibilidade_limiar_CWE.png"), painel_sensibilidade_cwe, width = 14, height = 7, dpi = 300, bg = "white")

# =========================================================
# 9. EXPORTAR TABELAS E CAMADAS DAS 4 VERSÕES DE HOTSPOT
# (enriquecidas com estado dominante e lista de espécies)
# =========================================================

# ---------------------------------------------------------
# 9.0 Função de categoria dominante por área (sempre redefinida
# aqui para garantir a versão corrigida, com st_transform do CRS)
# ---------------------------------------------------------
categoria_dominante_hex <- function(hex_sf, camada_sf, coluna_categoria) {
  camada_sf <- sf::st_transform(camada_sf, sf::st_crs(hex_sf))
  hex_valid    <- hex_sf    |> sf::st_make_valid() |> dplyr::filter(!sf::st_is_empty(geometry))
  camada_valid <- camada_sf |> sf::st_make_valid() |> dplyr::filter(!sf::st_is_empty(geometry))
  
  inter <- sf::st_intersection(hex_valid, camada_valid[, coluna_categoria]) |>
    dplyr::mutate(area_fragmento_m2 = as.numeric(sf::st_area(geometry)))
  
  area_total <- hex_sf |>
    dplyr::mutate(area_hex_m2 = as.numeric(sf::st_area(geometry))) |>
    sf::st_drop_geometry() |> dplyr::select(id_hex, area_hex_m2)
  
  inter |>
    sf::st_drop_geometry() |>
    dplyr::left_join(area_total, by = "id_hex") |>
    dplyr::mutate(proporcao = area_fragmento_m2 / area_hex_m2) |>
    dplyr::group_by(id_hex) |>
    dplyr::slice_max(proporcao, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(id_hex, dplyr::all_of(coluna_categoria), proporcao)
}

# ---------------------------------------------------------
# 9.1 Estado, bioma e domínio vegetacional dominantes por
# hexágono (mesmo critério de área para as três camadas)
# ---------------------------------------------------------
estado_dominante <- categoria_dominante_hex(hex_50km, estados_br, "name_state") |>
  dplyr::rename(estado = name_state)

bioma_dominante_hotspot <- categoria_dominante_hex(hex_50km, biomas_proj, "bioma")

dominio_dominante_hotspot <- categoria_dominante_hex(hex_50km, vegetacao_proj, "habitat") |>
  dplyr::rename(dominio = habitat)

# ---------------------------------------------------------
# 9.2 Lista de espécies por hexágono (nome e contagem)
# ---------------------------------------------------------
especies_por_hex <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::group_by(id_hex) |>
  dplyr::summarise(
    n_especies_hex = dplyr::n(),
    especies_hex = paste(sort(scientificName_det), collapse = "; "),
    .groups = "drop"
  )

# ---------------------------------------------------------
# 9.3 Função para enriquecer uma tabela de hotspot com estado,
# bioma, domínio vegetacional e lista de espécies
# ---------------------------------------------------------
enriquecer_hotspot <- function(hotspots_sf) {
  hotspots_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(id_hex = as.character(id_hex)) |>
    dplyr::left_join(estado_dominante |> dplyr::mutate(id_hex = as.character(id_hex)) |> dplyr::select(id_hex, estado), by = "id_hex") |>
    dplyr::left_join(bioma_dominante_hotspot |> dplyr::mutate(id_hex = as.character(id_hex)) |> dplyr::select(id_hex, bioma), by = "id_hex") |>
    dplyr::left_join(dominio_dominante_hotspot |> dplyr::mutate(id_hex = as.character(id_hex)) |> dplyr::select(id_hex, dominio), by = "id_hex") |>
    dplyr::left_join(especies_por_hex, by = "id_hex") |>
    dplyr::relocate(estado, bioma, dominio, n_especies_hex, especies_hex, .after = id_hex)
}

utils::write.csv(enriquecer_hotspot(res_riqueza_02$hotspots), file.path(dir_tabelas, "hotspots_riqueza_menos20pct.csv"), row.names = FALSE)
utils::write.csv(enriquecer_hotspot(res_riqueza_10$hotspots), file.path(dir_tabelas, "hotspots_riqueza_totalmente_fora.csv"), row.names = FALSE)
utils::write.csv(enriquecer_hotspot(res_cwe_02$hotspots), file.path(dir_tabelas, "hotspots_CWE_menos20pct.csv"), row.names = FALSE)
utils::write.csv(enriquecer_hotspot(res_cwe_10$hotspots), file.path(dir_tabelas, "hotspots_CWE_totalmente_fora.csv"), row.names = FALSE)

sf::st_write(res_riqueza_02$hotspots, file.path(dir_vetores, "hotspots_sensibilidade.gpkg"), layer = "riqueza_menos20pct", delete_layer = TRUE, quiet = TRUE)
sf::st_write(res_riqueza_10$hotspots, file.path(dir_vetores, "hotspots_sensibilidade.gpkg"), layer = "riqueza_totalmente_fora", delete_layer = TRUE, quiet = TRUE)
sf::st_write(res_cwe_02$hotspots, file.path(dir_vetores, "hotspots_sensibilidade.gpkg"), layer = "cwe_menos20pct", delete_layer = TRUE, quiet = TRUE)
sf::st_write(res_cwe_10$hotspots, file.path(dir_vetores, "hotspots_sensibilidade.gpkg"), layer = "cwe_totalmente_fora", delete_layer = TRUE, quiet = TRUE)

# =========================================================
# 10. ESPÉCIES PRESENTES NOS HOTSPOTS (usando o limiar principal = 20%)
# =========================================================
hotspots_riqueza <- res_riqueza_02$hotspots
hotspots_cwe <- res_cwe_02$hotspots
ids_hotspots_riqueza <- hotspots_riqueza$id_hex
ids_hotspots_cwe <- hotspots_cwe$id_hex

especies_hotspots_riqueza <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(id_hex %in% ids_hotspots_riqueza, !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::left_join(protecao_especies |> dplyr::select(scientificName_det, status_protecao_pontos, n_registros, n_registros_em_uc, percentual_registros_em_uc), by = "scientificName_det") |>
  dplyr::arrange(id_hex, scientificName_det)

especies_hotspots_cwe <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(id_hex %in% ids_hotspots_cwe, !is.na(scientificName_det), scientificName_det != "") |>
  dplyr::distinct(id_hex, scientificName_det) |>
  dplyr::left_join(protecao_especies |> dplyr::select(scientificName_det, status_protecao_pontos, n_registros, n_registros_em_uc, percentual_registros_em_uc), by = "scientificName_det") |>
  dplyr::arrange(id_hex, scientificName_det)

print(especies_hotspots_riqueza); print(especies_hotspots_cwe)

utils::write.csv(especies_hotspots_riqueza, file.path(dir_tabelas, "especies_hotspots_riqueza.csv"), row.names = FALSE)
utils::write.csv(especies_hotspots_cwe, file.path(dir_tabelas, "especies_hotspots_CWE.csv"), row.names = FALSE)

# =========================================================
# 11. MAPA CONJUNTO — HOTSPOTS DE RIQUEZA E CWE (limiar principal = 20%)
# =========================================================
mapas_hotspots <- res_riqueza_02$mapa + res_cwe_02$mapa + patchwork::plot_annotation(tag_levels = "A")
print(mapas_hotspots)
ggplot2::ggsave(file.path(dir_graficos, "mapas_hotspots_riqueza_CWE.png"), mapas_hotspots, width = 14, height = 7, dpi = 300, bg = "white")

# =========================================================
# 12. GRÁFICO — STATUS DE PROTEÇÃO POR CATEGORIA DE RARIDADE
# =========================================================
dados_raridade_protecao <- raridade_export |>
  dplyr::select(scientificName_det, classe_raridade) |>
  dplyr::left_join(protecao_especies |> dplyr::select(scientificName_det, status_protecao_pontos), by = "scientificName_det") |>
  dplyr::filter(!is.na(classe_raridade), !is.na(status_protecao_pontos))

resumo_raridade_protecao <- dados_raridade_protecao |>
  dplyr::count(classe_raridade, status_protecao_pontos, name = "n") |>
  dplyr::group_by(classe_raridade) |>
  dplyr::mutate(total_categoria = sum(n), proporcao = n / total_categoria) |>
  dplyr::ungroup()

grafico_raridade_protecao <- ggplot2::ggplot(
  resumo_raridade_protecao,
  ggplot2::aes(x = classe_raridade, y = proporcao, fill = status_protecao_pontos)
) +
  ggplot2::geom_col(position = "fill", width = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = n), position = ggplot2::position_fill(vjust = 0.5), size = 3) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Protection status across multidimensional rarity categories",
    x = "Rarity category", y = "Proportion of species", fill = "Protection status"
  )

print(grafico_raridade_protecao)
ggplot2::ggsave(file.path(dir_graficos, "grafico_raridade_status_protecao.png"), grafico_raridade_protecao, width = 10, height = 6, dpi = 300)

# =========================================================
# 13. SALVAR OBJETOS
# =========================================================
save(ucs, ucs_hex, ucs_union, hex_protecao, gap_riqueza, gap_endemismo,
     res_riqueza_02, res_riqueza_10, res_cwe_02, res_cwe_10,
     hotspots_riqueza, hotspots_cwe, protecao_especies,
     dados_raridade_protecao, resumo_raridade_protecao,
     file = file.path("outputs", "objetos", "conservacao.RData"))

sf::sf_use_s2(TRUE)
