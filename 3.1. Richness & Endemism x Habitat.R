# =========================================================
# 03.1 — COMPOSIÇÃO, ESFORÇO, RIQUEZA E ENDEMISMO POR HABITAT
#
# Perguntas:
# 1. How compositionally distinct are Eugenia assemblages among vegetation domains?
# 2. Are sampling effort, richness and endemism associated with particular vegetation domains?
#
# Metodologia: todo hexágono recebe um único "habitat dominante"
# (vegetação que ocupa a maior proporção da área do hexágono),
# evitando pseudo-replicação de hexágonos de borda entre categorias.
# =========================================================
library(dplyr); library(tidyr); library(sf); library(ggplot2); library(vegan); library(tibble); library(readr); library(FSA)
sf::sf_use_s2(FALSE)

# =========================================================
# 1. DIRETÓRIOS
# =========================================================
dir_objetos <- file.path("outputs", "objetos")
dir_tabelas <- file.path("outputs", "tabelas")
dir_vetores <- file.path("outputs", "vetores")
dir_saida   <- file.path("outputs", "figuras")
invisible(lapply(list(dir_objetos, dir_tabelas, dir_vetores, dir_saida), dir.create, recursive = TRUE, showWarnings = FALSE))

# =========================================================
# 2. CARREGAR A BASE
# =========================================================
load(file.path(dir_objetos, "base_analises.RData"))
ls()
table(dados_sf$habitat, useNA = "ifany")

# =========================================================
# 3. CATEGORIA DOMINANTE POR HEXÁGONO (função genérica)
# Aplicada às três camadas: vegetação, bioma e região
# =========================================================
categoria_dominante_hex <- function(hex_sf, camada_sf, coluna_categoria) {
  # Garantir mesmo CRS antes de qualquer operação espacial
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

categoria_dominante_hex <- function(hex_sf, camada_sf, coluna_categoria) {
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
# Corrige manualmente hexágonos sem interseção com a camada
# (adiciona se não existir, atualiza se já existir)
# ---------------------------------------------------------
corrigir_dominante <- function(df, ids, valor, coluna_categoria) {
  df <- df |> dplyr::mutate(!!coluna_categoria := dplyr::if_else(id_hex %in% ids, valor, .data[[coluna_categoria]]))
  faltantes <- setdiff(ids, df$id_hex)
  if (length(faltantes) > 0) {
    nova <- tibble::tibble(id_hex = faltantes, !!coluna_categoria := valor, proporcao = NA_real_)
    df <- dplyr::bind_rows(df, nova)
  }
  df
}

habitat_dominante <- categoria_dominante_hex(hex_50km, vegetacao_proj, "habitat") |>
  corrigir_dominante(c(187, 6224), "Floresta Ombrófila Densa", "habitat")

bioma_dominante <- categoria_dominante_hex(hex_50km, biomas_proj, "bioma")
# Confira se os mesmos 2 hexágonos também ficaram sem bioma/região —
# se sim, aplique a mesma correção manual com o bioma/região corretos:
bioma_dominante |> dplyr::filter(id_hex %in% c(187, 6224))
bioma_dominante <- bioma_dominante |>
  corrigir_dominante(187, "Amazônia", "bioma") |>
  corrigir_dominante(6224, "Mata Atlântica", "bioma")

# Garantir que ambos os objetos estejam no mesmo CRS
estados_br <- sf::st_transform(estados_br, sf::st_crs(hex_50km))

# Calcular a região dominante em cada hexágono
regiao_dominante <- categoria_dominante_hex(hex_50km, estados_br, "name_region") |> dplyr::rename(regiao = name_region)

# Checagens: cada tabela deve ter exatamente 4222 linhas, sem duplicatas de id_hex
purrr::walk(list(habitat_dominante, bioma_dominante, regiao_dominante), ~ print(nrow(.x)))
purrr::walk(list(habitat_dominante, bioma_dominante, regiao_dominante), ~ print(sum(duplicated(.x$id_hex))))

# ---------------------------------------------------------
# Unir as três categorizações por hexágono
# ---------------------------------------------------------
categorias_hex <- habitat_dominante |> dplyr::select(id_hex, habitat) |>
  dplyr::full_join(bioma_dominante |> dplyr::select(id_hex, bioma), by = "id_hex") |>
  dplyr::full_join(regiao_dominante |> dplyr::select(id_hex, regiao), by = "id_hex") |>
  dplyr::mutate(id_hex = as.character(id_hex))

# =========================================================
# 4. MATRIZ ESPÉCIE × HEXÁGONO (presença/ausência)
# =========================================================
dados_matriz <- dados_hex |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(id_hex), !is.na(scientificName_det)) |>
  dplyr::select(id_hex, scientificName_det) |> dplyr::distinct()

matriz_pa <- dados_matriz |>
  dplyr::mutate(presenca = 1) |>
  tidyr::pivot_wider(names_from = scientificName_det, values_from = presenca, values_fill = 0) |>
  tibble::column_to_rownames("id_hex") |> as.matrix()

dim(matriz_pa)

# =========================================================
# 5. ASSOCIAR HABITAT À MATRIZ (dados_composicao)
# =========================================================
dados_composicao <- habitat_dominante |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::filter(id_hex %in% rownames(matriz_pa)) |>
  dplyr::arrange(match(id_hex, rownames(matriz_pa)))

rownames(matriz_pa) <- as.character(rownames(matriz_pa))
hex_sem_habitat <- setdiff(rownames(matriz_pa), dados_composicao$id_hex)
cat("\nHexágonos na matriz:", nrow(matriz_pa), "| Com habitat:", nrow(dados_composicao), "| Sem habitat:", length(hex_sem_habitat), "\n")

ids_validos <- intersect(rownames(matriz_pa), dados_composicao$id_hex)
matriz_pa <- matriz_pa[ids_validos, , drop = FALSE]
dados_composicao <- dados_composicao[match(ids_validos, dados_composicao$id_hex), , drop = FALSE]
stopifnot(identical(rownames(matriz_pa), dados_composicao$id_hex))

idx <- !is.na(dados_composicao$habitat)
matriz_pa <- matriz_pa[idx, , drop = FALSE]
dados_composicao <- dados_composicao[idx, , drop = FALSE]

dados_composicao$riqueza <- rowSums(matriz_pa)

# =========================================================
# 6. ADICIONAR ESFORÇO AMOSTRAL (n_registros) E ENDEMISMO (WE, CWE)
# =========================================================
esforco_hex <- dados_hex |> sf::st_drop_geometry() |> dplyr::filter(!is.na(id_hex)) |>
  dplyr::count(id_hex, name = "n_registros") |> dplyr::mutate(id_hex = as.character(id_hex))

endemismo_tab <- hex_endemismo |> sf::st_drop_geometry() |> dplyr::select(id_hex, WE, CWE) |>
  dplyr::mutate(id_hex = as.character(id_hex))

dados_completo <- categorias_hex |>
  dplyr::left_join(esforco_hex, by = "id_hex") |>
  dplyr::left_join(dados_composicao |> dplyr::select(id_hex, riqueza), by = "id_hex") |>
  dplyr::left_join(endemismo_tab, by = "id_hex") |>
  dplyr::mutate(n_registros = dplyr::coalesce(n_registros, 0L), riqueza = dplyr::coalesce(riqueza, 0L))

head(dados_completo)

# =========================================================
# 7. FUNÇÕES GENÉRICAS (resumo, teste e gráfico por métrica)
# =========================================================
resumir_metrica <- function(df, coluna, categoria = "habitat") {
  df |> dplyr::filter(!is.na(.data[[coluna]])) |>
    dplyr::group_by(.data[[categoria]]) |>
    dplyr::summarise(n_hexagonos = dplyr::n(), media = mean(.data[[coluna]]), mediana = median(.data[[coluna]]), minimo = min(.data[[coluna]]), maximo = max(.data[[coluna]]), desvio_padrao = sd(.data[[coluna]]), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(media))
}

riqueza_total_por_categoria <- function(categoria) {
  dados_hex |> sf::st_drop_geometry() |>
    dplyr::filter(!is.na(id_hex), !is.na(scientificName_det), scientificName_det != "") |>
    dplyr::distinct(id_hex, scientificName_det) |>
    dplyr::mutate(id_hex = as.character(id_hex)) |>
    dplyr::inner_join(categorias_hex |> dplyr::select(id_hex, dplyr::all_of(categoria)), by = "id_hex") |>
    dplyr::filter(!is.na(.data[[categoria]])) |>
    dplyr::group_by(.data[[categoria]]) |>
    dplyr::summarise(riqueza_total = dplyr::n_distinct(scientificName_det), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(riqueza_total))
}

teste_kruskal_dunn <- function(df, coluna, categoria = "habitat") {
  formula_kw <- stats::as.formula(paste(coluna, "~", categoria))
  list(kruskal = kruskal.test(formula_kw, data = df), dunn = FSA::dunnTest(formula_kw, data = df, method = "bh")$res)
}

grafico_boxplot_metrica <- function(df, coluna, categoria = "habitat", ylab, fill = "gray85") {
  ggplot(df, aes(x = .data[[categoria]], y = .data[[coluna]])) +
    geom_boxplot(outlier.alpha = 0.3, fill = fill, color = "grey20", linewidth = 0.4) +
    theme_classic(base_size = 12) +
    labs(x = "Vegetation domain", y = ylab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"))
}

# =========================================================
# 8. RESUMOS POR HABITAT, BIOMA E REGIÃO — ESFORÇO, RIQUEZA, WE E CWE
# =========================================================
# --- Por região ---
resumo_esforco_regiao <- resumir_metrica(dados_completo, "n_registros", "regiao")
resumo_riqueza_regiao <- resumir_metrica(dados_completo, "riqueza", "regiao")
resumo_WE_regiao      <- resumir_metrica(dados_completo, "WE", "regiao")
resumo_CWE_regiao     <- resumir_metrica(dados_completo, "CWE", "regiao")
riqueza_total_regiao <- riqueza_total_por_categoria("regiao")

# --- Por bioma ---
resumo_esforco_bioma <- resumir_metrica(dados_completo, "n_registros", "bioma")
resumo_riqueza_bioma <- resumir_metrica(dados_completo, "riqueza", "bioma")
resumo_WE_bioma      <- resumir_metrica(dados_completo, "WE", "bioma")
resumo_CWE_bioma     <- resumir_metrica(dados_completo, "CWE", "bioma")
riqueza_total_bioma <- riqueza_total_por_categoria("bioma")

# --- Por habitat (domínio vegetacional) ---
resumo_esforco_habitat <- resumir_metrica(dados_completo, "n_registros", "habitat")
resumo_riqueza_habitat <- resumir_metrica(dados_completo, "riqueza", "habitat")
resumo_WE_habitat      <- resumir_metrica(dados_completo, "WE", "habitat")
resumo_CWE_habitat     <- resumir_metrica(dados_completo, "CWE", "habitat")
riqueza_total_habitat <- riqueza_total_por_categoria("habitat")

dados_comparacao <- bind_rows(
  dados_completo |> 
    filter(regiao == "Sudeste") |> 
    mutate(grupo = "Sudeste"),
  
  dados_completo |> 
    filter(bioma %in% c("Amazônia", "Caatinga")) |> 
    mutate(grupo = "Amazon + Caatinga")
) |> 
  group_by(grupo) |> 
  summarise(
    n_hexagonos = n(),
    n_hex_mais_5 = sum(n_registros > 5, na.rm = TRUE),
    porcentagem_mais_5 = 100 * n_hex_mais_5 / n_hexagonos
  ) |> 
  ungroup()

sprintf(
  "This disparity is illustrated by the fact that while %.1f%% of hexagons in the Southeast contained more than 5 records, only %.1f%% did so in the Amazon and Caatinga combined.",
  dados_comparacao$porcentagem_mais_5[dados_comparacao$grupo == "Sudeste"],
  dados_comparacao$porcentagem_mais_5[dados_comparacao$grupo == "Amazon + Caatinga"]
)

# =========================================================
# 9. TESTES DE KRUSKAL-WALLIS + DUNN — ESFORÇO, RIQUEZA, WE E CWE
# =========================================================
categorias <- c("habitat", "bioma", "regiao")
metricas   <- c("n_registros", "riqueza", "WE", "CWE")

todos_testes <- purrr::map(metricas, function(m) {
  purrr::map(categorias, ~ teste_kruskal_dunn(dados_completo, m, .x)) |> purrr::set_names(categorias)
}) |> purrr::set_names(metricas)

# Acesso: todos_testes[[métrica]][[categoria]]$kruskal
print(todos_testes$CWE$bioma$kruskal)
print(todos_testes$CWE$regiao$kruskal)
print(todos_testes$CWE$habitat$kruskal)

# =========================================================
# 10. GRÁFICOS BOXPLOT — ESFORÇO, RIQUEZA, WE E CWE
# =========================================================
grafico_esforco <- grafico_boxplot_metrica(dados_completo, "n_registros", ylab = "Number of records per 50-km hexagon")
grafico_riqueza_habitat <- grafico_boxplot_metrica(dados_completo, "riqueza", ylab = "Species richness per 50-km hexagon")
grafico_WE  <- grafico_boxplot_metrica(dados_completo, "WE",  ylab = "Weighted endemism (WE)")
grafico_CWE <- grafico_boxplot_metrica(dados_completo, "CWE", ylab = "Corrected weighted endemism (CWE)")

grafico_esforco; grafico_riqueza_habitat; grafico_WE; grafico_CWE

# =========================================================
# 11. DISTÂNCIA DE JACCARD, PERMANOVA E DISPERSÃO MULTIVARIADA
# =========================================================
dist_jaccard <- vegan::vegdist(matriz_pa, method = "jaccard", binary = TRUE)

set.seed(123)
permanova <- vegan::adonis2(dist_jaccard ~ habitat, data = dados_composicao, permutations = 999)
print(permanova)
permanova_tabela <- as.data.frame(permanova) |> tibble::rownames_to_column("termo")

dispersion <- vegan::betadisper(dist_jaccard, dados_composicao$habitat)
anova_disp <- anova(dispersion)
set.seed(123)
permutacao_disp <- vegan::permutest(dispersion, permutations = 999)
print(anova_disp); print(permutacao_disp)

distancias_centroide <- tibble::tibble(id_hex = as.character(names(dispersion$distances)), distancia_centroide = as.numeric(dispersion$distances)) |>
  dplyr::left_join(dados_composicao |> dplyr::select(id_hex, habitat), by = "id_hex")
write.csv(distancias_centroide, file.path(dir_tabelas, "dispersao_composicional_Eugenia.csv"), row.names = FALSE)

# =========================================================
# 12. NMDS E DISPERSÃO EM ESPAÇO DE ORDENAÇÃO
# =========================================================
set.seed(123)
nmds <- vegan::metaMDS(matriz_pa, distance = "jaccard", binary = TRUE, k = 2, trymax = 100, autotransform = FALSE)
print(nmds); cat("\nStress do NMDS:", round(nmds$stress, 4), "\n")

nmds_sites <- as.data.frame(vegan::scores(nmds, display = "sites"))
nmds_sites$id_hex <- rownames(nmds_sites)
nmds_sites <- nmds_sites |>
  dplyr::mutate(id_hex = as.character(id_hex)) |>
  dplyr::left_join(dados_composicao |> dplyr::mutate(id_hex = as.character(id_hex)) |> dplyr::select(id_hex, habitat, riqueza), by = "id_hex")
sum(is.na(nmds_sites$habitat)) # o normal é 0

resumo_dispersao <- nmds_sites |>
  dplyr::filter(!is.na(habitat)) |>
  dplyr::group_by(habitat) |>
  dplyr::summarise(n_hexagonos = dplyr::n(), riqueza_media = mean(riqueza, na.rm = TRUE), dispersao_nmds1 = sd(NMDS1, na.rm = TRUE), dispersao_nmds2 = sd(NMDS2, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(dispersao_nmds1))
print(resumo_dispersao)

cor_riqueza_dispersao <- cor.test(resumo_dispersao$riqueza_media, resumo_dispersao$dispersao_nmds1, method = "spearman")
print(cor_riqueza_dispersao)

# =========================================================
# 13. GRÁFICO NMDS FACETADO (nomes em inglês, 1 painel por domínio)
# =========================================================
traducao_dominios <- c("Campinarana" = "Campinarana", "Contato (Ecótono e Encrave)" = "Contact Zones", "Corpo d'água continental" = "Continental Inland Water Bodies", "Estepe" = "Steppe", "Floresta Estacional Decidual" = "Deciduous Seasonal Forest", "Floresta Estacional Semidecidual" = "Semideciduous Seasonal Forest", "Floresta Estacional Sempre-Verde" = "Evergreen Seasonal Forest", "Floresta Ombrófila Aberta" = "Open Ombrophilous Forest", "Floresta Ombrófila Densa" = "Dense Ombrophilous Forest", "Floresta Ombrófila Mista" = "Mixed Ombrophilous Forest", "Formação Pioneira" = "Pioneer Formations", "Savana" = "Savanna", "Savana-Estépica" = "Savanna-Steppe")

nmds_sites_en <- nmds_sites |> dplyr::mutate(habitat_en = dplyr::recode(habitat, !!!traducao_dominios))

hulls_nmds_en <- nmds_sites_en |>
  dplyr::filter(!is.na(habitat_en)) |>
  dplyr::group_by(habitat_en) |>
  dplyr::filter(dplyr::n() >= 3) |>
  dplyr::slice(chull(NMDS1, NMDS2)) |>
  dplyr::ungroup()

fundo_cinza <- nmds_sites_en |> dplyr::select(-habitat, -habitat_en)

ordem_paineis <- resumo_dispersao |> dplyr::mutate(habitat_en = dplyr::recode(habitat, !!!traducao_dominios)) |> dplyr::arrange(dplyr::desc(riqueza_media)) |> dplyr::pull(habitat_en)
nmds_sites_en <- nmds_sites_en |> dplyr::mutate(habitat_en = factor(habitat_en, levels = ordem_paineis))
hulls_nmds_en <- hulls_nmds_en |> dplyr::mutate(habitat_en = factor(habitat_en, levels = ordem_paineis))

grafico_nmds_facetado <- ggplot() +
  geom_point(data = fundo_cinza, aes(x = NMDS1, y = NMDS2), color = "gray88", size = 1.1, alpha = 0.8) +
  geom_polygon(data = hulls_nmds_en, aes(x = NMDS1, y = NMDS2, group = habitat_en), fill = "#2C5C8A", alpha = 0.18, color = "#2C5C8A", linewidth = 0.3) +
  geom_point(data = nmds_sites_en |> dplyr::filter(!is.na(habitat_en)), aes(x = NMDS1, y = NMDS2), color = "#2C5C8A", size = 1.6, alpha = 0.85) +
  facet_wrap(~ habitat_en, ncol = 4) + coord_fixed() +
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS ordination of Eugenia assemblage composition by vegetation domain", subtitle = paste0("Jaccard dissimilarity, stress = ", round(nmds$stress, 4))) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "gray95", color = "gray40", linewidth = 0.3), strip.text = element_text(face = "bold", size = 8.5, color = "gray15"), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "gray92", linewidth = 0.25), panel.spacing = unit(0.6, "lines"), axis.text = element_text(color = "black", size = 8), axis.title = element_text(size = 10.5, face = "bold"), plot.title = element_text(size = 12.5, face = "bold"), plot.subtitle = element_text(size = 10, color = "gray30"), plot.title.position = "plot")

grafico_nmds_facetado
ggsave(file.path(dir_saida, "nmds_dominios_vegetacionais.tiff"), grafico_nmds_facetado, width = 12, height = 10, dpi = 400, bg = "white", compression = "lzw")

# =========================================================
# 14. RELAÇÃO RIQUEZA × CWE
# =========================================================
riqueza_CWE <- dados_completo |> dplyr::select(id_hex, habitat, riqueza, CWE) |> dplyr::filter(!is.na(CWE))
cor_riqueza_CWE <- cor.test(riqueza_CWE$riqueza, riqueza_CWE$CWE, method = "spearman", exact = FALSE)
print(cor_riqueza_CWE)

cor_riqueza_CWE <- cor.test(riqueza_CWE$riqueza, riqueza_CWE$CWE, method = "spearman", exact = FALSE)
grafico_riqueza_CWE <- ggplot(riqueza_CWE, aes(x = riqueza, y = CWE, fill = habitat)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  theme_classic(base_size = 12) +
  labs(x = "Species richness per 50-km hexagon", y = "Corrected weighted endemism (CWE)", fill = "Vegetation domain")

grafico_riqueza_CWE

# =========================================================
# 15. EXPORTAR TABELAS
# =========================================================
write.csv(habitat_dominante, file.path(dir_tabelas, "habitat_dominante_hexagonos.csv"), row.names = FALSE)
write.csv(permanova_tabela, file.path(dir_tabelas, "PERMANOVA_composicao_Eugenia.csv"), row.names = FALSE)
write.csv(nmds_sites, file.path(dir_tabelas, "NMDS_composicao_Eugenia.csv"), row.names = FALSE)
write.csv(resumo_esforco, file.path(dir_tabelas, "esforco_amostral_por_habitat.csv"), row.names = FALSE)
write.csv(resumo_riqueza, file.path(dir_tabelas, "riqueza_Eugenia_por_habitat.csv"), row.names = FALSE)
write.csv(resumo_WE, file.path(dir_tabelas, "WE_por_habitat.csv"), row.names = FALSE)
write.csv(resumo_CWE, file.path(dir_tabelas, "CWE_por_habitat.csv"), row.names = FALSE)
write.csv(teste_esforco$dunn, file.path(dir_tabelas, "dunn_esforco_por_habitat.csv"), row.names = FALSE)
write.csv(teste_riqueza$dunn, file.path(dir_tabelas, "dunn_riqueza_por_habitat.csv"), row.names = FALSE)
write.csv(teste_WE$dunn, file.path(dir_tabelas, "dunn_WE_por_habitat.csv"), row.names = FALSE)
write.csv(teste_CWE$dunn, file.path(dir_tabelas, "dunn_CWE_por_habitat.csv"), row.names = FALSE)
write.csv(riqueza_CWE, file.path(dir_tabelas, "riqueza_CWE_por_hexagono.csv"), row.names = FALSE)

# =========================================================
# 16. EXPORTAR GRÁFICOS
# =========================================================
ggsave(file.path(dir_saida, "esforco_por_habitat.png"), grafico_esforco, width = 8, height = 6, dpi = 300)
ggsave(file.path(dir_saida, "riqueza_por_habitat.png"), grafico_riqueza_habitat, width = 8, height = 6, dpi = 300)
ggsave(file.path(dir_saida, "WE_por_habitat.png"), grafico_WE, width = 8, height = 6, dpi = 300)
ggsave(file.path(dir_saida, "CWE_por_habitat.png"), grafico_CWE, width = 8, height = 6, dpi = 300)
ggsave(file.path(dir_saida, "riqueza_CWE_Eugenia.png"), grafico_riqueza_CWE, width = 8, height = 6, dpi = 300)

# =========================================================
# 17. SALVAR OBJETOS
# =========================================================
save(matriz_pa, dados_composicao, habitat_dominante, hex_50km_habitat, dist_jaccard, permanova, dispersion, nmds, dados_completo, riqueza_CWE, resumo_esforco, resumo_riqueza, resumo_WE, resumo_CWE, teste_esforco, teste_riqueza, teste_WE, teste_CWE, file = file.path(dir_objetos, "analises_composicao_endemismo.RData"))

# =========================================================
# 18. RESUMO FINAL NO CONSOLE
# =========================================================
cat("\n==========================================================\nCOMPOSIÇÃO, ESFORÇO, RIQUEZA E ENDEMISMO POR HABITAT\n==========================================================\n")
cat("\nNúmero de hexágonos:", nrow(matriz_pa), "| Número de espécies:", ncol(matriz_pa), "\n")
cat("\n--- PERMANOVA ---\n"); print(permanova)
cat("\n--- DISPERSÃO MULTIVARIADA ---\n"); print(anova_disp)
cat("\n--- NMDS --- Stress =", round(nmds$stress, 4), "\n")
cat("\n--- KRUSKAL-WALLIS: ESFORÇO × HABITAT ---\n"); print(teste_esforco$kruskal)
cat("\n--- KRUSKAL-WALLIS: RIQUEZA × HABITAT ---\n"); print(teste_riqueza$kruskal)
cat("\n--- KRUSKAL-WALLIS: WE × HABITAT ---\n"); print(teste_WE$kruskal)
cat("\n--- KRUSKAL-WALLIS: CWE × HABITAT ---\n"); print(teste_CWE$kruskal)
cat("\n--- CORRELAÇÃO RIQUEZA × CWE ---\n"); print(cor_riqueza_CWE)
cat("\n==========================================================\nANÁLISES FINALIZADAS\n==========================================================\n")

sf::sf_use_s2(TRUE)