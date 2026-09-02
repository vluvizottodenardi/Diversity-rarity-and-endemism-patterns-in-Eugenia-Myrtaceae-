# =========================================================
# 08 — PRIORIZAÇÃO SISTEMÁTICA COM prioritizr — 4 CENÁRIOS
# =========================================================
library(dplyr)
library(sf)
library(ggplot2)
library(prioritizr)
library(patchwork)
library(geobr)

arquivo_base <- file.path("outputs", "objetos", "base_analises.RData")
arquivo_comp <- file.path("outputs", "objetos", "composicao.RData")
arquivo_cons <- file.path("outputs", "objetos", "conservacao.RData")
if (!file.exists(arquivo_base)) stop("Execute primeiro 00_preparacao_base.R")
if (!file.exists(arquivo_comp)) stop("Execute primeiro 06_beta_diversidade_NMDS.R")
if (!file.exists(arquivo_cons)) stop("Execute primeiro 04_conservacao_UCs_pontos_hotspots.R")
load(arquivo_base)
load(arquivo_comp)
load(arquivo_cons)

# ---------------------------------------------------------
# CENÁRIOS
# ---------------------------------------------------------
cenarios_prioritizr <- tibble::tribble(
  ~cenario, ~meta, ~limiar_uc,
  "A — Less conservative", 0.20, 0.25,
  "B — Intermediate",     0.30, 0.50,
  "C — Conservative",     0.40, 0.75,
  "D — Most conservative", 0.50, 1.00)

custo_padrao <- 1
usar_ucs_como_rede_existente <- TRUE

# ---------------------------------------------------------
# ESTADOS
# ---------------------------------------------------------
estados_br <- geobr::read_state(year = 2020, simplified = TRUE) |> sf::st_transform(4326)
# ---------------------------------------------------------
# UNIDADES DE PLANEJAMENTO
# ---------------------------------------------------------
pu_base <- hex_50km |>
  dplyr::filter(id_hex %in% as.integer(rownames(matriz_comunidade))) |>
  dplyr::arrange(match(id_hex, as.integer(rownames(matriz_comunidade))))

matriz_pu <- matriz_comunidade[
  match(pu_base$id_hex, as.integer(rownames(matriz_comunidade))), ,
  drop = FALSE
]

nomes_especies_originais <- colnames(matriz_pu)
nomes_especies_prioritizr <- make.names(nomes_especies_originais, unique = TRUE)
colnames(matriz_pu) <- nomes_especies_prioritizr

pu_base <- dplyr::bind_cols(pu_base,
  as.data.frame(matriz_pu, check.names = FALSE))

# ---------------------------------------------------------
# FUNÇÃO PARA EXECUTAR CADA CENÁRIO
# ---------------------------------------------------------
rodar_prioritizacao <- function(meta, limiar_uc, nome_cenario) {
  
  pu_priorizacao <- pu_base
  pu_priorizacao$custo <- custo_padrao
  
  if (usar_ucs_como_rede_existente && exists("hex_protecao") &&
      "proporcao_protegida" %in% names(hex_protecao)) {
    
    pu_priorizacao <- pu_priorizacao |>
      dplyr::left_join(
        hex_protecao |>
          sf::st_drop_geometry() |>
          dplyr::select(id_hex, proporcao_protegida),
        by = "id_hex") |>
      dplyr::mutate(
        locked_in = dplyr::coalesce(
          proporcao_protegida >= limiar_uc,
          FALSE),
        custo = dplyr::if_else(locked_in, 0, custo)
      )
    
  } else {
    pu_priorizacao$locked_in <- FALSE
  }
  
  problema <- prioritizr::problem(
    pu_priorizacao,
    features = nomes_especies_prioritizr,
    cost_column = "custo") |>
    prioritizr::add_min_set_objective() |>
    prioritizr::add_relative_targets(meta) |>
    prioritizr::add_binary_decisions()
  
  if (usar_ucs_como_rede_existente && any(pu_priorizacao$locked_in)) {
    problema <- problema |>
      prioritizr::add_locked_in_constraints("locked_in")
  }
  
  problema <- problema |>
    prioritizr::add_default_solver(verbose = FALSE)
  
  solucao <- tryCatch(
    solve(problema),
    error = function(e) {
      warning(
        paste(
          "Erro no cenário:", nome_cenario,
          "|", conditionMessage(e)))
      NULL
    })
  
  if (is.null(solucao)) return(NULL)
  
  nome_coluna_solucao <- grep(
    "^solution",
    names(solucao),
    value = TRUE
  )[1]
  
  if (is.na(nome_coluna_solucao)) {
    stop("A coluna de solução do prioritizr não foi encontrada.")}
  
  solucao <- solucao |>
    dplyr::mutate(
      selecionado = .data[[nome_coluna_solucao]] > 0.5,
      cenario = nome_cenario,
      meta = meta,
      limiar_uc = limiar_uc)
  
  solucao
}

# ---------------------------------------------------------
# EXECUTAR OS 4 CENÁRIOS
# ---------------------------------------------------------
solucoes_prioritizr <- lapply(
  seq_len(nrow(cenarios_prioritizr)),
  function(i) {
    rodar_prioritizacao(
      meta = cenarios_prioritizr$meta[i],
      limiar_uc = cenarios_prioritizr$limiar_uc[i],
      nome_cenario = cenarios_prioritizr$cenario[i])}
)

names(solucoes_prioritizr) <- cenarios_prioritizr$cenario

solucoes_prioritizr <- solucoes_prioritizr[
  !vapply(solucoes_prioritizr, is.null, logical(1))
]

# ---------------------------------------------------------
# MAPAS DOS 4 CENÁRIOS
# ---------------------------------------------------------
mapas_prioritizr <- lapply(
  names(solucoes_prioritizr),
  function(nome) {
    
    solucao <- solucoes_prioritizr[[nome]]
    
    ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = estados_br,
        fill = "gray95",
        color = "gray55",
        linewidth = 0.20
      ) +
      ggplot2::geom_sf(
        data = sf::st_transform(
          solucao |> dplyr::filter(selecionado),
          4326
        ),
        ggplot2::aes(fill = locked_in),
        color = "black",
        linewidth = 0.10
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = "red3",
          "TRUE" = "darkgreen"
        ),
        labels = c(
          "FALSE" = "New priority",
          "TRUE" = "Existing protected network"
        ),
        name = NULL
      ) +
      ggplot2::coord_sf(
        xlim = c(-74, -34),
        ylim = c(-34, 6),
        expand = FALSE
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = nome,
        subtitle = paste0(
          "Species target = ",
          unique(solucao$meta) * 100,
          "% | Existing network threshold = ",
          unique(solucao$limiar_uc) * 100,
          "%"
        ),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold")
      )
  }
)

# ---------------------------------------------------------
# PAINEL
# ---------------------------------------------------------
painel_prioritizr <- (
  mapas_prioritizr[[1]] |
    mapas_prioritizr[[2]]
) / (
    mapas_prioritizr[[3]] |
      mapas_prioritizr[[4]]
  )

print(painel_prioritizr)

# ---------------------------------------------------------
# RESUMO DOS RESULTADOS
# ---------------------------------------------------------
resumo_prioritizr <- dplyr::bind_rows(
  lapply(
    solucoes_prioritizr,
    function(x) {
      x |>
        sf::st_drop_geometry() |>
        dplyr::summarise(
          cenario = dplyr::first(cenario),
          meta = dplyr::first(meta),
          limiar_uc = dplyr::first(limiar_uc),
          n_hexagonos_selecionados = sum(selecionado),
          n_hexagonos_novos = sum(selecionado & !locked_in),
          n_hexagonos_rede_existente = sum(selecionado & locked_in)
        )
    }
  )
)

print(resumo_prioritizr)

# ---------------------------------------------------------
# SALVAR
# ---------------------------------------------------------
save(
  solucoes_prioritizr,
  mapas_prioritizr,
  painel_prioritizr,
  resumo_prioritizr,
  file = file.path("outputs", "objetos", "prioritizr_cenarios.RData")
)
