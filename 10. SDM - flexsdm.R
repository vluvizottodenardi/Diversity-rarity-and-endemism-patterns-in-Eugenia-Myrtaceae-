# ======================================================================
# SDM ROBUSTO E ADAPTATIVO COM flexsdm — CENÁRIO PRESENTE
# ======================================================================
# IMPORTANTE:
# - A escolha automática do protocolo é heurística.
# - O usuário pode substituir thinning e validação manualmente.
# - A distância da costa é opcional e entra como variável contínua.
# 1. PACOTES
# ======================================================================
library(readxl)
library(dplyr)
library(sf)
library(terra)
library(spThin)
library(flexsdm)
library(readr)
library(openxlsx)
library(rnaturalearth)
library(ggplot2)
library(tidyterra)
library(viridisLite)
library(patchwork)

sf::sf_use_s2(FALSE)

# ======================================================================
# 2. CONFIGURAÇÕES EDITÁVEIS
# ======================================================================
arquivo_dados <- "Teste.xlsx"
sp_nome <- "Eugenia aff. prasina"

col_especie <- "scientificName_det"
col_longitude <- "longitude.gazetteer"
col_latitude <- "latitude.gazetteer"

# Protocolo
# "automatico": usa regras heurísticas conforme o número de ocorrências, escolhe 
# entre thinning; LOOCV, k-fold ou part_senv.
# "manual": usa as configurações manuais definidas abaixo.

modo_protocolo <- "automatico"
thinning_manual_km <- 5
validacao_manual <- "part_senv"  # "loocv", "kfold" ou "part_senv"
folds_manual <- 5
buffer_manual_km <- 300

# Variáveis climáticas locais
dir_clima <- paste0("C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/", "Variaveis climaticas brasil/Bio 30s")

# Elevação e declividade locais
usar_elevacao <- TRUE
usar_declividade <- FALSE
arquivo_elevacao <- paste0("C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/","Variaveis climaticas brasil/Elev 30s/wc2.1_30s_elev.tif")

# Distância da costa - Recomendável testar para espécies costeiras/restinga.
usar_distancia_costa <- TRUE
# A distância será calculada em quilômetros a partir da linha de costa.
# Ela ficará fora da PCA e entrará diretamente no MaxEnt, em centenas de km.

# Solo
usar_solo <- TRUE
arquivo_solo <- paste0("C:/Users/vinil/OneDrive/Área de Trabalho/USB/Vida Acadêmica/Mestrado/R/R/solo/","soilgrids_brasil_0_5cm.tif")

# Vegetação
# A vegetação é usada como máscara pós-modelagem, não entra na PCA.
usar_mascara_vegetacao <- TRUE
arquivo_vegetacao <- paste0("C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/","vege_area/vege_area.shp")
coluna_vegetacao <- "legenda_1"
arquivo_uc <- "C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/Uc_Brazil/UC_Brazil.shp"
# NULL = detectar classes intersectadas pelas ocorrências.
classes_vegetacao_manuais <- NULL
# Número máximo de classes (tipos de vegetações) detectadas automaticamente antes de emitir alerta.
max_classes_vegetacao_automaticas <- 3

# Background - faz 20% da quantidade de células dentro do M
n_background_max <- 20000
proporcao_background_celulas <- 0.20 #Pode ser que em áreas pequenas seja mais que 20%
n_background_min <- 1000

# MaxEnt
n_cores <- 5
threshold_modelo <- "equal_sens_spec"
metrica_selecao <- "TSS"

# Grade mais conservadora para reduzir sobreajuste.
grid_regmult_padrao <- c(1, 1.5, 2, 2.5, 3, 4, 5)
grid_classes_padrao <- c("l", "lq")

# Extrapolação
calcular_extrapolacao <- TRUE
agregacao_extrapolacao <- 5
truncar_extrapolacao <- FALSE
limiar_extrapolacao <- 50

# ---------------------------------------------------------
# Figuras
# ---------------------------------------------------------
dpi_figuras <- 600
# Suavização SOMENTE para a figura. Nunca é usada na binarização, cálculo de área ou análises.
suavizar_apenas_figura <- FALSE
janela_suavizacao <- 3
# Zoom automático em torno das ocorrências.
margem_zoom_graus <- 5

# ---------------------------------------------------------
# Reprodutibilidade e diretórios
# ---------------------------------------------------------
set.seed(123)

nome_sp_arquivo <- gsub("[^A-Za-z0-9_]+","_",sp_nome)
dir_sdm <- file.path("outputs","SDM_robusto",nome_sp_arquivo)
dir.create(dir_sdm,recursive = TRUE,showWarnings = FALSE)

preencher_na <- function(r){terra::focal(r, w = matrix(1,3,3),
    fun = mean, na.policy = "only", na.rm = TRUE)}

# ======================================================================
# 3. FUNÇÕES AUXILIARES - Rodar só 1x
# ======================================================================
normalizar_nome <- function(x) {x <- trimws(as.character(x))
  x[x == ""] <- NA_character_ 
  x}

escolher_protocolo_automatico <- function(n_occ) {          #função que recebe o n° de ocorrência e decide
  if (n_occ < 5) {stop(paste(
        "A espécie possui apenas",n_occ,"ocorrências utilizáveis.",
        "O workflow automático exige pelo menos 5."))}
  
  if (n_occ <= 20) {
    return(list(
        categoria = "rara",
        thinning_km = 0,
        validacao = "loocv",
        folds = NA_integer_,
        buffer_km = 300,
        descricao = paste(
          "Espécie rara: uma ocorrência por pixel, sem thinning adicional, LOOCV",
          "e buffer inicial de 300 km.")))}
  
  if (n_occ <= 49) {
    return(list(
        categoria = "intermediaria",
        thinning_km = 2,
        validacao = "part_senv",
        folds = 5L,
        buffer_km = 300,
        descricao = paste(
          "Amostra intermediária: thinning de 2 km", "e validação 5-fold.",
          "e buffer inicial de 300 km")))}
  
  if (n_occ <= 200) {
    return(list(
      categoria = "bem_amostrada",
      thinning_km = 5,
      validacao = "part_senv",
      folds = 5L,
      buffer_km = 300,
      descricao = paste(
        "Espécie bem amostrada: thinning de 5 km", "e validação 5-fold.",
        "e buffer inicial de 300 km")))}
  
  list(
    categoria = "super_amostrada",
    thinning_km = 10,
    validacao = "part_senv",
    folds = 5L,
    buffer_km = 300,
    descricao = paste(
      "Espécie super amostrada: thinning de 10 km,","validação ambiental e fallback 5-fold.",
      "e buffer inicial de 300 km"))}

obter_protocolo <- function(n_occ) {
  
  if (modo_protocolo == "automatico") {
    protocolo <- escolher_protocolo_automatico(n_occ)
    protocolo$origem <- "heuristica_automatica"
    return(protocolo)}
  
  if (modo_protocolo != "manual") {
    stop("modo_protocolo deve ser 'automatico' ou 'manual'.")}
  
  if (!validacao_manual %in% c(
    "loocv",
    "kfold",
    "part_senv")) {
    stop("validacao_manual deve ser 'loocv', 'kfold' ou 'part_senv'.")}
  
  list(
    categoria = "manual",
    thinning_km = thinning_manual_km,
    validacao = validacao_manual,
    folds = folds_manual,
    buffer_km = buffer_manual_km,
    descricao = "Protocolo definido manualmente pelo usuário.",
    origem = "manual")}

salvar_raster <- function(raster,nome,datatype = NULL) {
  
  argumentos <- list(
    x = raster,
    filename = file.path(
      dir_sdm,
      nome),
    overwrite = TRUE,
    wopt = list(
      gdal = c("COMPRESS=DEFLATE","BIGTIFF=YES")))
  
  if (!is.null(datatype)) {argumentos$datatype <- datatype}
  
  do.call(terra::writeRaster,argumentos)}

#Ler as variáveis
ler_bioclimaticas <- function(diretorio) {
  
  if (!dir.exists(diretorio)) {
    stop(paste(
        "A pasta climática não foi encontrada:", diretorio))}
  
  arquivos <- list.files(diretorio,
    pattern = "\\.(tif|tiff|asc)$",
    full.names = TRUE,
    ignore.case = TRUE)
  
  if (length(arquivos) == 0) {
    stop("Nenhum raster .tif, .tiff ou .asc foi encontrado na pasta climática.")}
  
  nomes_base <- tolower(basename(arquivos))
  
  numero_bio <- suppressWarnings(
    as.integer(sub(
        ".*bio[_-]?([0-9]{1,2}).*",
        "\\1",
        nomes_base)))
  
  reconhecido <- grepl("bio[_-]?[0-9]{1,2}",
    nomes_base)
  
  selecionar <- reconhecido &
    !is.na(numero_bio) &
    numero_bio >= 1 &
    numero_bio <= 19
  
  arquivos <- arquivos[selecionar]
  numero_bio <- numero_bio[selecionar]
  
  ordem <- order(numero_bio)
  
  arquivos <- arquivos[ordem]
  numero_bio <- numero_bio[ordem]
  
  if (
    length(arquivos) != 19 ||
    !identical(
      numero_bio,
      1:19)) {
    stop(
      paste("As 19 variáveis BIO não foram identificadas corretamente.",
        "Números encontrados:",
        paste(numero_bio, collapse = ", ")))}
  
  bio <- terra::rast(arquivos)
  names(bio) <- paste0("bio",1:19)
  
  tabela_arquivos <- data.frame(
    camada = names(bio),
    arquivo = basename(arquivos))
  print(tabela_arquivos)
  
  #readr::write_csv(tabela_arquivos,
  #file.path(dir_sdm,"arquivos_bioclimaticos_utilizados.csv"))
  bio}

#criar arquivo da distância da costa
criar_distancia_costa <- function(
    referencia,brasil_vect,
    fator_agregacao = 10) {
  
  cat("\nAgregando o raster para acelerar o cálculo...\n")
  
  referencia_grossa <- terra::aggregate(
    referencia, fact = fator_agregacao,
    fun = "mean",
    na.rm = TRUE)
  
  costa <- rnaturalearth::ne_coastline(
    scale = "medium",
    returnclass = "sf")
  
  costa <- sf::st_transform(
    costa, terra::crs(referencia_grossa))
  
  costa_vect <- terra::vect(costa)
  
  costa_raster <- terra::rasterize(
    costa_vect,
    referencia_grossa,
    field = 1,
    background = NA,
    touches = TRUE)
  
  cat("\nCalculando distância da costa...\n")
  
  distancia_grossa_km <- terra::distance(costa_raster) / 1000
  
  cat("\nReamostrando para a resolução original...\n")
  
  distancia_final_km <- terra::resample(
    distancia_grossa_km, referencia,
    method = "bilinear")
  
  distancia_final_km <- terra::mask(distancia_final_km, brasil_vect)
  
  distancia_final_100km <- distancia_final_km / 100
  names(distancia_final_100km) <- "dist_coast_100km"
  
  distancia_final_100km}

harmonizar_raster_continuo <- function(raster,referencia,nome) {
  
  names(raster) <- nome
  
  if (!terra::same.crs(
    raster, referencia)) {
    raster <- terra::project(
      raster,referencia,
      method = "bilinear")}
  
  if (!terra::compareGeom(
    raster,referencia,
    stopOnError = FALSE)) {
    raster <- terra::resample(
      raster,referencia,
      method = "bilinear")}
  
  if (!terra::compareGeom(
    raster,referencia, stopOnError = FALSE)) {
    stop(paste("Não foi possível harmonizar o raster:",nome))}
  
  raster}

criar_figura_publicacao <- function(raster,ocorrencias,arquivo,titulo,
    binario = FALSE,
    zoom = FALSE,
    suavizar = FALSE,
    mostrar_pontos = TRUE) {
    
  # Raster usado apenas na figura
  raster_figura <- raster
  if (suavizar && !binario) {
    janela <- matrix(1,
      janela_suavizacao,
      janela_suavizacao)
    raster_figura <- terra::focal(
      raster_figura,
      w = janela,
      fun = "mean",
      na.policy = "omit",
      fillvalue = NA)}
  
  # Preparar ocorrências e limites administrativos
  ocorrencias_sf <- sf::st_as_sf(ocorrencias,
    coords = c("x", "y"), crs = 4326,
    remove = FALSE)
  brasil_plot <- sf::st_transform(
    brasil, terra::crs(raster_figura))
  estados_plot <- sf::st_transform(
    estados_brasil,terra::crs(raster_figura))
  ocorrencias_sf <- sf::st_transform(
    ocorrencias_sf,terra::crs(raster_figura))
  
  # Limites do zoom
  if (zoom) {
    limite_x <- c(
      min(ocorrencias$x, na.rm = TRUE) - margem_zoom_graus,
      max(ocorrencias$x, na.rm = TRUE) + margem_zoom_graus)
    limite_y <- c(
      min(ocorrencias$y, na.rm = TRUE) - margem_zoom_graus,
      max(ocorrencias$y, na.rm = TRUE) + margem_zoom_graus)
    
  } else {
    
    limite_x <- NULL
    limite_y <- NULL}
  
  # Iniciar gráfico
  grafico <- ggplot()
  
  # Mapa contínuo
  if (!binario) {
    grafico <- grafico +
      tidyterra::geom_spatraster(
        data = raster_figura) +
      scale_fill_viridis_c(
        option = "D", limits = c(0, 1),
        oob = scales::squish, na.value = "white",
        name = "Suitability",
        guide = guide_colorbar(
          title.position = "top",
          barheight = grid::unit(5, "cm"),
          barwidth = grid::unit(0.5, "cm")))
    
  # Mapa binário categórico
  } else {
    
    dados_binarios <- as.data.frame(raster_figura, xy = TRUE, na.rm = FALSE)
    names(dados_binarios)[3] <- "valor"
    dados_binarios$classe <- factor(
      dados_binarios$valor,
      levels = c(0, 1),
      labels = c(
        "Unsuitable",
        "Suitable"))
    
    grafico <- grafico +
      geom_raster(data = dados_binarios,
        aes(
          x = x,
          y = y,
          fill = classe)) +
      scale_fill_manual(
        values = c("Unsuitable" = "purple4", "Suitable" = "yellow2"),
        name = "Suitability",
        drop = FALSE,
        na.value = "white")}
  
  # Estados e contorno do Brasil
  grafico <- grafico + geom_sf(
      data = estados_plot, fill = NA,
      colour = "black", linewidth = 0.18)
  
  # Pontos de ocorrência
  if (mostrar_pontos) {
    grafico <- grafico + geom_sf(data = ocorrencias_sf,
        shape = 21,
        fill = "#D73027",
        colour = "white",
        stroke = 0.3,
        size = 1.4)}
  
  # Aparência geral
  grafico <- grafico +
    coord_sf(
      xlim = limite_x,
      ylim = limite_y,
      expand = FALSE) +
    labs(title = titulo,
      x = NULL,
      y = NULL,
      caption = if (suavizar) {
        paste(
          "Suavização aplicada apenas à figura;",
          "análises realizadas no raster original.")
      } else {
        NULL
      }) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "italic", size = 15, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"))
  
  # Salvar
  ggsave(filename = file.path(dir_sdm, arquivo),
         plot = grafico, width = 8.5, height = 7,
         units = "in", dpi = dpi_figuras, bg = "white")
  return(grafico)}

# ======================================================================
# 4. IMPORTAÇÃO E LIMPEZA DAS OCORRÊNCIAS
# ======================================================================
dados <- readxl::read_xlsx(arquivo_dados)

# Verificar colunas obrigatórias
colunas_obrigatorias <- c(
  col_especie,
  col_longitude,
  col_latitude,
  "recordedBy",
  "recordNumber")

# Criar pasta para tabelas de diagnóstico.
#dir.create(file.path(dir_sdm,"tabelas"), recursive = TRUE, showWarnings = FALSE)

# Padronizar campos e converter coordenadas
dados <- dados |>
  mutate(.especie = normalizar_nome(
      .data[[col_especie]]),
    recordedBy = trimws(
      tolower(as.character(recordedBy))),
    recordNumber = trimws(
      tolower(as.character(recordNumber))),
    .longitude = suppressWarnings(
      as.numeric(.data[[col_longitude]])),
    .latitude = suppressWarnings(
      as.numeric(.data[[col_latitude]]))
  ) |>
  filter(
    !is.na(.especie),
    !is.na(.longitude),
    !is.na(.latitude),
   
    .longitude >= -180,
    .longitude <= 180,
    .latitude >= -90,
    .latitude <= 90)

# ---------------------------------------------------------
# Remover duplicatas de coletor + número de coleta
# ---------------------------------------------------------
# Registros com "s.n." não são agrupados entre si. Cada ocorrência s.n. recebe um identificador próprio.
dados <- dados |>
  mutate(id_coleta = case_when(
      is.na(recordNumber) |
        recordNumber == "" |
        recordNumber %in% c(
          "s.n.",
          "s.n",
          "sn",
          "sem número",
          "sem numero") ~ paste0(.especie,
          "sem_numero_",
          row_number()),
      TRUE ~ paste(
        .especie,
        recordedBy,
        recordNumber,
        sep = "_"))) |>
  group_by(
    id_coleta) |>
  slice(1) |>
  ungroup()

# ---------------------------------------------------------
# Remover coordenadas suspeitas: centroides aproximados do Brasil e dos estados
# ---------------------------------------------------------
centroides_suspeitos <- data.frame(
  localidade = c(
    "Brasil","Acre","Alagoas","Amapa","Amazonas","Bahia",
    "Ceara","Distrito_Federal","Espirito_Santo","Goias",
    "Maranhao","Mato_Grosso","Mato_Grosso_do_Sul","Minas_Gerais",
    "Para","Paraiba","Parana","Pernambuco","Piaui","Rio_de_Janeiro",
    "Rio_Grande_do_Norte","Rio_Grande_do_Sul","Rondonia","Roraima",
    "Santa_Catarina","Sao_Paulo","Sergipe","Tocantins"),
  
  longitude = c(
    -53.08032013,-70.0,-36.5,-52.0,-63.0,-41.5,-39.5,-47.9,
    -40.5,-49.5,-45.0,-56.0,-54.5,-44.0,-52.0,-36.0,-51.5,
    -37.5,-43.0,-43.5,-36.5,-53.0,-63.5,-61.0,-50.5,-48.0,
    -37.5,-48.0),
  
  latitude = c(-10.74257447,-9.0,-9.5,1.0,-4.0,-12.5,-5.0,
    -15.8,-19.5,-16.0,-5.0,-13.0,-20.5,-18.5,-4.0,-7.0,-24.5,
    -8.5,-7.0,-22.5,-5.5,-30.0,-11.0,2.0,-27.0,-22.5,-10.5,-10.0))

dados$coord_centroide <- FALSE
dados$centroide_detectado <- NA_character_

for (i in seq_len(
  nrow(centroides_suspeitos))) {
  
  teste_centroide <-
    abs(dados$.longitude -
        centroides_suspeitos$longitude[i]) <= 0.02 &
    abs(dados$.latitude -
        centroides_suspeitos$latitude[i]) <= 0.02
  
  teste_centroide[
    is.na(teste_centroide)] <- FALSE
  
  dados$coord_centroide <-
    dados$coord_centroide |
    teste_centroide
  
  preencher_nome <- teste_centroide & is.na(dados$centroide_detectado)
  
  dados$centroide_detectado[
    preencher_nome] <- centroides_suspeitos$localidade[i]}

# Salvar registros excluídos por centroides
coordenadas_excluidas <- dados |> filter(coord_centroide)
cat("\nCoordenadas excluídas por coincidência com centroides:",nrow(coordenadas_excluidas),"\n")

# Remover centroides suspeitos.
dados <- dados |>
  filter(!coord_centroide) |>
  select(
    -coord_centroide,
    -centroide_detectado)

# Remover coordenadas duplicadas da mesma espécie
dados <- dados |>
  distinct(
    .especie,
    .longitude,
    .latitude,
    .keep_all = TRUE)

# ---------------------------------------------------------
# Selecionar espécie para o SDM
# ---------------------------------------------------------
sp <- dados |> filter(.especie == sp_nome)
n_occ_inicial <- nrow(sp)
if (n_occ_inicial == 0) {stop(paste("Nenhuma ocorrência foi encontrada para:", sp_nome))}
cat("\nEspécie:",sp_nome,"\nOcorrências válidas após a limpeza:",n_occ_inicial,"\n")

# ======================================================================
# 5. VARIÁVEIS AMBIENTAIS - Rodar apenas na 1x
# ======================================================================
bio <- ler_bioclimaticas(dir_clima)

brasil <- rnaturalearth::ne_countries(
  country = "Brazil",
  scale = "medium",
  returnclass = "sf") |>
  sf::st_make_valid()

estados_brasil <- rnaturalearth::ne_states(country = "Brazil", returnclass = "sf") |> sf::st_make_valid()

#Deixa o brasil e as bios no mesmo "formato"
brasil_vect <- terra::vect(sf::st_transform(brasil, terra::crs(bio)))
bio <- terra::crop(bio,brasil_vect) |> terra::mask(brasil_vect)
preditores_continuos <- bio

if (
  usar_elevacao ||
  usar_declividade) {
  
  if (!file.exists(arquivo_elevacao)) {
    stop(paste("O arquivo local de elevação não foi encontrado:",arquivo_elevacao))}
  
  elev <- terra::rast(arquivo_elevacao)
  elev <- harmonizar_raster_continuo(elev,bio[[1]],"elevation")
  elev <- terra::crop(elev, brasil_vect) |>terra::mask(brasil_vect)
  
  if (usar_elevacao) {preditores_continuos <- c(preditores_continuos,elev)}
  if (usar_declividade) {
    slope <- terra::terrain(elev,
      v = "slope",
      unit = "degrees")
    names(slope) <- "slope"
    preditores_continuos <- c(preditores_continuos,slope)}
  }

arquivo_distancia_costa <- paste0("C:/Users/vinil/OneDrive/Área de Trabalho/","Variaveis climaticas brasil/Distancia costa 30s/","distancia_costa_100km.tif")

if (usar_distancia_costa) {
  
  if (file.exists(arquivo_distancia_costa)) {
    cat("\nLendo distância da costa já existente...\n")
    distancia_costa <- terra::rast(arquivo_distancia_costa)
    
  } else {
    
    dir.create(dirname(arquivo_distancia_costa),recursive = TRUE, showWarnings = FALSE)
    distancia_costa <- criar_distancia_costa(
      referencia = bio[[1]], brasil_vect = brasil_vect,
      fator_agregacao = 10)  #distânica da costa 10x10 pixel pra ser mais rápido
    
    terra::writeRaster(distancia_costa,
      arquivo_distancia_costa,
      overwrite = TRUE,
      wopt = list(
        gdal = c(
          "COMPRESS=DEFLATE",
          "PREDICTOR=3",
          "BIGTIFF=YES")))}
  
  distancia_costa <- harmonizar_raster_continuo(distancia_costa,bio[[1]],"dist_coast_100km")
}

#teste pra ver se a distância da costa esta correta
if (usar_distancia_costa) {
  teste_distancia <- terra::extract(distancia_costa,
    terra::vect(data.frame(
        longitude = c(-46.63, -60.02), #O 1° é SP e o 2° Manaus - se tiver certo o 2° > 1°
        latitude = c(-23.55, -3.12)),
      geom = c("longitude","latitude"),
      crs = "EPSG:4326"))
  print(teste_distancia)}

# ======================================================================
# SOLO — HARMONIZAÇÃO SEM GRANULADO
# ======================================================================
if (usar_solo) {
  
  if (!file.exists(arquivo_solo)) {
    stop(paste("usar_solo = TRUE, mas o arquivo não foi encontrado:",arquivo_solo))
  }
  
  # 1. CARREGAR O SOLO ORIGINAL
  
  solo_original <- terra::rast(arquivo_solo)
  
  if (terra::nlyr(solo_original) != 3) {
    stop(paste0("O arquivo de solo precisa conter exatamente 3 camadas.\n",
        "Foram encontradas: ", terra::nlyr(solo_original)))
  }
  
  names(solo_original) <- c("sand", "clay", "phh2o")
  
  cat(
    "\n============================================================",
    "\nDIAGNÓSTICO DO SOLO ORIGINAL",
    "\n============================================================",
    "\nCRS:\n", terra::crs(solo_original, proj = TRUE),
    "\nResolução:", paste(terra::res(solo_original), collapse = " × "),
    "\nDimensão:", terra::nrow(solo_original),
    "linhas ×", terra::ncol(solo_original),
    "colunas\n")
  
  print(terra::global(solo_original,c("min", "max", "mean", "sd"),na.rm = TRUE))
  
  # 2. VERIFICAR A RELAÇÃO ENTRE AS RESOLUÇÕES
  
  resolucao_solo <- terra::res(solo_original)
  resolucao_clima <- terra::res(bio[[1]])
  
  cat(
    "\nResolução do solo:",paste(resolucao_solo, collapse = " × "),
    "\nResolução do clima:",paste(resolucao_clima, collapse = " × "),
    "\n")
  
  # 3. PROJETAR DIRETAMENTE PARA A GRADE DO CLIMA
  # --------------------------------------------------------------------
  # "average" é mais adequado quando diversos pixels mais finos do solo
  # precisam representar um pixel climático maior.
  # Isso evita que o valor final dependa apenas de uma pequena vizinhança
  # e reduz padrões quadriculados por subamostragem.
  # --------------------------------------------------------------------
  
  solo <- terra::project(solo_original,bio[[1]],method = "average")
  
  if (
    !terra::compareGeom(solo[[1]], bio[[1]], stopOnError = FALSE) ) {
    stop(paste("O solo não ficou alinhado à grade climática após",
        "terra::project(method = 'average')."))
  }
  
  # 4. RECORTAR E MASCARAR
  
  solo <- terra::crop(solo,brasil_vect)
  solo <- terra::mask(solo,brasil_vect)
  
  mascara_solo_valido <- terra::ifel(!is.na(solo[[1]]) & !is.na(solo[[2]]) & !is.na(solo[[3]]),
    1,
    NA
  )
  
  # 5. REDUZIR VARIAÇÃO DE ALTA FREQUÊNCIA

  janela_solo <- matrix(1, nrow = 3, ncol = 3)
  
  solo_media_local <- terra::focal(
    solo,
    w = janela_solo,
    fun = "mean",
    na.rm = TRUE,
    na.policy = "omit")
  
  # Medida da diferença entre o pixel e sua vizinhança.
  diferenca_solo <- abs(solo - solo_media_local)
  
  # Desvio-padrão global de cada camada.
  dp_solo <- terra::global(
    solo,
    "sd",
    na.rm = TRUE)[, 1]
  
  # Criar uma camada por vez para evitar problemas de reciclagem.
  solo_regularizado <- solo
  
  for (i in seq_len(terra::nlyr(solo))) {
    
    # Um pixel é tratado apenas quando difere da média local
    # em mais de meio desvio-padrão da camada.
    limiar_ruido <- 0.5 * dp_solo[i]
    
    solo_regularizado[[i]] <- terra::ifel(
      !is.na(diferenca_solo[[i]]) &
        diferenca_solo[[i]] > limiar_ruido,
      solo_media_local[[i]],
      solo[[i]])
  }
  
  names(solo_regularizado) <- c(
    "sand",
    "clay",
    "phh2o")
  
  # Restaurar a máscara original.
  solo_regularizado <- terra::mask(solo_regularizado, mascara_solo_valido)
  
  solo <- solo_regularizado
  
  # 6. LIMITES FISICAMENTE PLAUSÍVEIS
  
  valores_solo <- terra::global(solo, c("min", "max"), na.rm = TRUE)
  print(valores_solo)
  
  # Detectar automaticamente uma possível escala inteira do SoilGrids.
  # SoilGrids pode armazenar algumas propriedades com fator de escala.
  max_sand <- valores_solo["sand", "max"]
  max_clay <- valores_solo["clay", "max"]
  max_ph <- valores_solo["phh2o", "max"]
  
  # Areia e argila acima de 100 provavelmente estão em g/kg.
  if (is.finite(max_sand) &&
    max_sand > 100) {
    cat(
      "\nA camada sand parece estar em g/kg.",
      "\nConvertendo para porcentagem: sand / 10.\n")
    
    solo[["sand"]] <- solo[["sand"]] / 10}
  
  if (
    is.finite(max_clay) &&
    max_clay > 100) {
    cat(
      "\nA camada clay parece estar em g/kg.",
      "\nConvertendo para porcentagem: clay / 10.\n")
    
    solo[["clay"]] <- solo[["clay"]] / 10}
  
  # pH SoilGrids pode estar armazenado como pH × 10.
  if (
    is.finite(max_ph) &&
    max_ph > 14) {
    cat(
      "\nA camada phh2o parece estar multiplicada por 10.",
      "\nConvertendo para unidades de pH: phh2o / 10.\n")
    
    solo[["phh2o"]] <- solo[["phh2o"]] / 10}
  
  # Garantir limites possíveis depois da conversão.
  solo[["sand"]] <- terra::clamp(
    solo[["sand"]],
    lower = 0,
    upper = 100,
    values = TRUE)
  
  solo[["clay"]] <- terra::clamp(
    solo[["clay"]],
    lower = 0,
    upper = 100,
    values = TRUE)
  
  solo[["phh2o"]] <- terra::clamp(
    solo[["phh2o"]],
    lower = 0,
    upper = 14,
    values = TRUE)
  
  # 7. DIAGNÓSTICO DA REGULARIZAÇÃO
  cat(
    "\n============================================================",
    "\nDIAGNÓSTICO DO SOLO HARMONIZADO",
    "\n============================================================\n")
  
  print(terra::global(solo, c("min", "max", "mean", "sd"), na.rm = TRUE))
  
  cat(
    "\nGeometria igual à BIO:",
    terra::compareGeom(
      solo[[1]],
      bio[[1]],
      stopOnError = FALSE),"\n")
  
  # 8. REMOVER VERSÕES ANTERIORES DO SOLO
  nomes_solo <- c(
    "sand",
    "clay",
    "phh2o",
    "sand_mean_0-5cm",
    "clay_mean_0-5cm",
    "phh2o_mean_0-5cm")
  
  manter_preditores <- !names(preditores_continuos) %in% nomes_solo
  preditores_continuos <- preditores_continuos[[manter_preditores]]
  
  # 10. ADICIONAR O SOLO UMA ÚNICA VEZ
  preditores_continuos <- c(preditores_continuos,solo)
  
  cat(
    "\nCamadas de solo adicionadas:", paste(names(solo), collapse = ", "),
    "\nNúmero total de preditores:", terra::nlyr(preditores_continuos),
    "\nPreditores:",paste(
      names(preditores_continuos),
      collapse = ", "),"\n")
  }

preditores_continuos <- flexsdm::homogenize_na(preditores_continuos)

#conferir se não tem muitos pixels NA
terra::global(!is.na(preditores_continuos[[1]]),"sum",na.rm = TRUE)
resolucao_real <- terra::res(bio)
cat("\nResolução real dos rasters:", resolucao_real[1],"x",resolucao_real[2],"graus\n")
cat("\nPreditores contínuos:",paste(names(preditores_continuos),collapse = ", "),"\n")

# ======================================================================
# 6. FILTRAGEM AMBIENTAL, RECUPERAÇÃO COSTEIRA E FILTRO POR PIXEL
# ======================================================================
# Distância máxima permitida até uma célula ambiental válida.
# Para rasters WorldClim de 30 segundos, 1500 m corresponde aproximadamente a uma ou duas células.
raio_recuperacao_m <- 1500

# 6.1 PREPARAR AS OCORRÊNCIAS
sp_pts <- data.frame(
  id_ocorrencia = seq_len(nrow(sp)),
  scientificName_det = sp_nome,
  recordedBy = sp$recordedBy,
  recordNumber = sp$recordNumber,
  longitude.gazetteer = sp$.longitude,
  latitude.gazetteer = sp$.latitude,
  stringsAsFactors = FALSE)

coords_sp <- as.matrix(sp_pts[,c("longitude.gazetteer","latitude.gazetteer"), drop = FALSE])
storage.mode(coords_sp) <- "numeric"
if (any(!is.finite(coords_sp))) {stop("Existem coordenadas ausentes ou não numéricas em sp_pts.")}
if (nrow(coords_sp) != nrow(sp_pts)) {stop("O número de coordenadas não corresponde ao número de ocorrências.")}

# 6.2 IDENTIFICAR A CÉLULA ORIGINAL
sp_pts$cell_id_original <- terra::cellFromXY(preditores_continuos[[1]],coords_sp)

# 6.3 EXTRAIR OS VALORES AMBIENTAIS ORIGINAIS
valores_originais <- terra::extract(preditores_continuos,coords_sp)
if ("ID" %in% names(valores_originais)) {valores_originais$ID <- NULL}

# Garantir que as colunas correspondem aos preditores
colunas_ausentes_extracao <- setdiff(names(preditores_continuos), names(valores_originais))

if (length(colunas_ausentes_extracao) > 0) {
  stop(paste0("As seguintes variáveis não foram retornadas por terra::extract(): ",
      paste(colunas_ausentes_extracao, collapse = ", ")))}

valores_originais <- valores_originais[,names(preditores_continuos),drop = FALSE]

if (nrow(valores_originais) != nrow(sp_pts)) {stop(paste0(
      "A extração ambiental retornou ",nrow(valores_originais),
      " linhas, mas existem ",nrow(sp_pts)," ocorrências."))}

ambiente_original_completo <- complete.cases(valores_originais)

cat(
  "\nOcorrências antes do filtro ambiental:",nrow(sp_pts),
  "\nOcorrências com ambiente na célula original:", sum(ambiente_original_completo),
  "\nOcorrências inicialmente com NA:", sum(!ambiente_original_completo),"\n")

# 6.4 PREPARAR A RECUPERAÇÃO DOS PONTOS COSTEIROS
indices_na <- which(!ambiente_original_completo)

# Começar com os valores ambientais originais
valores_ambientais_sp <- valores_originais

# Colunas de diagnóstico
sp_pts$recuperada_celula_proxima <- FALSE
sp_pts$distancia_recuperacao_m <- NA_real_
sp_pts$cell_id_recuperada <- NA_real_
sp_pts$longitude_celula_recuperada <- NA_real_
sp_pts$latitude_celula_recuperada <- NA_real_

# 6.5 RECUPERAR SOMENTE AS OCORRÊNCIAS COM NA
if (length(indices_na) > 0) {
  
  # Coordenadas originais dos pontos com NA
  coords_na <- coords_sp[indices_na, ,drop = FALSE]
  
  # Converter as ocorrências problemáticas para SpatVector
  pontos_na_vect <- terra::vect(coords_na, type = "points",
    crs = terra::crs(preditores_continuos))
  
  # Nome da camada utilizada como referência
  nome_camada_referencia <- names(preditores_continuos)[1]
  
  # Procurar a célula válida mais próxima
  busca_proxima <- terra::extract(
    preditores_continuos[[1]],pontos_na_vect,
    ID = FALSE, xy = TRUE,
    search_radius = raio_recuperacao_m)
  
  cat("\nColunas retornadas pela busca:", paste(names(busca_proxima), collapse = ", "),"\n")
  
  # Verificar se as coordenadas da célula encontrada foram retornadas
  if (!all(c("x", "y") %in% names(busca_proxima))) {
    stop(paste0(
        "terra::extract() não retornou as coordenadas ",
        "da célula ambiental encontrada. ",
        "Colunas retornadas: ",
        paste(names(busca_proxima), collapse = ", ")))}
  
  # Verificar se a camada de referência foi retornada
  if (
    !nome_camada_referencia %in%
    names(busca_proxima)) {
    stop(paste0("A camada ambiental de referência não foi retornada. ","Camada esperada: ",
        nome_camada_referencia))}
  
  if (nrow(busca_proxima) != length(indices_na)) {
    stop(paste0("A busca por células próximas retornou ", nrow(busca_proxima),
        " linhas para ", length(indices_na), " ocorrências com NA."))}
  
  # Coordenadas das células encontradas
  coords_encontradas <- as.matrix(busca_proxima[, c("x", "y"), drop = FALSE])
  storage.mode(coords_encontradas) <- "numeric"
  
  # Busca válida:
  # 1. x e y existem;
  # 2. a camada de referência possui valor.
  busca_valida <-
    complete.cases(coords_encontradas) &
    !is.na(busca_proxima[[nome_camada_referencia]])
  
  # Objetos de diagnóstico locais
  celulas_locais <- rep(NA_real_, length(indices_na))
  distancias_locais <- rep(NA_real_, length(indices_na))
  
  if (any(busca_valida)) {
    
    # Identificar o número das células encontradas
    celulas_locais[busca_valida] <-
      terra::cellFromXY(
        preditores_continuos[[1]],
        coords_encontradas[
          busca_valida,
          ,
          drop = FALSE])
    
    # Pontos originais com busca válida
    pontos_originais_validos <- terra::vect(
      coords_na[busca_valida,
        ,
        drop = FALSE],
      type = "points",
      crs = terra::crs(preditores_continuos))
    
    # Centros das células encontradas
    pontos_encontrados_validos <- terra::vect(
      coords_encontradas[
        busca_valida,
        ,
        drop = FALSE],
      type = "points",
      crs = terra::crs(preditores_continuos))
    
    # Matriz de distâncias entre os pontos
    matriz_distancias <- terra::distance(pontos_originais_validos,pontos_encontrados_validos)
    matriz_distancias <- as.matrix(matriz_distancias)
    
    # Os pares correspondentes estão na diagonal
    distancias_locais[busca_valida] <- diag(matriz_distancias)}
  
  # Considerar recuperável somente quando:
  # 1. a busca foi válida;
  # 2. uma célula foi identificada;
  # 3. a distância foi calculada;
  # 4. a distância está dentro do limite.
  recuperavel_local <-
    busca_valida &
    !is.na(celulas_locais) &
    !is.na(distancias_locais) &
    distancias_locais <= raio_recuperacao_m
  
  indices_recuperados <- indices_na[recuperavel_local]
  
  cat(
    "\nPontos inicialmente com NA:",length(indices_na),
    "\nCélulas válidas encontradas:",sum(busca_valida),
    "\nPontos dentro do raio permitido:",sum(recuperavel_local),
    "\nPontos sem recuperação inicial:",sum(!recuperavel_local),"\n")
  
  if (length(indices_recuperados) > 0) {
    
    celulas_recuperadas <- celulas_locais[recuperavel_local]
    
    # Extrair todas as variáveis ambientais das células próximas encontradas
    valores_recuperados <- terra::extract(preditores_continuos,celulas_recuperadas)
    
    if ("ID" %in% names(valores_recuperados)) {valores_recuperados$ID <- NULL}
    
    colunas_ausentes_recuperacao <- setdiff(
      names(preditores_continuos),
      names(valores_recuperados))
    
    if (
      length(colunas_ausentes_recuperacao) > 0) {
      stop( paste0("As seguintes variáveis não foram retornadas ",
          "na recuperação ambiental: ", paste(colunas_ausentes_recuperacao,
            collapse = ", ")))}
    
    valores_recuperados <- valores_recuperados[
      ,
      names(preditores_continuos),
      drop = FALSE]
    
    if (
      nrow(valores_recuperados) !=
      length(indices_recuperados)) {
      stop(paste0(
          "Foram encontradas ", length(indices_recuperados),
          " ocorrências recuperáveis, mas a extração retornou ", nrow(valores_recuperados),
          " linhas."))}
    
    # Verificar se todas as variáveis estão disponíveis
    recuperacao_completa <- complete.cases(valores_recuperados)
    
    # Posições locais dentro do conjunto dos pontos com NA
    posicoes_recuperadas <- which(recuperavel_local)
    
    if (any(!recuperacao_completa)) {
      
      warning(
        sum(!recuperacao_completa),
        paste0(" células próximas ainda possuem NA ",
          "em pelo menos uma variável ambiental."))
      
      indices_recuperados <- indices_recuperados[recuperacao_completa]
      celulas_recuperadas <- celulas_recuperadas[recuperacao_completa]
      valores_recuperados <- valores_recuperados[recuperacao_completa,
        ,drop = FALSE]
      posicoes_recuperadas <- posicoes_recuperadas[recuperacao_completa]
    }
    
    if (length(indices_recuperados) > 0) {
      
      # Inserir os valores ambientais recuperados
      valores_ambientais_sp[
        indices_recuperados,
        names(preditores_continuos)] <- valores_recuperados[
        , names(preditores_continuos),drop = FALSE]
      
      # Registrar as ocorrências recuperadas
      sp_pts$recuperada_celula_proxima[indices_recuperados] <- TRUE
      sp_pts$distancia_recuperacao_m[indices_recuperados] <- distancias_locais[posicoes_recuperadas]
      sp_pts$cell_id_recuperada[indices_recuperados] <- celulas_recuperadas
      
      # Coordenadas dos centros das células usadas
      xy_recuperadas <- terra::xyFromCell(preditores_continuos[[1]], celulas_recuperadas)
      sp_pts$longitude_celula_recuperada[indices_recuperados] <- xy_recuperadas[, 1]
      sp_pts$latitude_celula_recuperada[indices_recuperados] <- xy_recuperadas[, 2]
    }
  }
}

# 6.6 DEFINIR A CÉLULA EFETIVAMENTE UTILIZADA
# Por padrão, utilizar a célula original
sp_pts$cell_id <- sp_pts$cell_id_original

# Para as ocorrências recuperadas, utilizar a célula próxima
indices_celula_recuperada <- which(sp_pts$recuperada_celula_proxima & !is.na(sp_pts$cell_id_recuperada))

if (length(indices_celula_recuperada) > 0) {
  sp_pts$cell_id[indices_celula_recuperada] <- sp_pts$cell_id_recuperada[indices_celula_recuperada]}

# 6.7 FILTRO AMBIENTAL FINAL
ambiente_completo <- complete.cases(valores_ambientais_sp)

if (
  length(ambiente_completo) !=
  nrow(sp_pts)) {stop(paste0(
      "O vetor ambiente_completo possui ",length(ambiente_completo),
      " valores, mas sp_pts possui ",nrow(sp_pts)," linhas."))}

cat(
  "\nOcorrências inicialmente completas:",sum(ambiente_original_completo),
  "\nOcorrências costeiras recuperadas:",sum(sp_pts$recuperada_celula_proxima),
  "\nOcorrências completas após recuperação:",sum(ambiente_completo),
  "\nOcorrências ainda com ambiente incompleto:",sum(!ambiente_completo),
  "\nOcorrências sem cell_id válido:",sum(is.na(sp_pts$cell_id)),"\n")

# 6.8 SALVAR DIAGNÓSTICOS
ocorrencias_recuperadas_ambiente <- sp_pts[sp_pts$recuperada_celula_proxima, ,drop = FALSE]

ocorrencias_excluidas_ambiente <- sp_pts[!ambiente_completo | is.na(sp_pts$cell_id),
  ,drop = FALSE]

if (
  nrow(ocorrencias_recuperadas_ambiente) > 0) {
  readr::write_csv(ocorrencias_recuperadas_ambiente,
    file.path(dir_sdm,"ocorrencias_costeiras_recuperadas.csv"))}

if (
  nrow(ocorrencias_excluidas_ambiente) > 0) {
  readr::write_csv(ocorrencias_excluidas_ambiente,
    file.path(dir_sdm,"ocorrencias_excluidas_por_na_ambiental.csv"))}

# 6.9 MANTER SOMENTE OCORRÊNCIAS VÁLIDAS
manter_ocorrencias <- ambiente_completo & !is.na(sp_pts$cell_id)

# Filtrar sp_pts e valores ambientais simultaneamente
sp_pts <- sp_pts[manter_ocorrencias,,drop = FALSE]
valores_ambientais_sp <- valores_ambientais_sp[manter_ocorrencias,,drop = FALSE]

if (
  nrow(sp_pts) !=
  nrow(valores_ambientais_sp)) {
  stop("sp_pts e valores_ambientais_sp ficaram desalinhados após o filtro ambiental.")}

# 6.10 MANTER UMA OCORRÊNCIA POR CÉLULA
# Identificar a primeira ocorrência de cada célula
manter_pixel <- !duplicated(sp_pts$cell_id)

# Salvar duplicatas removidas, caso existam
ocorrencias_duplicadas_pixel <- sp_pts[!manter_pixel,,drop = FALSE]

if (
  nrow(ocorrencias_duplicadas_pixel) > 0) {
  readr::write_csv(ocorrencias_duplicadas_pixel,
    file.path(dir_sdm,"ocorrencias_duplicadas_por_pixel.csv"))}

# Aplicar o mesmo filtro aos dois objetos
sp_pts <- sp_pts[manter_pixel,,drop = FALSE]
valores_ambientais_sp <- valores_ambientais_sp[manter_pixel,,drop = FALSE]

rownames(sp_pts) <- NULL
rownames(valores_ambientais_sp) <- NULL

if (
  nrow(sp_pts) !=
  nrow(valores_ambientais_sp)) {
  stop("sp_pts e valores_ambientais_sp ficaram desalinhados após o filtro por pixel.")}

# 6.11 ATUALIZAR AS COORDENADAS FINAIS
coords_sp <- as.matrix(sp_pts[,c(
      "longitude.gazetteer",
      "latitude.gazetteer"),
    drop = FALSE])

storage.mode(coords_sp) <- "numeric"

# 6.12 DIAGNÓSTICO FINAL
n_occ_pixel <- nrow(sp_pts)

cat(
  "\n============================================================",
  "\nRESULTADO FINAL DA FILTRAGEM",
  "\n============================================================",
  "\nOcorrências originais:",nrow(sp),
  "\nOcorrências completas na célula original:",sum(ambiente_original_completo),
  "\nOcorrências costeiras recuperadas:",sum(sp_pts$recuperada_celula_proxima),
  "\nOcorrências removidas por ambiente ou célula inválida:",sum(!manter_ocorrencias),
  "\nOcorrências duplicadas por pixel removidas:",sum(!manter_pixel),
  "\nOcorrências finais:",n_occ_pixel,
  "\n============================================================\n")

if (n_occ_pixel == 0) {
  stop(paste0(
      "Nenhuma ocorrência permaneceu após a filtragem ambiental ",
      "e a remoção de duplicatas por pixel."))}

# 6.13 DEFINIR O PROTOCOLO DE VALIDAÇÃO
protocolo <- obter_protocolo(n_occ_pixel)
n_occ_pixel <- nrow(sp_pts)
protocolo <- obter_protocolo(n_occ_pixel) #se tiver no automático escolhe entre o LOOCV, 5-fold, part_senv

if (protocolo$thinning_km > 0) {
  
  thin_result <- spThin::thin(
    loc.data = sp_pts,
    lat.col = "latitude.gazetteer",
    long.col = "longitude.gazetteer",
    spec.col = "scientificName_det",
    thin.par = protocolo$thinning_km,
    reps = 100,
    locs.thinned.list.return = TRUE,
    write.files = FALSE,
    write.log.file = FALSE,
    verbose = FALSE
  )
  
  if (
    is.null(thin_result) ||
    length(thin_result) == 0
  ) {
    stop("O thinning não retornou nenhuma repetição.")
  }
  
  # Quantidade de ocorrências mantidas em cada repetição
  tamanhos <- vapply(
    thin_result,
    nrow,
    integer(1)
  )
  
  if (all(tamanhos == 0)) {
    stop("Todas as repetições do thinning ficaram vazias.")
  }
  
  # Escolher a repetição que manteve mais ocorrências
  indice_melhor <- which.max(tamanhos)
  
  sp_thin <- thin_result[[indice_melhor]]
  
  # Criar occ somente depois de escolher a melhor repetição
  occ <- sp_thin |>
    dplyr::transmute(
      scientificName_det = sp_nome,
      x = as.numeric(Longitude),
      y = as.numeric(Latitude))
  
  cat(
    "\nMelhor repetição do thinning:",indice_melhor,
    "\nOcorrências mantidas:",nrow(sp_thin),
    "\nIntervalo entre as 100 repetições:",min(tamanhos),
    "a", max(tamanhos),"\n")
  
} else {
  
  occ <- sp_pts |>
    dplyr::transmute(
      scientificName_det = sp_nome,
      x = longitude.gazetteer,
      y = latitude.gazetteer)
}

occ <- occ |> dplyr::distinct(x,y, .keep_all = TRUE)

occ$pr_ab <- 1
n_occ_thin <- nrow(occ)
n_removidas_pixel <- n_occ_inicial - n_occ_pixel
n_removidas_thinning <- n_occ_pixel - n_occ_thin

cat("\nEspécie:", sp_nome,
  "\nOcorrências iniciais:", n_occ_inicial,
  "\nApós filtro por pixel:", n_occ_pixel,
  "\nRemovidas por pixel:", n_removidas_pixel,
  "\nBuffer da área M:", protocolo$buffer_km, "km",
  "\nApós thinning:", n_occ_thin,
  "\nRemovidas pelo thinning:", n_removidas_thinning,
  "\nThinning utilizado:", protocolo$thinning_km, "km",
  "\nProtocolo:", protocolo$origem,
  "\nDescrição:", protocolo$descricao,
  "\nValidação escolhida:", protocolo$validacao,"\n")

#plot(preditores_continuos[["sand"]],main = "Ocorrências removidas por NA do solo")
#points(ocorrencias_excluidas_ambiente$longitude.gazetteer, ocorrencias_excluidas_ambiente$latitude.gazetteer, pch = 21,bg = "red",col = "white")

# ======================================================================
# 7. ÁREA DE CALIBRAÇÃO M
# ======================================================================
ca <- flexsdm::calib_area(data = occ,
  x = "x", y = "y",
  method = c("buffer", width = protocolo$buffer_km * 1000),
  crs = terra::crs(preditores_continuos))
plot(preditores_continuos[[1]],main = "Área de calibração")
plot(ca,add = TRUE, border = "red",lwd = 2)
points(occ$x, occ$y,pch = 20,col = "blue")

# ======================================================================
# 8. PCA RESTRITA À ÁREA M
# ======================================================================
env_pca <- tryCatch(
  {flexsdm::correct_colinvar(
      env_layer = preditores_continuos,
      method = "pca",
      restric_to_region = ca)},
  error = function(e) {
    stop(paste("Falha ao calcular a PCA dentro da área M.",
        "A análise foi interrompida para evitar",
        "usar uma PCA nacional sem intenção.",
        "Mensagem:", conditionMessage(e)))})

env_pca_modelo <- env_pca$env_layer

if (
  is.null(env_pca_modelo) ||
  terra::nlyr(env_pca_modelo) == 0) {
  stop("A PCA não retornou componentes ambientais.")}

names(env_pca_modelo) <- paste0("PC", seq_len(terra::nlyr(env_pca_modelo)))
pc_names <- names(env_pca_modelo)

# Preditores finais: PCs + distância da costa diretamente
if (usar_distancia_costa) {
  if (!terra::compareGeom(
    env_pca_modelo[[1]],
    distancia_costa,
    stopOnError = FALSE)) {
    stop(paste("A distância da costa não possui",
        "a mesma geometria dos componentes da PCA."))}}

if (usar_distancia_costa) {
  env_modelo <- c(
    env_pca_modelo,
    distancia_costa)
  preditores_modelo <- c(
    pc_names,
    "dist_coast_100km")
} else {
  env_modelo <- env_pca_modelo
  preditores_modelo <- pc_names}

cat("\nPCs utilizadas:", paste(pc_names,collapse = ", "),
    "\nPreditores finais do modelo:",paste(preditores_modelo, collapse = ", "),"\n")

dir_pca <- file.path(dir_sdm,"PCA")
dir.create( dir_pca,recursive = TRUE,showWarnings = FALSE)

#terra::writeRaster(env_modelo, filename = file.path(dir_pca,
#    paste0(nome_sp_arquivo,"_preditores_finais.tif")),overwrite = TRUE)

if (!is.null(env_pca$coefficients)) {readr::write_csv(as.data.frame(env_pca$coefficients),
    file.path(dir_pca,"pca_coefficients.csv"))}
if (!is.null(env_pca$cumulative_variance)) {readr::write_csv(as.data.frame(env_pca$cumulative_variance),
    file.path(dir_pca,"pca_variance.csv"))}

env_M <- terra::mask(terra::crop(preditores_continuos, ca), ca)

teste_mean <- terra::global(
  env_M,
  "mean",
  na.rm = TRUE
)

teste_sd <- terra::global(
  env_M,
  "sd",
  na.rm = TRUE
)

env_pca_means <- as.numeric(teste_mean$mean)
env_pca_sds   <- as.numeric(teste_sd$sd)

names(env_pca_means) <- rownames(teste_mean)
names(env_pca_sds)   <- rownames(teste_sd)

# ======================================================================
# 9. BACKGROUND DENTRO DE M
# ======================================================================
masked_rast <- terra::mask(env_pca_modelo[[1]],ca)
n_celulas_validas <- as.integer(terra::global(!is.na(masked_rast),"sum",na.rm = TRUE)[1, 1])
n_background <- floor(n_celulas_validas * proporcao_background_celulas)

n_background <- min(n_background_max,n_background,n_celulas_validas)
n_background <- max(n_background, min(n_background_min, n_celulas_validas))
fracao_background_real <- n_background / n_celulas_validas

cat(
  "\nCélulas válidas em M:", n_celulas_validas,
  "\nBackground selecionado:", n_background,
  "\nFração real:", round(fracao_background_real, 4),
  "\nPercentual real:", round(fracao_background_real * 100, 2), "%",
  "\n")

bg <- flexsdm::sample_background(
  data = occ,
  x = "x",
  y = "y",
  n = n_background,
  method = "random",
  rlayer = env_pca_modelo[[1]],
  calibarea = ca,
  sp_name = NULL)

bg$pr_ab <- 0

plot(masked_rast, main = paste(sp_nome,"\nÁrea M e Background"))
plot(ca, add = TRUE, border = "blue",lwd = 2)
points(bg$x,bg$y,pch = 16,cex = 0.15, col = rgb(0,0,0,0.25))
points(occ$x,occ$y,pch = 21,bg = "red",col = "white",cex = 0.9)

# ======================================================================
# 10. EXTRAÇÃO AMBIENTAL
# ======================================================================
occ_vars <- flexsdm::sdm_extract(
  data = occ,
  x = "x",
  y = "y",
  env_layer = env_modelo,
  variables = preditores_modelo,
  filter_na = TRUE)

bg_vars <- flexsdm::sdm_extract(
  data = bg,
  x = "x",
  y = "y",
  env_layer = env_modelo,
  variables = preditores_modelo,
  filter_na = TRUE)

occ_vars$pr_ab <- 1
bg_vars$pr_ab <- 0

if (usar_distancia_costa) {
  cat("\nDistância da costa nas presenças, em km:\n")
  print(summary(occ_vars$dist_coast_100km * 100))
  cat("\nDistância da costa no background, em km:\n")
  print(summary(bg_vars$dist_coast_100km * 100))}

occ_bg_vars <- dplyr::bind_rows(occ_vars,bg_vars)
n_pres_modelo <- sum(occ_bg_vars$pr_ab == 1)

if (n_pres_modelo < 5) {stop("Restaram menos de 10 presenças após a extração ambiental.")}

# Reavaliar somente no modo automático.
if (
  modo_protocolo == "automatico" &&
  n_pres_modelo != n_occ_thin) {
  warning( paste(n_occ_thin - n_pres_modelo,
      "ocorrências foram removidas por valores ambientais ausentes.",
      "Revise se o protocolo de validação continua adequado."))}

# ======================================================================
# 11. PARTICIONAMENTO
# ======================================================================
if (protocolo$validacao == "loocv") {
  
  # -------------------------------------------------------
  # LOOCV
  # -------------------------------------------------------
  
  part_raw <- flexsdm::part_random(
    data = occ_bg_vars,
    pr_ab = "pr_ab",
    method = c(
      method = "loocv"))
  
  validacao_usada <- "LOOCV"
  
} else if (protocolo$validacao == "kfold")  {
  
  # K-fold
  
  part_raw <- flexsdm::part_random(
    data = occ_bg_vars,
    pr_ab = "pr_ab",
    method = c(
      method = "kfold",
      folds = protocolo$folds))
  
  validacao_usada <- paste0(
    protocolo$folds,
    "-fold")
  
} else if (protocolo$validacao == "part_senv") {
  
  # Primeira tentativa de part_senv: 3 a 5 grupos, incluindo coordenadas
  
  cat("\nTentando part_senv com 3 a 5 grupos","e coordenadas incluídas...\n")
  
  part_result <- tryCatch(
    {
      flexsdm::part_senv(
        env_layer = env_modelo,
        data = occ_bg_vars,
        x = "x",
        y = "y",
        pr_ab = "pr_ab",
        min_n_groups = 3,
        max_n_groups = 5,
        min_occ = 3,
        prop = 0.2,
        include_coords = TRUE)},
    error = function(e) {
      
      warning(paste("A primeira tentativa de part_senv gerou erro:",
          conditionMessage(e)))
      
      NULL})
  
  # Verificar se a primeira tentativa retornou
  # a estrutura correta.
  part_senv_valido <-
    is.list(part_result) &&
    "part" %in% names(part_result) &&
    !is.null(part_result$part) &&
    is.data.frame(part_result$part) &&
    nrow(part_result$part) > 0 &&
    ".part" %in% names(part_result$part)
  
  # Segunda tentativa de part_senv
  
  if (!part_senv_valido) {
    
    warning(
      paste(
        "A primeira tentativa de part_senv",
        "não encontrou uma partição adequada.",
        "Será feita uma segunda tentativa",
        "sem incluir as coordenadas."))
    
    cat(
      "\nTentando part_senv com 2 a 5 grupos", "e sem coordenadas...\n")
    
    part_result <- tryCatch(
      {
        flexsdm::part_senv(
          env_layer = env_modelo,
          data = occ_bg_vars,
          x = "x",
          y = "y",
          pr_ab = "pr_ab",
          min_n_groups = 2,
          max_n_groups = 5,
          min_occ = 3,
          prop = 0.2,
          include_coords = FALSE)},
      error = function(e) {
        
        warning(
          paste(
            "A segunda tentativa de part_senv gerou erro:",
            conditionMessage(e)))
        
        NULL
      })
    
    part_senv_valido <-
      is.list(part_result) &&
      "part" %in% names(part_result) &&
      !is.null(part_result$part) &&
      is.data.frame(part_result$part) &&
      nrow(part_result$part) > 0 &&
      ".part" %in% names(part_result$part)}
  
  # -------------------------------------------------------
  # Se part_senv falhou nas duas tentativas: usar 5-fold
  # -------------------------------------------------------
  
  if (!part_senv_valido) {
    
    warning(
      paste("O part_senv não encontrou uma partição adequada",
        "nas duas tentativas.",
        "Será utilizada validação 5-fold."))
    
    part_raw <- flexsdm::part_random(
      data = occ_bg_vars,
      pr_ab = "pr_ab",
      method = c(
        method = "kfold",
        folds = 5))
    
    validacao_usada <- "5-fold fallback"
    
  } else {
    
    # -----------------------------------------------------
    # part_senv funcionou
    # -----------------------------------------------------
    
    part_raw <- part_result$part
    validacao_usada <- "part_senv"
    
    #if (
    #  "best_part_info" %in% names(part_result) &&
    #  !is.null(part_result$best_part_info)) {
      
    #  readr::write_csv(
    #    as.data.frame(
    #      part_result$best_part_info),
    #    file.path(dir_sdm,
    #      "part_senv_diagnostico.csv"))}
    }
  
} else {
  
  stop(
    paste("Método de validação não reconhecido:",
      protocolo$validacao))} 

validacao_usada
names(part_raw)
table(part_raw$.part, part_raw$pr_ab)

# Guardar somente a informação da partição.
particoes <- part_raw |> dplyr::select(x,y,pr_ab,.part) |> dplyr::distinct()

# Reunir partições e preditores uma única vez.
occ_partition_vars <- occ_bg_vars |>
  dplyr::left_join(particoes,
    by = c("x", "y", "pr_ab"),
    relationship = "many-to-one")

if (anyNA(occ_partition_vars$.part)) {stop("Algumas linhas ficaram sem partição após o join.")}

if (!all(
  preditores_modelo %in%
  names(occ_partition_vars))) {
  stop("Os componentes da PCA não estão presentes após o particionamento.")}

# Manter somente partições presentes em ambos os grupos.
part_pres <- sort(unique(occ_partition_vars$.part[occ_partition_vars$pr_ab == 1]))
part_bg <- sort(unique(occ_partition_vars$.part[occ_partition_vars$pr_ab == 0]))
part_comuns <- intersect(part_pres, part_bg)
occ_partition_vars <- occ_partition_vars |> filter(.part %in% part_comuns)
background_vars <- occ_partition_vars |> filter(pr_ab == 0)

cat(
  "\nValidação usada:",validacao_usada,"
  \nPartições:",paste(part_comuns, collapse = ", "),"\n")

#write.table(occ_partition_vars,
#  file = file.path(dir_sdm, paste0(nome_sp_arquivo,"_dados_modelagem.txt")),sep = "\t",
#  row.names = FALSE,quote = FALSE)

occ_qgis <- occ |> dplyr::select(x, y)
write.table(occ_qgis,file = file.path(dir_sdm,paste0(nome_sp_arquivo, "_ocorrencias_finais_QGIS.txt")),
  sep = "\t",row.names = FALSE,quote = FALSE)

# ======================================================================
# 12. TUNING DO MAXENT
# ======================================================================
gridtest <- expand.grid(regmult = grid_regmult_padrao, classes = grid_classes_padrao)
print(gridtest)

occ_partition_vars$.part <-match(occ_partition_vars$.part,sort(unique(occ_partition_vars$.part)))
background_vars$.part <-match(background_vars$.part,sort(unique(background_vars$.part)))

max_tun <- flexsdm::tune_max(
  data = occ_partition_vars,
  response = "pr_ab",
  predictors = preditores_modelo,
  predictors_f = NULL,
  partition = ".part",
  background = background_vars,
  grid = gridtest,
  thr = threshold_modelo,
  metric = metrica_selecao,
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = n_cores)

 print(max_tun$performance)

#openxlsx::write.xlsx(max_tun$performance, file = file.path(
#    dir_sdm, paste0(nome_sp_arquivo, "_performance.xlsx")), overwrite = TRUE)

# ======================================================================
# 13. PREDIÇÃO NACIONAL
# ======================================================================
pred_nacional <- flexsdm::sdm_predict(
  models = max_tun,
  pred = env_modelo,
  thr = threshold_modelo,
  con_thr = FALSE,
  predict_area = NULL)

raster_cont_nacional <- pred_nacional$max[["max"]]
raster_bin_nacional <- pred_nacional$max[[threshold_modelo]]
cat("\nPixels adequados:",sum(terra::values(raster_bin_nacional)==1,na.rm=TRUE),"\n")

salvar_raster(raster_cont_nacional, paste0(nome_sp_arquivo, "_Brasil_continuo.tif"))
salvar_raster(raster_bin_nacional, paste0(nome_sp_arquivo,"_Brasil_binario.tif"),datatype = "INT1U")

#Área adequada
area_pixel <- terra::cellSize(raster_bin_nacional,unit="km")
area_total <- global(area_pixel * raster_bin_nacional,"sum",na.rm=TRUE)
cat("\nÁrea prevista:",round(area_total[[1]],0),"km²\n")

# ======================================================================
# 14. PREDIÇÃO RESTRITA AO BUFFER
# ======================================================================
raster_cont_buffer <- terra::mask(terra::crop(raster_cont_nacional, ca),ca)
raster_bin_buffer <- terra::mask(terra::crop(raster_bin_nacional,ca),ca)

salvar_raster(raster_cont_buffer,paste0(nome_sp_arquivo,"_Buffer_continuo.tif"))
salvar_raster(raster_bin_buffer,paste0(nome_sp_arquivo,"_Buffer_binario.tif"), datatype = "INT1U")

# ======================================================================
# 15. MÁSCARA DE VEGETAÇÃO
# ======================================================================
# Valores padrão caso a máscara não seja produzida.
classes_vegetacao_usadas <- NA_character_
area_adequada_vegetacao <- NA_real_
area_total_vegetacao <- NA_real_
percentual_adequado_vegetacao <- NA_real_

if (usar_mascara_vegetacao) {
  
  # Verificar o arquivo de vegetação
  if (!file.exists(arquivo_vegetacao)) {
    
    warning(paste("Shapefile de vegetação não encontrado:",arquivo_vegetacao))
    
  } else {
    
    # Ler e validar o shapefile
    vegetacao <- sf::st_read(arquivo_vegetacao,quiet = TRUE) |> sf::st_make_valid()
    
    if (!coluna_vegetacao %in% names(vegetacao)) {
      
      warning(paste("Coluna de vegetação não encontrada:", coluna_vegetacao))
      
    } else {
      
      # Padronizar a coluna das classes
      vegetacao$classe_modelo <- as.character(vegetacao[[coluna_vegetacao]])
      vegetacao <- vegetacao |> dplyr::mutate(classe_modelo = trimws(classe_modelo)) |>
        dplyr::filter(!is.na(classe_modelo),classe_modelo != "")
      
      if (nrow(vegetacao) == 0) {
        
        warning(paste("Nenhum polígono possui uma classe válida","na coluna:",
            coluna_vegetacao))
        
      } else {
        
        # Preparar vegetação e ocorrências para interseção
        # CRS métrico para análises em escala brasileira.
        crs_intersecao <- 5880
        
        vegetacao_intersecao <- vegetacao |> sf::st_transform(crs_intersecao)
        
        ocorrencias_sf <- sf::st_as_sf(occ,
          coords = c("x","y"),
          crs = 4326,
          remove = FALSE) |>
          sf::st_transform(crs_intersecao)
        
        # Identificar a vegetação de cada ocorrência
        
        ocorrencias_vegetacao <- sf::st_join(
          ocorrencias_sf,
          vegetacao_intersecao |> dplyr::select(classe_modelo),
          join = sf::st_intersects,
          left = TRUE)
    
        # Frequência das classes nas ocorrências
        
        frequencia_vegetacao <- ocorrencias_vegetacao |>
          sf::st_drop_geometry() |> dplyr::count(classe_modelo,sort = TRUE,
            name = "n_ocorrencias")
        
        #readr::write_csv(
        #  frequencia_vegetacao,
        #  file.path(dir_sdm,"frequencia_vegetacao_ocorrencias.csv"))
        
        cat("\nFrequência das classes de vegetação:\n")
        
        print(frequencia_vegetacao)
        
        # Escolher classes automáticas ou manuais
        if (is.null(classes_vegetacao_manuais)) {
          
          classes_vegetacao_usadas <-
            frequencia_vegetacao |>
            dplyr::filter(!is.na(classe_modelo)) |>
            dplyr::pull(classe_modelo) |>
            unique() |>
            sort()
          
          if (
            is.null(max_classes_vegetacao_automaticas) ||
            length(max_classes_vegetacao_automaticas) == 0 ||
            is.na(max_classes_vegetacao_automaticas)) {
            
            max_classes_vegetacao_automaticas <- 3
            
            warning(paste("max_classes_vegetacao_automaticas",
                "não estava definido.", "Foi utilizado o valor padrão de 3."))}
          
          if (
            length(classes_vegetacao_usadas) >
            max_classes_vegetacao_automaticas) {
            
            warning(paste("Foram detectadas",length(classes_vegetacao_usadas),
                "classes de vegetação.", "Revise o arquivo",
                "frequencia_vegetacao_ocorrencias.csv",
                "e considere definir","classes_vegetacao_manuais."))}
          
        } else {
          
          classes_vegetacao_usadas <- as.character(classes_vegetacao_manuais)
          classes_vegetacao_usadas <- trimws(classes_vegetacao_usadas)
          
          classes_vegetacao_usadas <- classes_vegetacao_usadas[
              !is.na(classes_vegetacao_usadas) &
                classes_vegetacao_usadas != ""]
          
          classes_vegetacao_usadas <- sort(
            unique(classes_vegetacao_usadas))}
        
        # Mostrar as classes selecionadas
        cat("\nClasses de vegetação selecionadas:\n")
        
        if (length(classes_vegetacao_usadas) > 0) {
          
          cat("- ",paste(classes_vegetacao_usadas,
              collapse = "\n- "),"\n")
          
        } else {
          
          cat("Nenhuma classe selecionada.\n")}
        
        # Conferir se classes manuais existem no shapefile
        classes_disponiveis <- sort(
          unique(vegetacao$classe_modelo))
        
        classes_nao_encontradas <- setdiff(
          classes_vegetacao_usadas, classes_disponiveis)
        
        if (length(classes_nao_encontradas) > 0) {
          
          warning(paste(
              "As seguintes classes não foram encontradas",
              "no shapefile:", paste(classes_nao_encontradas,collapse = ", ")))}
        
        # Manter apenas classes realmente existentes.
        classes_vegetacao_usadas <- intersect(
          classes_vegetacao_usadas, classes_disponiveis)
        
        # Selecionar os polígonos correspondentes
        vegetacao_selecionada <- vegetacao |>
          dplyr::filter(classe_modelo %in% classes_vegetacao_usadas)
        
        if (length(classes_vegetacao_usadas) == 0) {
          
          warning(paste("A máscara de vegetação não será produzida",
              "porque nenhuma classe válida foi selecionada."))
          
        } else if (nrow(vegetacao_selecionada) == 0) {
          
          warning(paste(
              "Nenhum polígono corresponde às classes",
              "de vegetação selecionadas."))
          
        } else {
          
          # Salvar a lista das classes utilizadas
          #readr::write_csv( data.frame(classe_modelo =
          #      classes_vegetacao_usadas),
          #  file.path(dir_sdm,"classes_vegetacao_usadas.csv"))
          
          # Converter a vegetação para SpatVector
          vegetacao_vect <- terra::vect(
            sf::st_transform(vegetacao_selecionada,
              terra::crs(raster_cont_nacional)))
          
          # Criar a máscara raster
          mascara_vegetacao <- terra::rasterize(
            vegetacao_vect,
            raster_cont_nacional,
            field = 1,
            background = NA,
            touches = TRUE)
          
          # Aplicar a máscara às predições
          raster_cont_vegetacao <- terra::mask(
            raster_cont_nacional, mascara_vegetacao)
          
          raster_bin_vegetacao <- terra::mask(
            raster_bin_nacional, mascara_vegetacao)
          
          # Verificar se a máscara possui pixels válidos
          n_pixels_vegetacao <- terra::global(
            !is.na(mascara_vegetacao),
            "sum", na.rm = TRUE)[1, 1]
          
          if (
            is.na(n_pixels_vegetacao) ||
            n_pixels_vegetacao == 0) {
            
            warning(paste(
                "A máscara de vegetação foi criada,",
                "mas não possui pixels válidos."))
            
          } else {
            
            # Salvar os rasters
            salvar_raster(raster_cont_vegetacao,
              paste0(nome_sp_arquivo, "_Vegetacao_continuo.tif"))
            
            salvar_raster(raster_bin_vegetacao,
              paste0(nome_sp_arquivo,"_Vegetacao_binario.tif"),
              datatype = "INT1U")
            
            salvar_raster(mascara_vegetacao,
              paste0(nome_sp_arquivo,"_mascara_vegetacao.tif"),datatype = "INT1U")
            
            # Calcular áreas em km²
             area_pixel_vegetacao <- terra::cellSize(
              raster_bin_vegetacao, unit = "km")
            
            area_adequada_vegetacao <- terra::global(
              area_pixel_vegetacao *
                raster_bin_vegetacao,"sum", na.rm = TRUE)[1, 1]
            
            area_total_vegetacao <- terra::global(
              terra::ifel(!is.na(raster_bin_vegetacao),
                area_pixel_vegetacao,NA),
              "sum",na.rm = TRUE)[1, 1]
            
            percentual_adequado_vegetacao <-
              area_adequada_vegetacao /
              area_total_vegetacao * 100
            
            cat("\nMáscara de vegetação concluída.",
              "\nNúmero de classes usadas:", length(classes_vegetacao_usadas),
              "\nPixels válidos na máscara:", n_pixels_vegetacao,
              "\nÁrea total das classes selecionadas:",round(
                area_total_vegetacao, 2),
              "km²", "\nÁrea adequada dentro das classes:",
              round(area_adequada_vegetacao,2),
              "km²","\nPercentual adequado:",
              round(percentual_adequado_vegetacao,2),"%\n")}}}}}}

# ======================================================================
# 16. EXTRAPOLAÇÃO AMBIENTAL
# ======================================================================
if (calcular_extrapolacao) {
  extrapolacao <- flexsdm::extra_eval(
    training_data = occ_partition_vars,
    pr_ab = "pr_ab",
    projection_data = env_modelo,
    metric = "mahalanobis",
    univar_comb = FALSE,
    aggreg_factor = agregacao_extrapolacao)
  
  if (truncar_extrapolacao) {
    
    raster_cont_truncado <- flexsdm::extra_truncate(
      suit = raster_cont_nacional,nextra = extrapolacao,threshold = limiar_extrapolacao,
      trunc_value = NA)
    
  # salvar_raster(raster_cont_truncado,
  #  paste0(nome_sp_arquivo,"_Brasil_continuo_sem_extrapolacao.tif"))
    }}

# ======================================================================
# 17.1 PREENCHIMENTO APENAS PARA FIGURAS
# ======================================================================
#contínuo
preencher_buracos_figura <- function(
    r, mascara, n_iter = 15, janela = 3){r_plot <- r
  for(i in seq_len(n_iter)){
    media <- terra::focal(r_plot,w = matrix(1, janela, janela),fun = mean, na.rm = TRUE, na.policy = "only")
    r_plot <- terra::cover(r_plot, media)}
  r_plot <- terra::mask(r_plot, mascara)
  r_plot}

raster_cont_figura <- preencher_buracos_figura(raster_cont_nacional, brasil_vect)

threshold_final <- as.numeric(max_tun$performance$thr_value[1])
#binário
preencher_buracos_binario_figura <- function(
    r, mascara, n_iter = 15,janela = 3){r_plot <- r
  moda_fun <- function(x, ...){x <- x[!is.na(x)]
    if(length(x) == 0) return(NA)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]}
  for(i in seq_len(n_iter)){
    moda <- terra::focal(
      r_plot, w = matrix(1, janela, janela),
      fun = moda_fun,
      na.policy = "only")
    r_plot <- terra::cover(r_plot, moda)}
  r_plot <- terra::mask(r_plot, mascara)
  
  return(r_plot)
}

raster_bin_figura <- terra::ifel(raster_cont_figura >= threshold_final,1,0)

#salvar os rasters
arquivo_cont <- file.path(dir_sdm,paste0(nome_sp_arquivo, "_Brasil_continuo_figura.tif"))
arquivo_bin <- file.path(dir_sdm,paste0(nome_sp_arquivo, "_Brasil_binario_figura.tif"))

terra::writeRaster(raster_cont_figura, arquivo_cont, overwrite = TRUE)
terra::writeRaster(raster_bin_figura, arquivo_bin, overwrite = TRUE)

# ======================================================================
# 17.2 MAPAS DE PUBLICAÇÃO
# ======================================================================
# CORREÇÃO POSTERIORI - BMCP
buffer_bmcp_km <- 300
buffer_bmcp_m  <- buffer_bmcp_km * 1000

# Preparar dados para o BMCP
occ_bg_bmcp <- occ_bg_vars |>
  dplyr::mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    pr_ab = as.numeric(pr_ab)) |>
  dplyr::filter(
    is.finite(x),
    is.finite(y),
    pr_ab %in% c(0, 1))

# Executar correção posteriori
resultado_bmcp <- flexsdm::msdm_posteriori(
  records = occ_bg_bmcp,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  cont_suit = raster_cont_figura,
  method = "bmcp",
  thr = threshold_final,
  buffer = buffer_bmcp_m,
  con_thr = TRUE,
  crs = terra::crs(raster_cont_figura))

if (!"msdm_cont" %in% names(resultado_bmcp)) {stop("A camada 'msdm_cont' não foi retornada pelo BMCP.")}
raster_bmcp <- resultado_bmcp[["msdm_cont"]]
names(raster_bmcp) <- "Adequabilidade_BMCP"

# ======================================================================
# FIGURA DE PUBLICAÇÃO
# ======================================================================
criar_figura_publicacao(
  raster = raster_bmcp,
  ocorrencias = occ,
  arquivo = paste0(nome_sp_arquivo,"_BMCP_", buffer_bmcp_km, "km_publicacao.png"),
  titulo = sp_nome,
  binario = FALSE,
  zoom = FALSE,
  suavizar = suavizar_apenas_figura,
  mostrar_pontos = TRUE)

terra::plot(raster_cont_nacional)
salvar_raster(raster_cont_obr,paste0(nome_sp_arquivo,"_Brasil_continuo_OBR.tif"))

criar_figura_publicacao(raster_cont_figura,
  occ, paste0(nome_sp_arquivo, "_Brasil_continuo_publicacao.png"),
  paste(sp_nome),
  binario = FALSE, 
  zoom = FALSE,
  suavizar = suavizar_apenas_figura)

criar_figura_publicacao(raster_bin_figura,
  occ, paste0(nome_sp_arquivo,"_Brasil_binario_publicacao.png"),
  paste(sp_nome),
  binario = TRUE,
  zoom = FALSE,
  suavizar = FALSE)

p_cont <- criar_figura_publicacao(
  raster_cont_figura, occ, "teste.png", sp_nome,
  binario = FALSE, zoom = FALSE,
  suavizar = suavizar_apenas_figura)
p_cont

p_bin <- criar_figura_publicacao(
  raster_bin_figura, occ, "teste1.png", sp_nome,
  binario = TRUE, zoom = FALSE,
  suavizar = FALSE)
p_bin

fig_sdm <- p_cont + p_bin +
  patchwork::plot_annotation(
    title = sp_nome,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "italic",
        size = 16,
        hjust = 0.5)))
fig_sdm

ggsave(
  file.path(dir_sdm,paste0(nome_sp_arquivo, "_SDM_publicacao.png")),
  fig_sdm,
  width = 14,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white")

# ======================================================================
# ÁREA DE UC DENTRO DA ÁREA ADEQUADA
# ======================================================================
ucs <- sf::st_read(arquivo_uc, quiet = TRUE)
sdm_bin_sf <- terra::as.polygons(raster_bin_nacional, dissolve = TRUE, values = TRUE) |> sf::st_as_sf()
sdm_bin_sf <- sdm_bin_sf |> dplyr::filter(equal_sens_spec == 1)
sdm_bin_sf <- sf::st_transform(sdm_bin_sf, 5880)
ucs <- sf::st_transform(ucs, 5880)
ucs <- sf::st_make_valid(ucs)
uc_sdm <- sf::st_intersection(ucs, sdm_bin_sf)
uc_sdm_union <- sf::st_union(uc_sdm)
area_adequada_nacional <- terra::global(terra::cellSize(raster_bin_nacional, unit = "km") * raster_bin_nacional, "sum", na.rm = TRUE)[1, 1]
area_uc_km2 <- as.numeric(sf::st_area(uc_sdm_union)) / 1e6
percentual_uc <- (area_uc_km2 / area_adequada_nacional) * 100
cat("\nÁrea protegida dentro do SDM:", round(area_uc_km2, 2), "km²", "\nPercentual protegido:", round(percentual_uc, 2), "%\n")

# ======================================================================
# 18. ÁREAS E RESUMO METODOLÓGICO
# ======================================================================
threshold_final
model_eval <- flexsdm::sdm_summarize(models = list(max_tun))
area_M_km2 <- terra::expanse(ca, unit = "km")

resumo_metodologico <- data.frame(species = sp_nome,
  
  n_occ_initial = n_occ_inicial,
  n_occ_after_pixel_filter = n_occ_pixel,
  n_occ_after_thinning = n_occ_thin,
  n_occ_model = n_pres_modelo,
  
  validation = validacao_usada,
  thinning_km = protocolo$thinning_km,
  buffer_km = protocolo$buffer_km,
  
  threshold_method = threshold_modelo,
  threshold_value = threshold_final,
  
  auc = max(model_eval$AUC_mean, na.rm = TRUE),
  tss = max(model_eval$TSS_mean, na.rm = TRUE),
  #boyce = max(model_eval$Boyce_mean, na.rm = TRUE),
  
  area_M_km2 = area_M_km2,
  
  suitable_area_Brazil_km2 = area_adequada_nacional,
  suitable_area_vegetation_km2 = area_adequada_vegetacao,
  
  protected_area_km2 = area_uc_km2,
  protected_area_percent = area_uc_km2 / area_adequada_nacional * 100
)

readr::write_csv(resumo_metodologico,
  file.path(dir_sdm,paste0(nome_sp_arquivo, "_resumo_modelo.csv")))

cat("\n====================================================",
  "\nSDM robusto finalizado",
  "\nEspécie:",sp_nome,
  "\nProtocolo:",protocolo$origem,
  "\nValidação:",validacao_usada,
  "\nPresenças finais:",n_pres_modelo,
  "\nÁrea adequada nacional (km²):",round(area_adequada_nacional,2),
  "\nResultados:",dir_sdm,
  "\n====================================================\n")
