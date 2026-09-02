# ======================================================================
# SDM FUTURO COM flexsdm — UM CENÁRIO POR VEZ — BAIXO USO DE RAM
# ======================================================================
library(flexsdm)
library(terra)
library(dplyr)
library(readr)
library(openxlsx)
set.seed(123)

# ======================================================================
# 1. ESCOLHA DO CENÁRIO
# ======================================================================
cenario_escolhido <- "2021_2040_SSP370"

# Opções:
# "2021_2040_SSP245"
# "2021_2040_SSP370"
# "2081_2100_SSP245"
# "2081_2100_SSP370"

# ======================================================================
# 2. DIRETÓRIO DOS CENÁRIOS
# ======================================================================

dir_clima_futuro <- "C:/Users/vinil/OneDrive/Área de Trabalho/Qgis/Variaveis climaticas futuras"

# ======================================================================
# 3. TABELA DOS CENÁRIOS
# ======================================================================

cenarios <- data.frame(
  id = c("2021_2040_SSP245","2021_2040_SSP370","2081_2100_SSP245","2081_2100_SSP370"),
  periodo = c("2021-2040","2021-2040","2081-2100","2081-2100"),
  ssp = c("SSP2-4.5","SSP3-7.0","SSP2-4.5","SSP3-7.0"),
  pasta = c("2021_2040_SSP245","2021_2040_SSP370","2081_2100_SSP245","2081_2100_SSP370"),
  stringsAsFactors = FALSE
)

if (!cenario_escolhido %in% cenarios$id) stop("Cenário inválido: ", cenario_escolhido)

cenario <- cenarios[cenarios$id == cenario_escolhido, , drop = FALSE]

# ======================================================================
# 4. VERIFICAR OBJETOS DO MODELO PRESENTE
# ======================================================================

objetos_necessarios <- c("max_tun","env_pca","preditores_continuos","env_modelo","preditores_modelo","ca","sp_nome","dir_sdm","threshold_modelo")

objetos_faltantes <- objetos_necessarios[!vapply(objetos_necessarios, exists, logical(1), inherits = TRUE)]

if (length(objetos_faltantes) > 0) {
  stop("Objetos ausentes no ambiente do R: ", paste(objetos_faltantes, collapse = ", "), "\nExecute primeiro o script 'SDM - flexsdm.R'.")
}

if (!exists("brasil_vect")) stop("O objeto 'brasil_vect' não existe.")

if (!dir.exists(dir_clima_futuro)) stop("Diretório de cenários futuros não encontrado: ", dir_clima_futuro)

nome_sp_arquivo <- gsub("[^A-Za-z0-9_]+", "_", sp_nome)

# ======================================================================
# 5. DIRETÓRIOS
# ======================================================================

dir_futuro <- file.path(dir_sdm,"Cenarios_futuros")
dir_proj_pca <- file.path(dir_futuro,"PCA_projetada")
dir_preditores_futuro <- file.path(dir_futuro,"Preditores_completos")
dir_temp <- file.path(dir_futuro,"TEMP_RASTERS")

dir.create(dir_futuro, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_proj_pca, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_preditores_futuro, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

dir_entrada <- file.path(dir_clima_futuro,cenario$pasta)
dir_preditores <- file.path(dir_preditores_futuro,cenario$id)
dir_pca <- file.path(dir_proj_pca,cenario$id)

dir.create(dir_preditores, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_pca, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(dir_entrada)) stop("Pasta do cenário não encontrada: ", dir_entrada)

# ======================================================================
# 6. CONFIGURAÇÕES
# ======================================================================
restric_pca_proj_futuro <- TRUE
calcular_extrapolacao_futura <- TRUE

aplicar_mascara_vegetacao_futuro <- if (exists("usar_mascara_vegetacao")) {
  isTRUE(usar_mascara_vegetacao)
} else {
  FALSE
}

salvar_pcs_projetados <- TRUE
salvar_continuos <- TRUE
salvar_binarios <- TRUE

# ======================================================================
# 7. CONFIGURAÇÕES DO TERRA PARA REDUZIR USO DE RAM
# ======================================================================
terra::terraOptions(tempdir = dir_temp, memfrac = 0.30, progress = 1)
gc()

cat(
  "\n==========================================================",
  "\nPROJEÇÃO FUTURA — BAIXO USO DE RAM",
  "\n==========================================================",
  "\nEspécie: ", sp_nome,
  "\nCenário: ", cenario$id,
  "\nPeríodo: ", cenario$periodo,
  "\nSSP: ", cenario$ssp,
  "\nMáscara de vegetação: ", aplicar_mascara_vegetacao_futuro,
  "\nExtrapolação: ", calcular_extrapolacao_futura,
  "\nterra memfrac: 0.30\n",
  sep = "")


# ======================================================================
# 8. FUNÇÕES
# ======================================================================
ler_bioclimaticas_futuras <- function(diretorio) {
  arquivos <- list.files(diretorio, pattern="\\.(tif|tiff)$", full.names=TRUE, ignore.case=TRUE)
  if (length(arquivos)==0) stop("Nenhum TIFF encontrado em: ",diretorio)
  if (length(arquivos)>1) stop("Foram encontrados vários TIFFs em ",diretorio,". O script espera UM TIFF contendo as 19 BIOCLIM.")
  bio <- terra::rast(arquivos[1])
  if (terra::nlyr(bio)!=19) stop("O TIFF futuro precisa conter exatamente 19 BIOCLIM. Encontradas: ",terra::nlyr(bio))
  names(bio) <- paste0("bio",1:19)
  bio
}

salvar_futuro <- function(r,arquivo,datatype=NULL) {
  wopt <- list(gdal=c("COMPRESS=DEFLATE","BIGTIFF=YES"))
  if (!is.null(datatype)) wopt$datatype <- datatype
  terra::writeRaster(r,arquivo,overwrite=TRUE,wopt=wopt)
}

calcular_area_binaria <- function(r) {
  if (is.null(r)) return(NA_real_)
  area_pixel <- terra::cellSize(r,unit="km")
  resultado <- terra::global(area_pixel*r,"sum",na.rm=TRUE)
  rm(area_pixel)
  gc()
  as.numeric(resultado[1,1])
}

harmonizar_para_futuro <- function(r,referencia,nome,method="bilinear",arquivo=NULL) {
  names(r) <- nome
  if (!terra::same.crs(r,referencia)) r <- terra::project(r,referencia,method=method,filename=arquivo,overwrite=TRUE)
  if (!terra::compareGeom(r,referencia,stopOnError=FALSE)) r <- terra::resample(r,referencia,method=method,filename=arquivo,overwrite=TRUE)
  if (!terra::compareGeom(r,referencia,stopOnError=FALSE)) stop("Não foi possível harmonizar: ",nome)
  r
}

# ======================================================================
# 9. IDENTIFICAR PCA PRESENTE
# ======================================================================
if (!all(preditores_modelo %in% names(env_modelo))) {
  stop("Preditores do modelo ausentes em env_modelo. Esperados: ",paste(preditores_modelo,collapse=", "))
}

nomes_bio <- paste0("bio",1:19)
preditores_estaticos <- setdiff(names(preditores_continuos),nomes_bio)
nomes_pca <- as.character(env_pca$coefficients$variable)
means_presente <- env_pca_means
stds_presente  <- env_pca_sds

if (!exists("means_presente") || !exists("stds_presente")) {
  stop("Os objetos 'means_presente' e 'stds_presente' precisam existir.")
}

if (!identical(names(means_presente),nomes_pca) || !identical(names(stds_presente),nomes_pca)) {
  stop("A ordem de means_presente/stds_presente não coincide com a PCA presente.")
}

coef_pca <- as.matrix(env_pca$coefficients[,-1,drop=FALSE])
rownames(coef_pca) <- nomes_pca

pcs_modelo <- preditores_modelo[grepl("^PC[0-9]+$",preditores_modelo)]

if (length(pcs_modelo)==0) stop("Nenhum PC foi encontrado em preditores_modelo.")

pc_indices <- match(pcs_modelo,colnames(coef_pca))

if (anyNA(pc_indices)) {
  stop("PCs solicitados pelo modelo não existem na PCA: ",paste(pcs_modelo[is.na(pc_indices)],collapse=", "))
}

# ======================================================================
# 10. LER BIOCLIM E RECORTAR DIRETAMENTE PARA DISCO
# ======================================================================
cat("\n[1/6] Lendo BIOCLIM...\n")

bio_original <- ler_bioclimaticas_futuras(dir_entrada)

arquivo_bio_cropped <- file.path(dir_temp,paste0(cenario$id,"_bio_brasil.tif"))

cat("[1/6] Recortando BIOCLIM para o Brasil...\n")

bio_crop <- terra::crop(
  bio_original,
  brasil_vect,
  filename=arquivo_bio_cropped,
  overwrite=TRUE,
  wopt=list(datatype="FLT4S",gdal=c("COMPRESS=DEFLATE","BIGTIFF=YES"))
)

rm(bio_original)
gc()

cat("[1/6] Aplicando máscara do Brasil...\n")

arquivo_bio_mask <- file.path(dir_temp,paste0(cenario$id,"_bio_brasil_mask.tif"))

bio_fut <- terra::mask(
  bio_crop,
  brasil_vect,
  filename=arquivo_bio_mask,
  overwrite=TRUE,
  wopt=list(datatype="FLT4S",gdal=c("COMPRESS=DEFLATE","BIGTIFF=YES"))
)

rm(bio_crop)
gc()

# ======================================================================
# 11. PREPARAR PREDITORES FUTUROS NO DISCO
# ======================================================================
cat("\n[2/6] Preparando preditores futuros...\n")

nomes_esperados <- names(preditores_continuos)
arquivos_preditores <- character(length(nomes_esperados))
names(arquivos_preditores) <- nomes_esperados

# BIOCLIM
for (j in seq_len(terra::nlyr(bio_fut))) {
  nome_var <- names(bio_fut)[j]
  arquivo_var <- file.path(dir_preditores,paste0(nome_var,".tif"))
  salvar_futuro(bio_fut[[j]],arquivo_var)
  arquivos_preditores[nome_var] <- arquivo_var
  cat("  Salvo: ",nome_var,"\n",sep="")
}

rm(bio_fut)
gc()

# Variáveis estáticas
if (length(preditores_estaticos)>0) {
  for (nome_var in preditores_estaticos) {
    cat("  Preparando: ",nome_var,"\n",sep="")
    arquivo_var <- file.path(dir_preditores,paste0(nome_var,".tif"))
    r_static <- preditores_continuos[[nome_var]]
    r_ref <- terra::rast(arquivos_preditores["bio1"])
    r_static <- harmonizar_para_futuro(r_static,r_ref,nome_var,"bilinear")
    salvar_futuro(r_static,arquivo_var)
    arquivos_preditores[nome_var] <- arquivo_var
    rm(r_static,r_ref)
    gc()
  }
}

faltantes <- nomes_esperados[!file.exists(arquivos_preditores)]

if (length(faltantes)>0) stop("Preditores futuros ausentes: ",paste(faltantes,collapse=", "))

# ======================================================================
# 12. PROJEÇÃO PCA DIRETAMENTE DOS TIFFs
# ======================================================================
cat("\n[3/6] Projetando PCA em blocos...\n")

env_fut_pca <- terra::rast(arquivos_preditores[nomes_pca])
names(env_fut_pca) <- nomes_pca

arquivo_pca_saida <- file.path(dir_pca,paste0(nome_sp_arquivo,"_",cenario$id,"_PCs.tif"))

vars <- rownames(coef_pca)
loadings <- coef_pca[,pc_indices,drop=FALSE]
means_vec <- as.numeric(means_presente[vars])
sds_vec <- as.numeric(stds_presente[vars])

if (any(!is.finite(means_vec)) || any(!is.finite(sds_vec)) || any(sds_vec==0)) {
  stop("Médias/SDs inválidos na PCA presente.")
}

pca_fun <- function(v) {
  if (all(is.na(v)) || anyNA(v)) return(rep(NA_real_,ncol(loadings)))
  z <- (v-means_vec)/sds_vec
  as.numeric(crossprod(z,loadings))
}

terra::app(
  env_fut_pca,
  fun=pca_fun,
  cores=1,
  filename=arquivo_pca_saida,
  overwrite=TRUE,
  wopt=list(
    datatype="FLT4S",
    gdal=c("COMPRESS=DEFLATE","BIGTIFF=YES")
  )
)

pca_fut <- terra::rast(arquivo_pca_saida)
names(pca_fut) <- pcs_modelo

rm(env_fut_pca)
gc()

cat("[3/6] PCA concluída.\n")

# ======================================================================
# 13. LIBERAR OS PREDITORES ANTES DO MAXENT
# ======================================================================

rm(pca_fun,loadings,means_vec,sds_vec)
gc()

# ======================================================================
# 14. AMBIENTE DO MAXENT
# ======================================================================

cat("\n[4/6] Preparando ambiente do MaxEnt...\n")

if ("dist_coast_100km" %in% preditores_modelo) {
  if (!exists("distancia_costa")) stop("O modelo usa dist_coast_100km, mas distancia_costa não existe.")
  arquivo_dist <- file.path(dir_temp,paste0(cenario$id,"_dist_coast_100km.tif"))
  distancia_fut <- harmonizar_para_futuro(distancia_costa,pca_fut[[1]],"dist_coast_100km","bilinear",arquivo_dist)
  salvar_futuro(distancia_fut,arquivo_dist)
  env_modelo_fut <- terra::rast(c(arquivo_pca_saida,arquivo_dist))
  names(env_modelo_fut) <- c(pcs_modelo,"dist_coast_100km")
  rm(distancia_fut)
} else {
  env_modelo_fut <- pca_fut
}

faltantes_modelo <- setdiff(preditores_modelo,names(env_modelo_fut))

if (length(faltantes_modelo)>0) stop("Preditores necessários pelo MaxEnt ausentes: ",paste(faltantes_modelo,collapse=", "))

env_modelo_fut <- env_modelo_fut[[match(preditores_modelo,names(env_modelo_fut))]]
names(env_modelo_fut) <- preditores_modelo

gc()

# ======================================================================
# 15. MAXENT
# ======================================================================

cat("\n[5/6] Executando sdm_predict()...\n")

pred_fut <- flexsdm::sdm_predict(
  models=max_tun,
  pred=env_modelo_fut,
  thr=threshold_modelo,
  con_thr=FALSE,
  predict_area=NULL,
  clamp=TRUE,
  pred_type="cloglog"
)

cat("[5/6] sdm_predict() concluído.\n")

if (!"max" %in% names(pred_fut)) stop("sdm_predict() não retornou 'max'.")

raster_cont_fut <- pred_fut$max[["max"]]
raster_bin_fut <- pred_fut$max[[threshold_modelo]]

if (is.null(raster_cont_fut) || is.null(raster_bin_fut)) stop("Predição contínua/binária ausente.")

names(raster_cont_fut) <- "suitability"
names(raster_bin_fut) <- "suitable"

# ======================================================================
# 16. ÁREA NACIONAL SEM values()
# ======================================================================

cat("\nCalculando área adequada...\n")

area_nacional <- calcular_area_binaria(raster_bin_fut)

# ======================================================================
# 17. BUFFER
# ======================================================================

arquivo_buffer_cont <- file.path(dir_temp,paste0(cenario$id,"_buffer_cont.tif"))
arquivo_buffer_bin <- file.path(dir_temp,paste0(cenario$id,"_buffer_bin.tif"))

raster_cont_buffer_fut <- terra::mask(
  terra::crop(raster_cont_fut,ca),
  ca,
  filename=arquivo_buffer_cont,
  overwrite=TRUE,
  wopt=list(datatype="FLT4S",gdal=c("COMPRESS=DEFLATE","BIGTIFF=YES"))
)

raster_bin_buffer_fut <- terra::mask(
  terra::crop(raster_bin_fut,ca),
  ca,
  filename=arquivo_buffer_bin,
  overwrite=TRUE,
  wopt=list(datatype="INT1U",gdal=c("COMPRESS=DEFLATE","BIGTIFF=YES"))
)

area_buffer <- calcular_area_binaria(raster_bin_buffer_fut)

gc()

# ======================================================================
# 18. MÁSCARA DE VEGETAÇÃO
# ======================================================================
raster_cont_vegetacao_fut <- NULL
raster_bin_vegetacao_fut <- NULL
area_vegetacao <- NA_real_

if (aplicar_mascara_vegetacao_futuro && exists("mascara_vegetacao")) {
  cat("\nAplicando máscara de vegetação...\n")
  mascara_fut <- harmonizar_para_futuro(mascara_vegetacao,raster_cont_fut,"mask_vegetacao","near")
  raster_cont_vegetacao_fut <- terra::mask(raster_cont_fut,mascara_fut)
  raster_bin_vegetacao_fut <- terra::mask(raster_bin_fut,mascara_fut)
  area_vegetacao <- calcular_area_binaria(raster_bin_vegetacao_fut)
  rm(mascara_fut)
  gc()
}

# ======================================================================
# 19. EXTRAPOLAÇÃO
# ======================================================================
extrapolacao_fut <- NULL

if (calcular_extrapolacao_futura && exists("occ_partition_vars") && exists("agregacao_extrapolacao")) {
  cat("\nCalculando extrapolação...\n")
  extrapolacao_fut <- tryCatch(
    flexsdm::extra_eval(
      training_data=occ_partition_vars,
      pr_ab="pr_ab",
      projection_data=env_modelo_fut,
      metric="mahalanobis",
      univar_comb=FALSE,
      aggreg_factor=agregacao_extrapolacao
    ),
    error=function(e) {
      warning("extra_eval não foi calculado: ",conditionMessage(e))
      NULL
    }
  )
}

# ======================================================================
# 20. SALVAR RESULTADOS
# ======================================================================

cat("\n[6/6] Salvando resultados...\n")

if (salvar_pcs_projetados) {salvar_futuro(pca_fut,arquivo_pca_saida)}

if (salvar_continuos) {
  salvar_futuro(raster_cont_fut,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_Brasil_continuo.tif")))
  #salvar_futuro(raster_cont_buffer_fut,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_Buffer_continuo.tif")))
}

if (salvar_binarios) {
  salvar_futuro(raster_bin_fut,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_Brasil_binario.tif")),"INT1U")
  #salvar_futuro(raster_bin_buffer_fut,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_Buffer_binario.tif")),"INT1U")
}

if (!is.null(raster_cont_vegetacao_fut) && salvar_continuos) {
  salvar_futuro(raster_cont_vegetacao_fut,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_Vegetacao_continuo.tif")))
}

if (!is.null(raster_bin_vegetacao_fut) && salvar_binarios) {
  salvar_futuro(raster_bin_vegetacao_fut,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_Vegetacao_binario.tif")),"INT1U")
}

# ======================================================================
# 21. THRESHOLD
# ======================================================================

threshold_value <- if (!is.null(max_tun$performance) && "thr_value" %in% names(max_tun$performance)) as.numeric(max_tun$performance$thr_value[1]) else NA_real_

# IMPORTANTE: não usar terra::values(), pois isso pode carregar todo o raster na RAM.
n_suitable_pixels <- as.numeric(terra::global(raster_bin_fut,"sum",na.rm=TRUE)[1,1])

# ======================================================================
# 22. RESUMO
# ======================================================================

resumo_futuro <- data.frame(
  species=sp_nome,
  scenario_id=cenario$id,
  period=cenario$periodo,
  SSP=cenario$ssp,
  threshold=threshold_modelo,
  threshold_value=threshold_value,
  suitable_area_Brazil_km2=area_nacional,
  #suitable_area_buffer_km2=area_buffer,
  suitable_area_current_vegetation_km2=area_vegetacao,
  n_suitable_pixels=n_suitable_pixels,
  stringsAsFactors=FALSE
)

# ======================================================================
# 23. ÁREAS DO PRESENTE
# ======================================================================

if (!exists("raster_bin_nacional") || !exists("raster_bin_buffer")) {
  stop("Os rasters binários do presente não existem no ambiente do R.")
}

area_presente_nacional <- calcular_area_binaria(raster_bin_nacional)
area_presente_buffer <- calcular_area_binaria(raster_bin_buffer)
area_presente_vegetacao <- if (exists("raster_bin_vegetacao")) calcular_area_binaria(raster_bin_vegetacao) else NA_real_

sum(terra::values(raster_bin_nacional) == 1, na.rm = TRUE)
n_suitable_pixels_presente <- as.numeric(terra::global(raster_bin_nacional, "sum", na.rm = TRUE)[1,1])

# ======================================================================
# 24. COMPARAÇÃO PRESENTE × CENÁRIO
# ======================================================================
resumo_comparacao <- dplyr::bind_rows(
  data.frame(
    species = sp_nome,
    scenario_id = "PRESENTE",
    period = "Presente",
    SSP = "Presente",
    threshold = threshold_modelo,
    threshold_value = threshold_value,
    suitable_area_Brazil_km2 = area_presente_nacional,
    #suitable_area_buffer_km2 = area_presente_buffer,
    suitable_area_current_vegetation_km2 = area_presente_vegetacao,
    n_suitable_pixels = sum(terra::values(raster_bin_nacional) == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  ),
  resumo_futuro
)

resumo_comparacao$change_Brazil_percent <- ((resumo_comparacao$suitable_area_Brazil_km2 / area_presente_nacional) - 1) * 100
resumo_comparacao$change_buffer_percent <- ((resumo_comparacao$suitable_area_buffer_km2 / area_presente_buffer) - 1) * 100
resumo_comparacao$change_vegetation_percent <- if (is.finite(area_presente_vegetacao) && area_presente_vegetacao > 0) ((resumo_comparacao$suitable_area_current_vegetation_km2 / area_presente_vegetacao) - 1) * 100 else NA_real_
resumo_comparacao$change_pixels_percent <- ((resumo_comparacao$n_suitable_pixels / sum(terra::values(raster_bin_nacional) == 1, na.rm = TRUE)) - 1) * 100
resumo_comparacao$change_Brazil_percent[resumo_comparacao$scenario_id == "PRESENTE"] <- 0
resumo_comparacao$change_buffer_percent[resumo_comparacao$scenario_id == "PRESENTE"] <- 0
resumo_comparacao$change_vegetation_percent[resumo_comparacao$scenario_id == "PRESENTE"] <- 0
resumo_comparacao$change_pixels_percent[resumo_comparacao$scenario_id == "PRESENTE"] <- 0

# ======================================================================
# 25. SALVAR TABELAS
# ======================================================================

arquivo_resumo <- file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_resumo.csv"))
arquivo_comparacao <- file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_comparacao_presente.csv"))
arquivo_excel <- file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_comparacao_presente.xlsx"))

readr::write_csv(resumo_futuro,arquivo_resumo)
readr::write_csv(resumo_comparacao,arquivo_comparacao)
openxlsx::write.xlsx(resumo_comparacao,file=arquivo_excel,overwrite=TRUE)

# ======================================================================
# 26. WORKFLOW
# ======================================================================

info_workflow <- data.frame(
  species=sp_nome,
  scenario_id=cenario$id,
  period=cenario$periodo,
  SSP=cenario$ssp,
  PCA_present_restricted_to_M=TRUE,
  PCA_projection_restricted_to_M=restric_pca_proj_futuro,
  number_PCs=length(pcs_modelo),
  PCs_used_in_model=paste(pcs_modelo,collapse="; "),
  predictors_used_in_model=paste(preditores_modelo,collapse="; "),
  threshold=threshold_modelo,
  threshold_value=threshold_value,
  model_algorithm="MaxEnt tuned with flexsdm::tune_max",
  future_model_fitted_again=FALSE,
  PCA_refitted_for_future=FALSE,
  coast_distance_outside_PCA="dist_coast_100km" %in% preditores_modelo,
  current_vegetation_mask_applied_to_future=aplicar_mascara_vegetacao_futuro,
  extrapolation_calculated=calcular_extrapolacao_futura,
  terra_memfrac=0.30,
  stringsAsFactors=FALSE
)

readr::write_csv(info_workflow,file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_workflow.csv")))

# ======================================================================
# 27. MAPA
# ======================================================================
arquivo_mapa <- file.path(dir_futuro,paste0(nome_sp_arquivo,"_",cenario$id,"_conferencia.png"))

png(filename=arquivo_mapa,width=1800,height=800,res=200)
par(mfrow=c(1,2),mar=c(3,3,3,5))
plot(raster_cont_fut,main=paste(sp_nome,"\n",cenario$periodo,cenario$ssp))
plot(raster_bin_fut,main=paste("Binário -",cenario$periodo,cenario$ssp))
dev.off()

# ======================================================================
# 28. RESULTADO
# ======================================================================
cat("\n============================================================\n")
cat("CENÁRIO FINALIZADO COM SUCESSO\n")
cat("============================================================\n")
cat("Espécie: ",sp_nome,"\n",sep="")
cat("Cenário: ",cenario$id,"\n",sep="")
cat("Período: ",cenario$periodo,"\n",sep="")
cat("SSP: ",cenario$ssp,"\n",sep="")
cat("Área adequada nacional: ",round(area_nacional,2)," km²\n",sep="")
cat("Área adequada no buffer: ",round(area_buffer,2)," km²\n",sep="")
if (is.finite(area_vegetacao)) cat("Área adequada na vegetação atual: ",round(area_vegetacao,2)," km²\n",sep="")
cat("Arquivos salvos em: ",dir_futuro,"\n",sep="")
cat("============================================================\n")

print(resumo_comparacao)

# ======================================================================
# 29. LIMPEZA FINAL
# ======================================================================

cat("\nLiberando memória...\n")

objetos_limpar <- c(
  "pca_fut",
  "pred_fut",
  "env_modelo_fut",
  "raster_cont_fut",
  "raster_bin_fut",
  "raster_cont_buffer_fut",
  "raster_bin_buffer_fut",
  "raster_cont_vegetacao_fut",
  "raster_bin_vegetacao_fut",
  "distancia_fut"
)

objetos_limpar <- objetos_limpar[vapply(objetos_limpar,exists,logical(1))]

if (length(objetos_limpar)>0) rm(list=objetos_limpar)

gc()

cat("\nMemória liberada.\n")
cat("\n============================================================\n")
cat("FIM DA EXECUÇÃO\n")
cat("============================================================\n")
cat("Para executar outro cenário, altere apenas:\n")
cat("cenario_escolhido <- \"2021_2040_SSP370\"\n")
cat("============================================================\n")

