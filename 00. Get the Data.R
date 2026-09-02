library("remotes")
install_github("LimaRAF/plantR", ref = "dev")
library("plantR")

usethis::edit_r_environ()

# ========================
# Função para baixar Eugenia
# ========================
baixar_eugenia <- function() {
  
  # pegar taxonKey para Eugenia
  key <- name_backbone(name = "Eugenia L.")$usageKey
  if (is.null(key)) stop("Taxon Eugenia não encontrado no backbone do GBIF")
  
  # iniciar o download
  job_key <- occ_download(
    pred("taxonKey", key),
    pred("basisOfRecord", "PRESERVED_SPECIMEN"),
    format = "DWCA"
  )
  
  cat("Download solicitado para Eugenia\n")
  cat("Chave do job:", job_key, "\n")
  
  # aguardar até o download ficar pronto
  status <- occ_download_meta(job_key)$status
  while (status != "SUCCEEDED") {
    cat("Status:", status, "- aguardando 30 segundos...\n")
    Sys.sleep(30)
    status <- occ_download_meta(job_key)$status
  }
  cat("✅ Download pronto!\n")
  
  # baixar o arquivo
  res_get <- occ_download_get(job_key, overwrite = TRUE)
  
  # importar os dados
  dados <- occ_download_import(res_get)
  
  cat("Total de registros importados:", nrow(dados), "\n")
  
  # salvar CSV
  write.csv(dados, "Eugenia_GBIF.csv", row.names = FALSE)
  cat("📂 Dados salvos em: Eugenia_GBIF.csv\n")
  
  return(dados)
}

# ========================
# Rodar para Eugenia
# ========================
dados_eugenia <- baixar_eugenia()

# conferir primeiras linhas
head(dados_eugenia)


dados <- read.csv("Eugenia_GBIF.csv", sep = ",", stringsAsFactors = FALSE)
View(dados)

# ver nomes atuais
names(dados)

# renomear apenas uma coluna
names(dados)[names(dados) == "countryCode"] <- "country"

# conferir
names(dados)

occs <- formatDwc(gbif_data = dados, drop = TRUE)
occs <- getCode(occs)

#Arrumar coletores
occs$recordedBy.new <- prepName(occs$recordedBy,
                                output = "first",
                                sep.out = "; ")
occs$recordedBy.aux <- prepName(occs$recordedBy,
                                output = "aux",
                                sep.out = "; ")

#Arrumar determinadores
occs$identifiedBy.new <- prepName(occs$identifiedBy,
                                  output = "first",
                                  sep.out = "; ")
occs$identifiedBy.aux <- prepName(occs$identifiedBy,
                                  output = "aux",
                                  sep.out = "; ")

#Arrumar quando não tem coletor ou determinador
occs$recordedBy.new <- missName(occs$recordedBy.new,
                                type = "collector", noName = "s.n.")
occs$identifiedBy.new <- missName(occs$identifiedBy.new,
                                type = "identificator", noName = "s.n.")

#Número do coletor
occs$recordNumber.new <- colNumber(occs$recordNumber, noNumb = "s.n.")

occs$year.new <- getYear(occs$year, noYear = "n.d.")
occs$yearIdentified.new <- getYear(occs$dateIdentified, noYear = "n.d.")
View(occs)

#Começando a arrumar a localização
occs <- fixLoc(occs,
               loc.levels = c("country", "stateProvince", "municipality", "locality"),
               scrap = TRUE)
locs <- strLoc(occs)
locs$loc.string <- prepLoc(locs$loc.string)
locs$loc.string1 <- prepLoc(locs$loc.string1)
locs$loc.string2 <- prepLoc(locs$loc.string2)

locs <- getLoc(locs, gazet = "plantR", gazet.names = c("loc", "loc.correct","latitude.gazetteer", "longitude.gazetteer",
                                                       "resolution.gazetteer"),
               orig.names = FALSE)

occs <- cbind.data.frame(occs, locs[, c("loc","loc.correct","latitude.gazetteer","longitude.gazetteer", "resolution.gazetteer")])

head(getAdmin(occs))

#Coordenadas
occs <- prepCoord(occs)
table(!is.na(occs$decimalLatitude.new))/dim(occs)[1]

occs <- getCoord(occs,
                lat.orig = "decimalLatitude",
                lon.orig = "decimalLongitude",
                lat.gazet = "latitude.gazetteer",
                lon.gazet = "longitude.gazetteer",
                res.gazet = "resolution.gazetteer",
                lat.new = "decimalLatitude.new",
                lon.new = "decimalLongitude.new",
                rm.gazet = FALSE)

table(!is.na(occs$decimalLatitude.new))/dim(occs)[1]

#Informações taxônomicas
occs <- fixSpecies(occs)
sort(table(occs$scientificName))
sort(table(occs$scientificName.new))
sort(table(occs$scientificNameStatus))

occs <- prepSpecies(occs, db = c("bfo"), sug.dist = 0.85)
table(occs$suggestedName)
table(occs$tax.notes)
View(occs)

#Validação de dados
occs <- validateLoc(occs, res.orig = "resol.orig", res.gazet = "resolution.gazetteer")
unique(occs$loc.check)

names(occs)[duplicated(names(occs))]
names(occs) <- make.unique(names(occs))

occs <- checkCoord(occs, keep.cols = c("geo.check", "NAME_0", "country.gazet"))
unique(occs$geo.check)           

occs <- checkBorders(occs, output = 'same.col')
unique(occs$geo.check)

occs <- checkShore(occs, output = 'same.col')
unique(occs$geo.check)

occs <- checkInverted(occs, output = 'same.col')
unique(occs$geo.check)

#Ver se tem materiais cultivados
occs <- getCult(occs)
table(occs$cult.check)

#Outilier espaciais
occs <- checkOut(occs, tax.name = "scientificName.new")

#Confiança de det
occs <- validateTax(occs)
table(occs$tax.check)/dim(occs)[1]

occs1 <- validateTax(occs, miss.taxonomist = c("Myrtaceae_Maruyama, A.", "Myrtaceae_Snow", "Myrtaceae_Hatschbach, G.", "Myrtaceae_Souza, V.C.", "Myrtaceae_Fernandes, T.", "Myrtaceae_Maruyama"))
table(occs1$tax.check)/dim(occs)[1]

#Duplicatas
occs1$numTombo <- getTombo(occs1[, "collectionCode.new"],
                           occs1[, "catalogNumber"])

dups <- prepDup(occs1,
                comb.fields = list(c("family", "col.last.name", "col.number", "col.loc"), c("family",
                                                                                            "col.year", "col.number", "col.loc"), c("species", "col.last.name", "col.number",
                                                                                                                                    "col.year"), c("col.year", "col.last.name", "col.number", "col.loc")))

dups <- getDup(dups)
head(dups[order(dups$dup.ID),],6)
table(dups$dup.prop)

occs1 <- cbind.data.frame(occs1,
                          dups[,c("dup.ID","dup.numb","dup.prop")],
                          stringsAsFactors = FALSE)

occs1 <- mergeDup(occs1)
table(occs1$tax.check, occs1$tax.check1)
View(occs1)

dim(occs1)
occs1 <- rmDup(occs1)

#Exportar
summ <- summaryData(occs1)
flags <- summaryFlags(occs1)
spp.check <- checkList(occs1)

saveData(occs, by = "family")

saveData(
  occs1,
  file.name = "output",
  dir.name = "",
  path = "",
  by = NULL,
  file.format = "csv",
  compress = FALSE,
  rm.dup = FALSE
)
