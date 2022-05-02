# ======================================================================
# ===    Realizar el ACP para todas las combinaciones de mosaicos    ===
# ======================================================================
#
library('raster')

# 1. Funciones para generar un mosaico por cada combinacion y por cada zona
# _________________________________________________________________________

# 1.1. Generar la direccion a las bandas de cada mosaico
# - data: Fila del csv con las combinaciones
# - zona: Nombre de la carpeta con las bandas
# - bandas: Named vector con el nombre de las bandas (sin prefijos ni sufijos)
getPaths <- function(data, zona, bandas){

  combi_id <- data['id']

  # Seleccionar los nombres de las bandas en la combinacion
  bandnames <- data[['band_name']]
  bandnames <- strsplit(bandnames, " ")

  # Comprobar pretratamiento
  # len = 2 - not thermal
  # len = 3 - thermal
  treat <- data[['band_pretreatmeant']]

  if(treat != "None"){
    prefix <- "clipped"
    suffix <- "_resize_clip.tif"
  } else {
    # Coger imgs originales, no hay termica en el mosaico
    prefix <- ""
    suffix <- ".tif"
  }

  # Vector donde almacenar las direcciones
  paths <- c()

  # Path generico para la zona
  input_path <- file.path('E:',"arqueologia_2022","celsa_dron",
                          "Yacimiento_Lepida_Celsa","Lepida Celsa",
                          "Mapas de reflectancia", zona, prefix)
  # Path especifico para cada banda
  for(band in bandnames){
    filename <- paste("Lepida_Celsa_",zona,"_",bandas[band],suffix,sep="")
    path <- file.path(input_path,filename)
    paths <- append(paths, path)
  }

  # Directorio de salida para el ACP
  output_path <- file.path('E:',"arqueologia_2022","celsa_dron",
                          "Yacimiento_Lepida_Celsa","Lepida Celsa",
                          "Mapas de reflectancia", zona, 'ACP')

  # Nombre del mosaico con los ACP y sus estadisticas
  filename <- paste0("Lepida_Celsa_",zona,"_ACP_",combi_id)

  return(list(
    "paths" = paths,
    "bandnames" = bandnames,
    "outputPath" = output_path,
    "outputFilename"=filename
  ))
}

# 1.2. Crea el mosaico con cada uno de los paths
# NOTA: En esta version se guarda el archivo con el mosaico en el directorio
# del script, pero se borra automaticamente despues de haber ejcutado el ACP
createMosaic <- function(paths){
  mosaic <- raster::stack(paths)
  # Output as RasterBrick
  return(brick(mosaic))
}

# 2. Funcion que realiza el ACP por cada mosaico
# _________________________________________________________________________
# Mismos parametros que la funcion pathResults
pcanalysis <- function(combis, zona, bandas){

  # Generar informacion a partir de las bandas integradas en la combinacion
  pathResults <- getPaths(combis, zona, bandas)

  bandnames <- pathResults[['bandnames']][[1]]
  output_path <- pathResults[['outputPath']]
  filename <- pathResults[['outputFilename']]

  cat("\n\nCreate ACP ",filename,
      "\n_________________________________________________________\n\n")

  cat("..Crear mosaico\n")
  # Generar mosaico
  mosaic <- createMosaic(pathResults[['paths']])

  # Transformar la imagen en un DF
  cat("..Transformarlo en un DF\n")
  mosaic.df <- as.data.frame(mosaic)
  names(mosaic.df) <- bandnames
  cat("..Convertir NAs en valores 0\n")
  mosaic.df[is.na(mosaic.df)] <- 0

  # Run PCA
  cat("..Realizar ACP\n")
  pca <- prcomp(mosaic.df, scale.=TRUE)

  # Obtener el numero de bandas de la imagen
  nBands <- length(bandnames)

  # Save PCs
  cat("..Guardar componentes\n")
  #mosaic.pc <- pca$x

  # Obtener un vector con el nombre de los PCs
  colNamesPC <- retrievePCnames(pca)

  # UPDATE: Los nombres de las columnas se guardan con un nombre estandar
  # al utilizar el metodo data.frame[[]]
  # Para conseguir que se renombren primero se guardan los nombres de las
  # bandas integradas en el mosaico, y luego se unen los de los PC
  ancientColumnNames <- names(mosaic)
  newColumnNames <- append(ancientColumnNames, colNamesPC)
  
  # Insert each PC inside the original raster
  for(pc in seq(1,nBands)){
    cat('....Componente ',pc,'\n')
    newColumn <- colNamesPC[pc]
    mosaic[[newColumn]] <- pca$x[,pc]
  }
  
  names(mosaic) <- newColumnNames
  
  # Exportar estadisticas
  cat("..Exportar estadisticas\n")
  writeStats(pca, output_path, filename, colNamesPC)
  
  cat("..Eliminar objetos innecesarios\n")
  rm(pca, mosaic.df)
  
  # Nombrar al nuevo mosaico con los ACP
  out_filename = paste0(filename,'.tif')
  
  # Exportar imagen solo con los ACP
  cat("..Exportar imagen con los CPs\n")
  raster::writeRaster(
    subset(mosaic,(nBands+1):length(names(mosaic))),
    file.path(output_path,out_filename),
    bylayer = F, options=c("COMPRESS=LZW", "PREDICTOR=3")
  )
}

# 3. Funciones para exportar resultados
# _________________________________________________________________________
# 3.1. Obtener vector con los nombres de los PCs del mosaico
retrievePCnames <- function(pca){
  # Obtener numero de CP
  nPC <- length(pca$sdev)
  # Generar los nombres de las columnas a partir de ellos
  names <- c()
  for(i in seq(1,nPC)){
    name <- paste0("PC",i)
    names <- append(names, name)
  }

  return(names)
}

# 3.2. Exportar stats del ACP
writeStats <- function(pca,output_path,filename,statsColnames){

  # Create loading matrix
  loadings <- pca$rotation
  # Save loading matrix
  write.csv(loadings,
            file=file.path(output_path,
                           paste0(filename,"_loadings.csv")))

  # Obtener las stats:
  # Eigenvalues - (varianza de cada PC, obtenida elevando al cuadrado su desvest)
  eigs <- pca$sdev^2

  # Desviacion estandar - raiz cuadrada de la varianza
  SD = pca$sdev
  #Porcentaje de varianza
  proportion = eigs/sum(eigs)
  cumulative = cumsum(eigs)/sum(eigs)

  pca_stats <- rbind(
    SD = SD,
    Proportion = proportion,
    Cumulative = cumulative
  )

  # Renombrar columnas con los PCs
  colnames(pca_stats) <- statsColnames

  # Exportar
  write.csv(pca_stats,
            file=file.path(output_path,
                           paste0(filename,"_stats.csv")))

}

# 4. Crear el ACP por cada mosaico
# _________________________________________________________________________

# Importar el archivo csv
ids_path <- file.path('E:',"arqueologia_2022","celsa_dron","scripts","GDAL_Python",
                      "mosaics_ids.csv")
mosaic_ids <- read.csv(ids_path, header=T)

# Escribir las bandas
bandas = c(
  "MS_index_ndvi",
  "MS_transparent_reflectance_green",
  "MS_transparent_reflectance_red",
  "MS_transparent_reflectance_red_edge",
  "MS_transparent_reflectance_nir",
  "transparent_reflectance_thermal_ir_withoutOutlayers"
)

names(bandas) <- c("ndvi","green","red","red_edge","nir","thermal")

# Escribir las zonas
zonas = c("Yacimiento")
#zonas = c("Campos","Yacimiento")

for(zona in zonas){
  for(i in seq(42,length(mosaic_ids[,1]))){ 
    inicio <- Sys.time()
    # Save row
    row <- mosaic_ids[i,]
    # Compute PCA
    pcanalysis(row,zona,bandas)
    # Tiempo de ejecucion
    print(inicio-Sys.time())
    # Clean storage
    gc()
  }
  remove(row)
}
