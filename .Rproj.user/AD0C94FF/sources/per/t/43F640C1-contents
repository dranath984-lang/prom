source("instalar_paquetes.R")  # Ajusta la ruta si está en una subcarpeta

instalar_y_cargar(
  paquetes_cran = c("fs", "here", "tidyverse", "data.table", "pheatmap", "writexl", "readxl", "skimr", "ggvenn", "factoextra", "stats","corrplot", "DataExplorer", "janitor", "dplyr

", "stringr", "tidyr", "readr", "ggplot2", "conflicted"),
  paquetes_bioc = c("biomartr", "clusterProfiler")
)


#######################################################
#opcion para forzar la instalación 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("Biostrings")
########################resolucion de conflictos 
 


conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("summarise", "dplyr")



########################################################




getwd()
dir_create("datos")        # Para archivos TSV
dir_create("scripts")      # Para scripts .R
dir_create("salidas")      # Para resultados
dir_create("salidas_data") # Para datos procesados    


# Nombres de archivo (sin extensión)
archivos <- c(
  "DE_EV",
  "DE_PHR-L7-RNAi",
  "DE_PHR-L7-Ox",
  "DE_PHR-L7-RNAi-LPivsEVLPi",
  "DE_PHR-L7Ox-LPivsEVLPi",
  "DE_PHR-L7RNAi_OPivsEV_OPi",
  "DE_PHR-L7Ox_OPivsEV_OPi"
)

# Leer todos los archivos y guardarlos en una lista
datos_list <- lapply(archivos, function(nombre) {
  read_excel(paste0("datos/", nombre, ".xlsx")) 
})


# Crear nombres de objetos más cortos y seguros
nombres_limpios <- c(
  "ev",
  "rna",
  "ox",
  "rna_lp",
  "ox_lp",
  "rna_op",
  "ox_op"
)

# Asignar nombres a la lista
names(datos_list) <- nombres_limpios

###################

datos_list$ev       # DE_EV.xlsx
datos_list$rna      # DE_PHR-L7-RNAi.xlsx
datos_list$ox
datos_list$rna_lp 
datos_list$ox_lp    # DE_PHR-L7Ox-LPivsEVLPi.xlsx
datos_list$rna_op
datos_list$ox_op  
  
  

for (i in seq_along(nombres_limpios)) {
  assign(paste0("dat_", nombres_limpios[i]), datos_list[[i]])
}


data.frame(objeto = paste0("dat_", nombres_limpios), archivo = archivos)


#####################################

# Vector con nombres de tus objetos
objetos <- c("dat_ev", "dat_rna", "dat_ox", "dat_rna_lp", "dat_ox_lp", "dat_rna_op", "dat_ox_op")

# Recorrer cada objeto y mostrar información clave
for (obj in objetos) {
  cat("\n🔍 Explorando:", obj, "\n")
  
  # Obtener el objeto
  tabla <- get(obj)
  
  # Mostrar dimensiones
  cat("📐 Dimensiones:", dim(tabla)[1], "filas x", dim(tabla)[2], "columnas\n")
  
  # Mostrar nombres de columnas
  cat("🧬 Columnas:", paste(names(tabla), collapse = ", "), "\n")
  
  # Mostrar tipos de datos
  cat("🔤 Estructura:\n")
  str(tabla)
  
  # Mostrar primeras filas
  cat("👀 Primeras filas:\n")
  print(head(tabla))
  
  # Mostrar resumen estadístico
  cat("📊 Resumen:\n")
  print(summary(tabla))
  
  # Mostrar resumen extendido si skimr está disponible
  if ("skimr" %in% rownames(installed.packages())) {
    cat("📋 Resumen extendido con skimr:\n")
    print(skimr::skim(tabla))
  }
  
  # Abrir en pestaña si estás en RStudio
  if (interactive()) {
    View(tabla)
  }
}


#############################################
# Definir función para limpiar una tabla
clean_table <- function(tabla, origen) {
  # Paso 1: Renombrar columnas clave basado en posiciones observadas (col1 = GeneID, col8 = log2FC, col9 = padj)
  # Asumiendo estructura: col1 título/GeneID, col2-7 réplicas, col8 log2FC, col9 padj (ajusta si varía)
  nombres <- names(tabla)
  tabla <- tabla %>%
    rename(
      GeneID = 1,          # Primera columna: título largo, pero fila 2+ son IDs
      log2FC = 8,          # Columna de log2FoldChange
      padj = 9             # Columna de padj
    ) %>%
    # Paso 2: Convertir a numérico (están como chr)
    mutate(
      log2FC = as.numeric(log2FC),
      padj = as.numeric(padj)
    ) %>%
    # Paso 3: Drop réplicas y extras (mantener solo clave)
    select(GeneID, log2FC, padj) %>%
    # Paso 4: Filtrar filas inválidas (e.g., si primera fila es header texto, o NA)
    filter(!is.na(GeneID) & str_detect(GeneID, "^Phvul\\."), na.rm = TRUE) %>%  # Solo IDs válidos (empiezan con Phvul.)
    # Paso 5: Agregar origen
    mutate(tabla_origen = origen)
  
  return(tabla)
}

# Lista de objetos originales y sus orígenes (códigos cortos para simplicidad)
tablas <- list(
  dat_ev = "ev",          # Empty Vector OPi vs LPi
  dat_rna = "rna",        # RNAi OPi vs LPi
  dat_ox = "ox",          # OX OPi vs LPi
  dat_rna_lp = "rna_lp",  # RNAi LPi vs EV LPi
  dat_ox_lp = "ox_lp",    # OX LPi vs EV LPi
  dat_rna_op = "rna_op",  # RNAi OPi vs EV OPi
  dat_ox_op = "ox_op"     # OX OPi vs EV OPi
)

# Aplicar la función a cada tabla y guardar en lista
tablas_limpias <- map2(tablas, names(tablas), ~ clean_table(get(.y), .x))

# Nombrar la lista para fácil acceso (e.g., tablas_limpias$dat_ev)
names(tablas_limpias) <- names(tablas)

# Verificación rápida: Dimensiones y head de una (e.g., dat_ev limpia)
dim(tablas_limpias$dat_ev)
head(tablas_limpias$dat_ev)

# Guardar la lista completa como RDS en salidas_data (para reuse en próximos scripts)
saveRDS(tablas_limpias, here("salidas_data", "tablas_limpias.rds"))

# Opcional: Exportar cada limpia como CSV individual para chequeo manual
walk2(tablas_limpias, names(tablas_limpias), ~ write_csv(.x, here("salidas_data", paste0(.y, "_clean.csv"))))

###########################################################
# scripts/03_convert_csv_to_xlsx.R
# Objetivo: Convertir CSVs limpios a XLSX para legibilidad externa.

# Path a CSVs
###csv_dir <- here("salidas_data")

# Listar CSVs limpios
###csv_files <- dir_ls(csv_dir, regexp = "_clean\\.csv$")

# Convertir cada CSV a XLSX
###for (csv_path in csv_files) {
  # Cargar CSV
###  tabla <- read_csv(csv_path, show_col_types = FALSE)
  
  # Generar nombre XLSX (e.g., dat_ev_clean.xlsx)
###  xlsx_name <- str_replace(basename(csv_path), "\\.csv$", ".xlsx")
###  xlsx_path <- here("salidas_data", xlsx_name)
  
  # Escribir a XLSX
###  write_xlsx(tabla, xlsx_path)
  
###  cat("Convertido:", csv_path, "a", xlsx_path, "\n")
###}

##############################################################

# scripts/02_explore_clean_csvs.R
# Objetivo: Explorar características de los CSVs limpios en salidas_data (similar a loop original).

# Definir path a carpeta de CSVs
csv_dir <- here("salidas_data")

# Listar automáticamente los CSVs que terminan en "_clean.csv" (ignora otros como .rds)
csv_files <- fs::dir_ls(csv_dir, regexp = "_clean\\.csv$")  # Usa fs para listar

# Extraer nombres base (e.g., "dat_ev" de "dat_ev_clean.csv")
objetos <- basename(csv_files) %>% str_remove("_clean\\.csv$")

# Recorrer cada CSV y mostrar información clave
for (i in seq_along(csv_files)) {
  obj <- objetos[i]
  file_path <- csv_files[i]
  
  cat("\n🔍 Explorando CSV:", obj, "(de", file_path, ")\n")
  
  # Cargar el CSV
  tabla <- read_csv(file_path, show_col_types = FALSE)  # Silencia mensajes
  
  # Mostrar dimensiones
  cat("📐 Dimensiones:", nrow(tabla), "filas x", ncol(tabla), "columnas\n")
  
  # Mostrar nombres de columnas
  cat("🧬 Columnas:", paste(names(tabla), collapse = ", "), "\n")
  
  # Mostrar tipos de datos
  cat("🔤 Estructura:\n")
  str(tabla)
  
  # Mostrar primeras filas
  cat("👀 Primeras filas:\n")
  print(head(tabla))
  
  # Mostrar resumen estadístico
  cat("📊 Resumen:\n")
  print(summary(tabla))
  
  # Mostrar resumen extendido si skimr está disponible
  if ("skimr" %in% rownames(installed.packages())) {
    cat("📋 Resumen extendido con skimr:\n")
    print(skimr::skim(tabla))
  }
  
  # Abrir en pestaña si estás en RStudio
  if (interactive()) {
    View(tabla)
  }
}

###########################################################################################################################
# Paso a: Renombrar columnas (versión con bucle for a prueba de errores)

# 1. Crear una lista vacía para guardar los resultados
tablas_renombradas <- list()

# 2. Iterar sobre los NOMBRES de la lista original
for (nombre_tabla in names(tablas_limpias)) {
  
  cat("Procesando tabla:", nombre_tabla, "\n") # Mensaje para ver el progreso
  
  # 3. Obtener la tabla actual
  tabla_actual <- tablas_limpias[[nombre_tabla]]
  
  # 4. Renombrar y seleccionar columnas en esta tabla
  tabla_procesada <- tabla_actual %>%
    rename_with(~ paste0(., "_", nombre_tabla), c("log2FC", "padj")) %>%
    select(GeneID, starts_with("log2FC_"), starts_with("padj_"))
  
  # 5. Guardar la tabla procesada en la nueva lista
  tablas_renombradas[[nombre_tabla]] <- tabla_procesada
}

# El resto de tu script continúa exactamente igual...
df_unificado <- reduce(tablas_renombradas, full_join, by = "GeneID")










# Paso b: Merge secuencial con reduce y full_join por GeneID
df_unificado <- reduce(tablas_renombradas, full_join, by = "GeneID")

# Paso c: Manejar duplicados GeneID (si hay, aunque raro)
# Chequea y resume (e.g., media de log2FC si duplicado)
duplicados <- df_unificado %>%
  group_by(GeneID) %>%
  filter(n() > 1)

if (nrow(duplicados) > 0) {
  cat("⚠️ Duplicados encontrados:", nrow(duplicados), "- Resumiendo con media.\n")
  df_unificado <- df_unificado %>%
    group_by(GeneID) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))  # Media para num, o ajusta
} else {
  cat("✅ No duplicados en GeneID.\n")
}

# Verificación post-merge
cat("📐 Dimensiones unificado:", nrow(df_unificado), "filas x", ncol(df_unificado), "columnas\n")
cat("🧬 Columnas:", paste(names(df_unificado), collapse = ", "), "\n")
cat("👀 Primeras filas:\n")
print(head(df_unificado))
cat("📊 Resumen:\n")
print(summary(df_unificado))
# Opcional: skim(df_unificado) si skimr cargado
if ("skimr" %in% rownames(installed.packages())) {
  print(skimr::skim(df_unificado))
}

# Abrir en View para inspección
if (interactive()) {
  View(df_unificado)
}

# Paso d: Guardar unificado
saveRDS(df_unificado, here("salidas_data", "df_unificado.rds"))
write_csv(df_unificado, here("salidas_data", "df_unificado.csv"))  # Opcional CSV para backup
cat("💾 Guardado en salidas_data/df_unificado.rds y .csv\n")

##########################################################################################

sessionInfo()

##########################################################################################
#Diseño del Script para Crear la Matriz log2FC por GeneID/Tabla

# scripts/04_create_matrices.R
# Objetivo: Crear matriz log2FC (y opcional padj) por GeneID/tabla a partir de df_unificado.

# Cargar paquetes
library(tidyverse)  # Para select, mutate, etc.
library(here)       # Para paths

# Cargar df_unificado desde RDS
df_unificado <- readRDS(here("salidas_data", "df_unificado.rds"))

# Paso 1: Crear matriz log2FC
# Selecciona GeneID + todas columnas que empiezan con "log2FC_"
matriz_log2FC <- df_unificado %>%
  select(GeneID, starts_with("log2FC_")) %>%  # Matriz: rows=GeneID, cols=tablas
  # Opcional: Manejar NA (e.g., imputar 0 para no-DEGs, pero cuidado—puede bias; deja NA para ahora)
  # mutate(across(starts_with("log2FC_"), ~ replace_na(., 0))) %>%
  arrange(GeneID)  # Ordena por GeneID para consistencia

# Verificación
cat("📐 Dimensiones matriz_log2FC:", nrow(matriz_log2FC), "filas x", ncol(matriz_log2FC), "columnas\n")
cat("🧬 Columnas:", paste(names(matriz_log2FC), collapse = ", "), "\n")
cat("👀 Primeras filas:\n")
print(head(matriz_log2FC))
cat("📊 Resumen:\n")
print(summary(matriz_log2FC))
# Opcional skim
if ("skimr" %in% rownames(installed.packages())) {
  print(skimr::skim(matriz_log2FC))
}

# Abrir en View
if (interactive()) {
  View(matriz_log2FC)
}

# Opcional: Crear matriz padj similar (para filtrar overlaps por significancia)
matriz_padj <- df_unificado %>%
  select(GeneID, starts_with("padj_")) %>%
  # mutate(across(starts_with("padj_"), ~ replace_na(., 1))) %>%  # NA=1 (no sig) si quieres
  arrange(GeneID)

# Verificación similar para padj (repite si usas)
cat("📐 Dimensiones matriz_padj:", nrow(matriz_padj), "filas x", ncol(matriz_padj), "columnas\n")
# ... (agrega head/summary/skim similar)

# Guardar matrices
saveRDS(matriz_log2FC, here("salidas_data", "matriz_log2FC.rds"))
write_csv(matriz_log2FC, here("salidas_data", "matriz_log2FC.csv"))
saveRDS(matriz_padj, here("salidas_data", "matriz_padj.rds"))  # Opcional
write_csv(matriz_padj, here("salidas_data", "matriz_padj.csv"))  # Opcional
cat("💾 Matrices guardadas en salidas_data.\n")

##################################################################################
#####################################################################

# scripts/05_venn_overlaps.R
# Objetivo: Filtrar genes up, crear listas por grupo (LPi vs OPi), calcular overlaps/Venn.

# Cargar matrices desde RDS
matriz_log2FC <- readRDS(here("salidas_data", "matriz_log2FC.rds"))
matriz_padj <- readRDS(here("salidas_data", "matriz_padj.rds"))

# Definir thresholds para "up" genes (ajusta según necesidad, e.g., log2FC >1 = 2x up, padj <0.05)
fc_threshold <- 0.5    # log2FC >1
padj_threshold <- 0.05

# Definir grupos de tablas (basado en LPi vs OPi focus)
grupo_LPi <- c("log2FC_dat_ev", "log2FC_dat_rna_lp", "log2FC_dat_ox_lp")  # LPi-related
grupo_OPi <- c("log2FC_dat_rna_op", "log2FC_dat_ox_op")                   # OPi-related
# Opcional: Agrega más grupos, e.g., RNAi: "log2FC_dat_rna", "log2FC_dat_rna_lp", "log2FC_dat_rna_op"

# Función para filtrar genes up en un grupo de columnas
get_up_genes <- function(mat_fc, mat_padj, grupo_cols) {
  mat_fc %>%
    filter(if_any(all_of(grupo_cols), ~ . > fc_threshold) &  # Up en TODAS columnas del grupo (stricter; usa if_any para ANY)
             if_any(matches(str_replace(grupo_cols, "log2FC", "padj")), ~ . < padj_threshold)) %>%  # Sig en padj correspondientes
    pull(GeneID) %>%  # Extrae lista de GeneID
    unique()          # Evita dups
}

# Crear listas de genes up por grupo
genes_up_LPi <- get_up_genes(matriz_log2FC, matriz_padj, grupo_LPi)
genes_up_OPi <- get_up_genes(matriz_log2FC, matriz_padj, grupo_OPi)

# Verificación: Conteos iniciales
cat("Genes up en LPi grupo:", length(genes_up_LPi), "\n")
cat("Genes up en OPi grupo:", length(genes_up_OPi), "\n")

# Paso para overlaps/counts (simple intersección)
overlap_LPi_OPi <- intersect(genes_up_LPi, genes_up_OPi)  # Compartidos
unique_LPi <- setdiff(genes_up_LPi, genes_up_OPi)         # Únicos LPi
unique_OPi <- setdiff(genes_up_OPi, genes_up_LPi)         # Únicos OPi

cat("Overlap (compartidos):", length(overlap_LPi_OPi), "\n")
cat("Únicos LPi:", length(unique_LPi), "\n")
cat("Únicos OPi:", length(unique_OPi), "\n")

# Crear matriz binaria de presencia (1=up en tabla, 0=no) para análisis avanzado
# Usa matriz_log2FC como base, binariza basado en thresholds
matriz_binaria <- matriz_log2FC %>%
  mutate(across(starts_with("log2FC_"), 
                ~ ifelse(. > fc_threshold & 
                           matriz_padj[[str_replace(cur_column(), "log2FC", "padj")]] < padj_threshold, 1, 0)))

# Verificación matriz binaria
cat("📐 Dimensiones matriz_binaria:", nrow(matriz_binaria), "filas x", ncol(matriz_binaria), "columnas\n")
print(head(matriz_binaria))
print(summary(matriz_binaria))  # Debería mostrar 0/1/NA

# Venn plot simple (para 2 grupos; extiende a más con list)
venn_list <- list(LPi = genes_up_LPi, OPi = genes_up_OPi)
ggvenn(venn_list, fill_color = c("blue", "green"), show_percentage = TRUE) +
  ggtitle("Venn: Up genes LPi vs OPi")
ggsave(here("salidas", "venn_LPi_vs_OPi.png"))  # Guarda plot

# Guardar outputs
saveRDS(venn_list, here("salidas_data", "venn_lists.rds"))  # Listas para reuse
write_csv(matriz_binaria, here("salidas_data", "matriz_binaria.csv"))  # Matriz para análisis
cat("💾 Outputs guardados en salidas_data y salidas.\n")

###############################################################################

#Extraer listas de GeneID para promotores
write_csv(tibble(GeneID = unique_LPi), here("salidas_data", "candidatos_LPi_unicos_Venn.csv"))
write_csv(tibble(GeneID = overlap_LPi_OPi), here("salidas_data", "candidatos_compartidos_Venn.csv"))
write_csv(tibble(GeneID = unique_OPi), here("salidas_data", "candidatos_OPi_unicos_Venn.csv"))

################################################################################
# scripts/06_clustering.R
# Objetivo: Clustering jerárquico/k-means en matriz_binaria/log2FC para patrones up/down.

# Cargar paquetes
#library(factoextra) # Para fviz_nbclust (optimal k), fviz_dend, etc. (instala si no)
#library(stats)      # Para hclust, kmeans (base R)

# Cargar matrices
matriz_binaria <- read_csv(here("salidas_data", "matriz_binaria.csv"), show_col_types = FALSE)
# O usa matriz_log2FC si prefieres valores continuos (comenta binaria)
# matriz_log2FC <- readRDS(here("salidas_data", "matriz_log2FC.rds"))

# Preparación: Convertir a matrix numérica (rows=GeneID, cols=tablas)
# Maneja NA: Imputa 0 (asumir no up); filtra rows con >50% NA
mat_clust <- matriz_binaria %>%
  column_to_rownames(var = "GeneID") %>%  # Rows como GeneID
  mutate(across(everything(), ~ replace_na(., 0))) %>%  # Imputa NA=0 (no up)
  filter(rowSums(is.na(.)) / ncol(.) < 0.5) %>%  # Filtra rows >50% NA
  as.matrix()  # Para clustering

# Verificación
cat("📐 Dimensiones para clustering:", nrow(mat_clust), "filas x", ncol(mat_clust), "columnas\n")
print(head(mat_clust))

# Paso 1: Clustering Jerárquico (hclust)
dist_mat <- dist(mat_clust, method = "euclidean")  # Distancia (euclidean para binario)
hclust_res <- hclust(dist_mat, method = "ward.D2")  # Método ward para clusters compactos

# Plot dendrograma
fviz_dend(hclust_res, k = 5,  # Elige k inicial (ajusta)
          cex = 0.5, horiz = TRUE, main = "Dendrograma Jerárquico (k=5)") +
  theme_minimal()
ggsave(here("salidas", "dendrograma_hclust.png"))

# Extraer clusters (corta en k=5-10)
k_clusters <- 5  # Ajusta basado en dendrograma o elbow
clusters_h <- cutree(hclust_res, k = k_clusters)
# Agrega a df
df_clusters_h <- tibble(GeneID = rownames(mat_clust), cluster_h = clusters_h)

# Paso 2: K-means (alternativa no-jerárquica)
# Encuentra optimal k (elbow method)
fviz_nbclust(mat_clust, kmeans, method = "wss") +  # Within Sum Squares
  ggtitle("Elbow para Optimal K")
ggsave(here("salidas", "elbow_kmeans.png"))

# Aplica k-means con k elegido (e.g., del elbow ~4-6)
kmeans_res <- kmeans(mat_clust, centers = k_clusters, nstart = 25)  # nstart para estabilidad
clusters_k <- kmeans_res$cluster
df_clusters_k <- tibble(GeneID = rownames(mat_clust), cluster_k = clusters_k)

# Verificación clusters
cat("Conteos por cluster (hclust):\n")
print(table(df_clusters_h$cluster_h))
cat("Conteos por cluster (kmeans):\n")
print(table(df_clusters_k$cluster_k))

# Plot clusters (e.g., PCA colored by cluster)
fviz_cluster(kmeans_res, data = mat_clust, geom = "point", ellipse.type = "convex") +
  ggtitle("K-means Clusters")
ggsave(here("salidas", "plot_kmeans_clusters.png"))

# Guardar outputs
saveRDS(df_clusters_h, here("salidas_data", "clusters_hclust.rds"))
saveRDS(df_clusters_k, here("salidas_data", "clusters_kmeans.rds"))
write_csv(df_clusters_h, here("salidas_data", "clusters_hclust.csv"))  # GeneID por cluster
write_csv(df_clusters_k, here("salidas_data", "clusters_kmeans.csv"))
cat("💾 Clusters guardados en salidas_data.\n")


######################################################################################################
# scripts/07_correlations.R
# Objetivo: Correlaciones entre condiciones (e.g., tablas), identificar genes con FC invertido (RNAi vs OX).

# Cargar paquetes

#library(corrplot)   # Para heatmap correl (instala si no: install.packages("corrplot"))

# Cargar matrices y clusters
matriz_log2FC <- readRDS(here("salidas_data", "matriz_log2FC.rds"))
df_clusters_h <- read_csv(here("salidas_data", "clusters_hclust.csv"), show_col_types = FALSE)  # O kmeans

# Unir clusters a matriz para correl por grupo (opcional)
matriz_with_clusters <- matriz_log2FC %>%
  left_join(df_clusters_h, by = "GeneID")

# Paso 1: Matriz de correlaciones (Pearson entre columnas log2FC)
# Prepara numérica (imputa NA=0 para cor)
mat_cor <- matriz_log2FC %>%
  column_to_rownames("GeneID") %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  as.matrix()

cor_mat <- cor(mat_cor, method = "pearson")  # Correl entre tablas

# Plot heatmap correl
corrplot(cor_mat, method = "color", type = "upper", tl.cex = 0.8,
         title = "Heatmap Correlaciones log2FC", mar = c(0,0,1,0))
# Guarda manual o con png() si corrplot no ggsave

# Verificación
print(cor_mat)  # Matriz numérica (e.g., alta cor LPi-tablas = similares)

# Paso 2: Filtrar genes con FC invertido (e.g., up RNAi, down OX)
# Define grupos RNAi/OX (todas columnas relevantes)
grupo_RNAi <- c("log2FC_dat_rna", "log2FC_dat_rna_lp", "log2FC_dat_rna_op")
grupo_OX <- c("log2FC_dat_ox", "log2FC_dat_ox_lp", "log2FC_dat_ox_op")

genes_invertidos <- matriz_log2FC %>%
  filter(if_any(all_of(grupo_RNAi), ~ . > 1) &  # Up en al menos una RNAi
           if_any(all_of(grupo_OX), ~ . < -1) &  # Down en al menos una OX
           if_any(matches(str_replace(grupo_RNAi, "log2FC", "padj")), ~ . < 0.05) &  # Sig
           if_any(matches(str_replace(grupo_OX, "log2FC", "padj")), ~ . < 0.05)) %>%
  pull(GeneID)

# Verificación
cat("Genes con FC invertido (RNAi up, OX down):", length(genes_invertidos), "\n")
print(head(genes_invertidos))  # Lista IDs

# Opcional: Correl por cluster (e.g., en cluster 2)
for (clust in unique(df_clusters_h$cluster_h)) {
  mat_clust <- mat_cor[df_clusters_h$GeneID[df_clusters_h$cluster_h == clust], ]
  cor_clust <- cor(mat_clust)
  cat("Correl en cluster", clust, ":\n")
  print(cor_clust)
}

# Guardar outputs
saveRDS(list(genes_invertidos = genes_invertidos), here("salidas_data", "genes_invertidos.rds"))
write_csv(tibble(GeneID = genes_invertidos), here("salidas_data", "genes_invertidos.csv"))
saveRDS(cor_mat, here("salidas_data", "cor_matrix.rds"))
cat("💾 Outputs guardados.\n")

############################################### todo sirve

# scripts/08_annotate_geneids.R
# Objetivo: Anotar GeneIDs con funciones/GO (vía Phytozome/BioMart).

# Cargar paquetes
#library(biomartr)   # Para BioMart anot (instala BiocManager::install("biomartr") si no; alterna manual)

# Cargar datos previos (e.g., genes_invertidos, df_clusters_h)
genes_invertidos <- read_csv(here("salidas_data", "genes_invertidos.csv"), show_col_types = FALSE)$GeneID
df_clusters_h <- read_csv(here("salidas_data", "clusters_hclust.csv"), show_col_types = FALSE)

# Paso 1: Exportar listas GeneID para anotación manual/Phytozome
# e.g., sube a Phytozome > Search > Bulk download annotations (GO/functions por ID list)
write_lines(genes_invertidos, here("salidas_data", "geneids_invertidos.txt"))  # TXT para upload
write_lines(unique(df_clusters_h$GeneID[df_clusters_h$cluster_h == 2]), here("salidas_data", "geneids_cluster2.txt"))  # Por cluster

# Alternativa automatizada con biomartr (para Phaseolus vulgaris GO/functions; ajusta dataset)
# Nota: Phytozome no siempre en BioMart—si falla, usa manual
organism <- "Phaseolus vulgaris"  # Ajusta si en BioMart



#############################bien

# --- INICIA EL CÓDIGO FINAL DE ANOTACIÓN ---

## 1. Cargar el archivo de anotaciones descargado de Phytozome
# ----------------------------------------------------------------
# Asegúrate de que el archivo se llame 'phytozome_annotations.csv' y esté en la carpeta 'datos'.
annot_manual <- read_csv(here("datos", "phytozome_annotations.csv"))

# OJO: Es vital revisar que los nombres de las columnas se hayan leído bien.
# Ejecuta `colnames(annot_manual)` y verifica que coincidan con los que usamos abajo.
# Si se llaman "Gene Name" (con espacio), el código usará comillas invertidas (`Gene Name`).
print("Columnas cargadas:")
print(colnames(annot_manual))


## 2. Procesar y consolidar las anotaciones
# ----------------------------------------------------------------
# Agrupamos por transcrito y juntamos todas sus anotaciones en una sola línea.
annot_procesado <- annot_manual %>%
  group_by(`Transcript Name`) %>% # Agrupamos por el ID del transcrito
  summarise(
    # Juntamos todos los GO IDs únicos, separados por "; "
    GO_IDs = paste(unique(na.omit(`GO ID`)), collapse = "; "),
    # Hacemos lo mismo con las descripciones funcionales
    Descriptions = paste(unique(na.omit(Description)), collapse = " | ")
  ) %>%
  ungroup() # Es buena práctica desagrupar al final

cat("\n✅ Tabla de anotaciones procesada. Vista previa:\n")
print(head(annot_procesado))


## 3. Unir las anotaciones con tus datos originales
# ----------------------------------------------------------------
# Primero, con los datos de clusters
df_annotated_clusters <- df_clusters_h %>%
  # Unimos usando el ID del gen/transcrito.
  left_join(annot_procesado, by = c("GeneID" = "Transcript Name"))

# Segundo, con tu lista de genes invertidos
df_annotated_invertidos <- tibble(GeneID = genes_invertidos) %>%
  left_join(annot_procesado, by = c("GeneID" = "Transcript Name"))


## 4. Verificación y guardado final
# ----------------------------------------------------------------
cat("\n📐 Vista previa de los clusters anotados:\n")
# Usamos `glimpse` para ver todas las columnas de forma compacta
glimpse(df_annotated_clusters)

cat("\n💾 Guardando resultados anotados...\n")
write_csv(df_annotated_clusters, here("salidas_data", "annotated_clusters.csv"))
write_csv(df_annotated_invertidos, here("salidas_data", "annotated_invertidos.csv"))

cat("¡Proceso de anotación completado con éxito! ✨\n")

# --- FIN DEL CÓDIGO ---

# Filtra el data frame para ver solo las filas que SÍ tienen anotaciones
df_annotated_clusters %>% 
  filter(!is.na(GO_IDs))

# Esta tabla debería tener 40 filas y NO tener NAs (o muy pocos, si algunos genes no tenían anotación)
print(df_annotated_invertidos)

####################### 31-10-25

# =================================================================
# BLOQUE DE VERIFICACIÓN DE DATOS (CHEQUEO DE SANIDAD)
# =================================================================

cat("🤔 Iniciando chequeo de integridad de los datos de entrada...\n\n")

# --- 1. Verificación de 'df_annotated_clusters' ---
cat("📊 Verificando 'df_annotated_clusters'...\n")
glimpse(df_annotated_clusters)
genes_con_go_en_clusters <- sum(!is.na(df_annotated_clusters$GO_IDs) & df_annotated_clusters$GO_IDs != "")
cat(" -> Total de genes con anotación GO en clusters:", genes_con_go_en_clusters, "\n\n")

# --- 2. Verificación de 'df_annotated_invertidos' ---
cat("📊 Verificando 'df_annotated_invertidos'...\n")
glimpse(df_annotated_invertidos)
genes_con_go_en_invertidos <- sum(!is.na(df_annotated_invertidos$GO_IDs) & df_annotated_invertidos$GO_IDs != "")
cat(" -> Total de genes con anotación GO en invertidos:", genes_con_go_en_invertidos, "\n\n")

# --- 3. Verificación de 'genes_up_LPi' ---
cat("📊 Verificando 'genes_up_LPi'...\n")
cat(" -> Formato de los primeros IDs:", head(genes_up_LPi, 3), "\n")
cat(" -> Número total de genes en la lista 'up LPi':", length(genes_up_LPi), "\n")

# --- 4. Verificación CRÍTICA: Cobertura del universo ---
# Creamos el universo de genes que SÍ tienen anotación
universe <- unique(df_annotated_clusters$GeneID[!is.na(df_annotated_clusters$GO_IDs) & df_annotated_clusters$GO_IDs != ""])
cat(" -> Tamaño del universo (genes con GO):", length(universe), "\n")

# Comprobamos cuántos de nuestros genes de interés están en el universo
genes_up_LPi_en_universo <- sum(genes_up_LPi %in% universe)
cat(" -> De los", length(genes_up_LPi), "genes 'up LPi',", genes_up_LPi_en_universo, "están en el universo y pueden ser analizados.\n\n")

cat("✅ Chequeo de sanidad completado. Los datos parecen estar en orden.\n")
cat("=================================================================\n\n")

# --- El resto de tu script de enriquecimiento comenzaría aquí ---

##########################COTINUACION DEL CODIGO
####Importante descargar la siguiente libreria 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.At.thaliana.eg.db")

######################

# scripts/09_enrich_go_kegg.R (Versión Corregida)
# Objetivo: Enriquecer GO en subconjuntos usando anotaciones personalizadas.

library(tidyverse)
library(here)
library(clusterProfiler)

# --- 1. Cargar y Preparar Datos de Anotación Personalizados ---

# Cargar el archivo con TODAS las anotaciones GO para P. vulgaris
all_go_annotations <- read_csv(here("datos", "pvulgaris_all_go_annotations.csv"))

# Crear el archivo TERM2GENE: relaciona cada GO ID con su gen
# Esta es la base de datos para el enriquecimiento
# Crear el archivo TERM2GENE: relaciona cada GO ID con su gen
# Esta es la base de datos para el enriquecimiento
term2gene <- all_go_annotations %>%
  # Explicitamente le decimos a R que use la función select() del paquete dplyr
  dplyr::select(term = `GO ID`, gene = `Transcript Name`) %>%
  dplyr::filter(!is.na(term) & !is.na(gene)) %>%
  dplyr::distinct() # Nos aseguramos de que no haya filas duplicadas

cat(" -> Se creó el mapa TERM2GENE con", nrow(term2gene), "asociaciones.\n")

# (Opcional pero recomendado) Crear archivo TERM2NAME para tener descripciones legibles
# NOTA: Para esto necesitarás la base de datos de GO. La instalamos si no está.
if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db")
}
library(GO.db)

# Aquí especificamos que queremos usar la función select() del paquete dplyr
term2name <- as.data.frame(AnnotationDbi::select(GO.db, keys(GO.db), "TERM")) %>%
  dplyr::select(term = GOID, name = TERM)

cat(" -> Se cargó el mapa TERM2NAME con", nrow(term2name), "descripciones de GO.\n\n")


# --- 2. Cargar tus listas de genes de interés ---
df_annotated_clusters <- read_csv(here("salidas_data", "annotated_clusters.csv"))
venn_list <- readRDS(here("salidas_data", "venn_lists.rds"))
genes_up_LPi <- venn_list$LPi
genes_invertidos <- df_annotated_invertidos$GeneID


# --- 3. Realizar el Análisis de Enriquecimiento con enricher() ---

# Definir el universo: TODOS los genes que tienen al menos una anotación GO
universe <- unique(term2gene$gene)
cat(" -> Tamaño del universo corregido:", length(universe), "genes.\n")

# Enriquecimiento para los genes invertidos
enrich_invertidos <- enricher(
  gene = genes_invertidos,
  universe = universe,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# Enriquecimiento para los genes 'up LPi'
enrich_up_LPi <- enricher(
  gene = genes_up_LPi,
  universe = universe,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# --- 4. Visualizar y Guardar Resultados ---
cat("\n📊 Resultados de enriquecimiento para genes invertidos:\n")
print(as.data.frame(enrich_invertidos))

# Crear un dotplot si hay resultados significativos
if (nrow(as.data.frame(enrich_invertidos)) > 0) {
  dotplot(enrich_invertidos, showCategory = 15) + ggtitle("GO Enrichment - Genes Invertidos")
  ggsave(here("salidas", "dotplot_invertidos.png"))
} else {
  cat(" -> No se encontraron términos enriquecidos para los genes invertidos.\n")
}

write_csv(as.data.frame(enrich_invertidos), here("salidas_data", "enrich_go_invertidos.csv"))
cat("✅ Proceso de enriquecimiento completado.\n")
#################################
##REVISION DE LOS DATOS
# Revisa los resultados de tu lista de genes más grande
cat("\n📊 Resultados de enriquecimiento para genes 'up LPi':\n")
print(as.data.frame(enrich_up_LPi))

# Y visualízalos si hay resultados
if (nrow(as.data.frame(enrich_up_LPi)) > 0) {
  dotplot(enrich_up_LPi, showCategory = 15) + ggtitle("GO Enrichment - Genes Up LPi")
  ggsave(here("salidas", "dotplot_up_LPi.png"))
} else {
  cat(" -> No se encontraron términos enriquecidos para los genes 'up LPi'.\n")
}

############################
# --- 5. Enriquecimiento por Cluster (EL MÉTODO CORRECTO) ---

# Creamos una lista para guardar los resultados de cada cluster
enrich_clusters <- list()

# Hacemos un bucle que recorre cada número de cluster (del 1 al 5)
for (clust_num in unique(df_annotated_clusters$cluster_h)) {
  
  cat("\nAnalizando enriquecimiento para el Cluster", clust_num, "...\n")
  
  # Extraemos los genes que pertenecen a este cluster
  genes_del_cluster <- df_annotated_clusters %>%
    filter(cluster_h == clust_num) %>%
    pull(GeneID)
  
  # Realizamos el enriquecimiento para la lista de genes de este cluster
  resultado_enrich <- enricher(
    gene = genes_del_cluster,
    universe = universe,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  # Guardamos el resultado en nuestra lista
  enrich_clusters[[paste0("cluster_", clust_num)]] <- resultado_enrich
}

# Ahora puedes ver los resultados de un cluster específico, por ejemplo, el cluster 2
cat("\n📊 Resultados para el Cluster 2:\n")
print(as.data.frame(enrich_clusters$cluster_2))

#### hasta este punto el script sirve perfectamente 

# Imprime la tabla de resultados para el Cluster 3
print(as.data.frame(enrich_clusters$cluster_3))

# Crea un gráfico si hay resultados
if (nrow(as.data.frame(enrich_clusters$cluster_3)) > 0) {
  dotplot(enrich_clusters$cluster_3) + ggtitle("GO Enrichment - Cluster 3")
}


###################################asta este punto todo funciona correctamente 

# --- 6. Guardar Tablas y Gráficos para CADA Cluster ---

cat("\n💾 Iniciando el guardado de resultados por cluster...\n")

# Hacemos un bucle que recorre cada resultado en la lista 'enrich_clusters'
# names(enrich_clusters) nos da los nombres: "cluster_1", "cluster_2", etc.
for (nombre_cluster in names(enrich_clusters)) {
  
  # Extraemos el resultado de enriquecimiento para el cluster actual
  resultado_actual <- enrich_clusters[[nombre_cluster]]
  
  # Convertimos el resultado a una tabla (data frame)
  tabla_resultados <- as.data.frame(resultado_actual)
  
  cat(" -> Procesando", nombre_cluster, "...\n")
  
  # Verificamos si la tabla tiene resultados antes de intentar guardarla
  if (nrow(tabla_resultados) > 0) {
    
    # 1. GUARDAR LA TABLA COMO ARCHIVO CSV
    # ------------------------------------
    # Creamos un nombre de archivo dinámico, ej: "enrich_go_cluster_1.csv"
    nombre_archivo_csv <- here("salidas_data", paste0("enrich_go_", nombre_cluster, ".csv"))
    write_csv(tabla_resultados, nombre_archivo_csv)
    cat("    -> Tabla guardada en:", nombre_archivo_csv, "\n")
    
    # 2. CREAR Y GUARDAR EL GRÁFICO DOTPLOT
    # ------------------------------------
    # Creamos un nombre de archivo dinámico, ej: "dotplot_cluster_1.png"
    nombre_archivo_png <- here("salidas", paste0("dotplot_", nombre_cluster, ".png"))
    
    # Creamos el gráfico
    grafico <- dotplot(resultado_actual, showCategory = 15) + 
      ggtitle(paste("GO Enrichment -", gsub("_", " ", nombre_cluster))) # Título bonito
    
    # Guardamos el gráfico
    ggsave(nombre_archivo_png, plot = grafico, width = 8, height = 6)
    cat("    -> Gráfico guardado en:", nombre_archivo_png, "\n")
    
  } else {
    # Si no hay resultados, simplemente lo informamos
    cat("    -> Sin resultados de enriquecimiento significativos para este cluster.\n")
  }
}

cat("\n✅ ¡Proceso de guardado completado para todos los clusters!\n")

#######################asta este punto todo funciona correctamente 
# --- 7. (Opcional) Crear una Tabla Maestra con TODOS los resultados ---

# Usamos la función `bind_rows` para unir todas las tablas de la lista en una sola.
# El argumento .id = "cluster" crea una nueva columna que nos dice de qué cluster vino cada fila.

# Usamos map() para convertir cada enrichResult en un data.frame PRIMERO.
# Luego, el resultado (una lista de data.frames) se pasa a bind_rows().


tabla_maestra_clusters <- purrr::map(enrich_clusters, as.data.frame) %>%
  bind_rows(.id = "cluster")

# Verificamos cómo se ve la tabla maestra
cat("\n📊 Vista previa de la tabla maestra combinada:\n")
glimpse(tabla_maestra_clusters)

# Guardamos esta tabla maestra en un solo archivo CSV
write_csv(tabla_maestra_clusters, here("salidas_data", "enrich_go_TODOS_LOS_CLUSTERS.csv"))

cat("\n✅ Tabla maestra guardada exitosamente.\n")


############################3 hasta este punto sirve el el script 



######################## a partir de aqui 1-11-2025

#install.packages("zoo") Es necesario para que sirva este script 

# scripts/10_cruce_anotaciones.R
# Objetivo: Limpiar la tabla de genes simbióticos (S2) y cruzarla con nuestros resultados de enriquecimiento.
#library(janitor) # Muy útil para limpiar nombres de columnas

# --- 1. Cargar los datos que vamos a usar ---

# Tus resultados de enriquecimiento de todos los clusters
tabla_maestra_clusters <- read_csv(here("salidas_data", "enrich_go_TODOS_LOS_CLUSTERS.csv")) %>%
  as_tibble() 

tab_s2_cruda <- read_csv(here("datos", "tab_s2_early_signaling.csv"), col_names = FALSE)

#############################################
tab_s2_cruda <- read_csv(here("datos", "tab_s2_early_signaling.csv"), col_names = FALSE)
tab_s2_cruda <- as_tibble(tab_s2_cruda)  # Forzamos que sea un tibble


cat("Iniciando la limpieza de la Tabla S2 (versión CSV)...\n")

tab_s2_limpia <- read_csv(
  here("datos", "tab_s2_early_signaling.csv"),
  col_names = FALSE,
  skip = 2
) %>%
  row_to_names(row_number = 1) %>%
  clean_names() %>%
  mutate(across(everything(), as.character)) %>%
  mutate(
    across(contains("log2fc"), as.numeric),
    across(contains("padj"), as.numeric)
  )

cat("Tabla S2 limpiada correctamente.\n")
cat("Filas:", nrow(tab_s2_limpia), "| Columnas:", ncol(tab_s2_limpia), "\n")
print(head(tab_s2_limpia))


# Guardar la tabla limpia en la carpeta "salidas_data" (o donde quieras)
write_csv(
  tab_s2_limpia,
  here("salidas_data", "tab_s2_early_signaling_limpia.csv")
)

cat("Tabla S2 limpia guardada en: salidas_data/tab_s2_early_signaling_limpia.csv\n")

###################3hasta aqui todo funciono bien 
#library(readr)
#library(dplyr)      # ← Carga dplyr
#library(tidyr)      # ← Carga tidyr
#library(stringr)    # ← Carga stringr
#library(here)
#library(janitor)

cat("Iniciando cruce final...\n")

# --- CARGAR Y CONVERTIR ---
tab_s2_limpia <- read_csv(here("salidas_data", "tab_s2_early_signaling_limpia.csv")) %>%
  as_tibble()

tabla_maestra_clusters <- read_csv(here("salidas_data", "enrich_go_TODOS_LOS_CLUSTERS.csv")) %>%
  as_tibble()

# --- EXPANDIR geneID (con dplyr:: explícito) ---
clusters_expandido <- tabla_maestra_clusters %>%
  dplyr::select(cluster, ID, Description, RichFactor, FoldEnrichment, p.adjust, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(id_raw = geneID) %>%
  dplyr::mutate(
    id_clean = stringr::str_remove(id_raw, "\\.\\d+$") %>% stringr::str_trim()
  )

# --- LIMPIAR IDs en S2 ---
tab_s2_limpia <- tab_s2_limpia %>%
  dplyr::mutate(
    id_clean = stringr::str_remove(id, "\\.\\d+$") %>% stringr::str_trim()
  )

# --- CRUCE ---
cruce_final <- tab_s2_limpia %>%
  dplyr::left_join(clusters_expandido, by = "id_clean") %>%
  dplyr::select(-id_raw, -id_clean)

# --- GUARDAR ---
ruta_final <- here("salidas_data", "cruce_S2_con_clusters_ENRIQUECIDO.csv")
write_csv(cruce_final, ruta_final)

cat("¡CRUCE COMPLETADO!\n")
cat("Archivo guardado en:", ruta_final, "\n")
cat("Genes de S2 con GO asignado:", sum(!is.na(cruce_final$cluster)), "de", nrow(tab_s2_limpia), "\n")



########################### hasta aqui sirve el codigo

# --- SOLUCIÓN: Usa dplyr:: ---
genes_con_go <- cruce_final %>%
  dplyr::filter(!is.na(cluster)) %>%
  dplyr::select(name, id, cluster, Description, p.adjust) %>%
  dplyr::arrange(p.adjust)

# Ver resultado
print(genes_con_go)

# Guardar
write_csv(genes_con_go, here("salidas_data", "genes_S2_en_GO_enriquecido.csv"))

############################## hasta aqui todo bien 
#library(ggplot2)
#library(dplyr)

# Contar genes por cluster
resumen <- genes_con_go %>%
  count(cluster, name = "n_genes") %>%
  mutate(cluster = factor(cluster, levels = c("cluster_1", "cluster_5")))

# Gráfico
ggplot(resumen, aes(x = cluster, y = n_genes, fill = cluster)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_text(aes(label = n_genes), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("cluster_1" = "#1f78b4", "cluster_5" = "#33a02c")) +
  labs(
    title = "Genes simbiontes en clusters enriquecidos (GO)",
    subtitle = "16 de 194 genes (8.2%)",
    x = "Módulo de coexpresión",
    y = "Número de genes simbiontes",
    caption = "p.adjust < 0.05"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(here("salidas_data", "FIG_genes_simbiontes_en_clusters.png"), 
       width = 7, height = 6, dpi = 300, bg = "white")

#################### bien todo bien 
# Tabla bonita para paper
tabla_suplemento <- genes_con_go %>%
  distinct(name, id, .keep_all = TRUE) %>%
  select(name, id, cluster, Description, p.adjust) %>%
  arrange(cluster, p.adjust)

write_csv(tabla_suplemento, here("salidas_data", "Tabla_S2_genes_enriquecidos.csv"))


#####################################################TODOD BIEN HASTA

# Ejemplo: Log2FC de los 16 genes en cada contraste
expr_16 <- cruce_final %>%
  filter(!is.na(cluster)) %>%
  select(name, id, contains("log2fc")) %>%
  distinct(name, .keep_all = TRUE) %>%
  pivot_longer(cols = starts_with("log2fc"), names_to = "contraste", values_to = "log2fc") %>%
  mutate(contraste = case_when(
    contraste == "log2fc" ~ "PHR-L7OX vs EV (-Pi)",
    contraste == "log2fc_2" ~ "PHR-L7OX vs EV (+Pi)",
    contraste == "log2fc_3" ~ "PHR-L7RNAi vs EV (-Pi)",
    TRUE ~ contraste
  ))

# Heatmap
library(pheatmap)
pheatmap(
  expr_16 %>% pivot_wider(names_from = contraste, values_from = log2fc) %>% column_to_rownames("name"),
  scale = "row",
  main = "Expresión de 16 genes simbiontes bajo PHR-L7"
)










################################










# Ver algunos IDs limpios de cada tabla
cat("=== IDs en S2 (limpios) ===\n")
print(head(tab_s2_limpia$id_clean, 10))

cat("\n=== IDs en GO (limpios) ===\n")
print(head(clusters_expandido$id_clean, 10))

# ¿Hay intersección?
interseccion <- intersect(tab_s2_limpia$id_clean, clusters_expandido$id_clean)
cat("\nGenes en común:", length(interseccion), "\n")
if(length(interseccion) > 0) print(interseccion)


