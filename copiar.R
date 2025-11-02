source("instalar_paquetes.R")  # Ajusta la ruta si está en una subcarpeta

instalar_y_cargar(
  paquetes_cran = c("fs", "here", "tidyverse", "data.table", "pheatmap", "writexl", "readxl", "skimr", "ggvenn", "factoextra", "stats","corrplot"),
  paquetes_bioc = c("biomartr")
)
#######################################################
#opcion para forzar la instalación 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("Biostrings")







##############################################################
#############################################################
#importante para el analisis de las caracteristicas de un objeto
class(tabla_maestra_clusters)
str(tabla_maestra_clusters)        # Muestra la estructura: tipo de objeto, columnas, clases
glimpse(tabla_maestra_clusters)    # Similar a str(), pero más legible (requiere dplyr)
head(tabla_maestra_clusters)       # Primeras 6 filas
tail(tabla_maestra_clusters)       # Últimas 6 filas
summary(tabla_maestra_clusters)    # Estadísticas básicas por columna (solo para numéricas y factores)
names(tabla_maestra_clusters)           # Nombres de columnas
sapply(tabla_maestra_clusters, class)   # Clase de cada columna
table(tabla_maestra_clusters$columna)   # Frecuencia de valores (para factores o categóricas)
unique(tabla_maestra_clusters$columna)  # Valores únicos
sapply(tabla_maestra_clusters, typeof)
sapply(tabla_maestra_clusters, is.numeric)
glimpse(tabla_maestra_clusters)
sum(is.na(tabla_maestra_clusters))              # Total de NAs
colSums(is.na(tabla_maestra_clusters))          # NAs por columna
dim(tabla_maestra_clusters)           # Número de filas y columnas
nrow(tabla_maestra_clusters)          # Solo filas
ncol(tabla_maestra_clusters)          # Solo columnas
skim(tabla_maestra_clusters)          # Resumen detallado por columna
introduce(tabla_maestra_clusters)     # Estadísticas generales
plot_intro(tabla_maestra_clusters)    # Visualización de estructura



class(tabla_suplemento)
str(tabla_suplemento)        # Muestra la estructura: tipo de objeto, columnas, clases
glimpse(tabla_suplemento)    # Similar a str(), pero más legible (requiere dplyr)
head(tabla_suplemento)       # Primeras 6 filas
tail(tabla_suplemento)       # Últimas 6 filas
summary(tabla_suplemento)    # Estadísticas básicas por columna (solo para numéricas y factores)
names(tabla_suplemento)           # Nombres de columnas
sapply(tabla_suplemento, class)   # Clase de cada columna
table(tabla_suplemento$columna)   # Frecuencia de valores (para factores o categóricas)
unique(tabla_suplemento$columna)  # Valores únicos
sapply(tabla_suplemento, typeof)
sapply(tabla_suplemento, is.numeric)
glimpse(tabla_suplemento)
sum(is.na(tabla_suplemento))              # Total de NAs
colSums(is.na(tabla_suplemento))          # NAs por columna
dim(tabla_suplemento)           # Número de filas y columnas
nrow(tabla_suplemento)          # Solo filas
ncol(tabla_suplemento)          # Solo columnas
skim(tabla_suplemento)          # Resumen detallado por columna
introduce(tabla_suplemento)     # Estadísticas generales
plot_intro(tabla_suplemento)    # Visualización de estructura

############################################## todo sirve

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
