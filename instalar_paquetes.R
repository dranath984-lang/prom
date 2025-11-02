# instalar_paquetes.R

instalar_y_cargar <- function(paquetes_cran = c(), paquetes_bioc = c()) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in paquetes_cran) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("📦 Instalando paquete CRAN:", pkg, "...\n")
      tryCatch({
        install.packages(pkg)
        cat("✅ Instalado correctamente:", pkg, "\n")
      }, error = function(e) {
        cat("❌ Error al instalar", pkg, ":", e$message, "\n")
      })
    } else {
      cat("✔️ Ya está instalado:", pkg, "\n")
    }
    
    tryCatch({
      library(pkg, character.only = TRUE)
      cat("📚 Cargado:", pkg, "\n")
    }, error = function(e) {
      cat("⚠️ No se pudo cargar", pkg, ":", e$message, "\n")
    })
  }
  
  for (pkg in paquetes_bioc) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("🧬 Instalando paquete Bioconductor:", pkg, "...\n")
      tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
        cat("✅ Instalado correctamente:", pkg, "\n")
      }, error = function(e) {
        cat("❌ Error al instalar", pkg, ":", e$message, "\n")
      })
    } else {
      cat("✔️ Ya está instalado:", pkg, "\n")
    }
    
    tryCatch({
      library(pkg, character.only = TRUE)
      cat("📚 Cargado:", pkg, "\n")
    }, error = function(e) {
      cat("⚠️ No se pudo cargar", pkg, ":", e$message, "\n")
    })
  }
}
  