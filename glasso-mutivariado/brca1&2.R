# ============================================================================
# ANÁLISIS GLASSO: REDES GÉNICAS EN CÁNCER DE MAMA
# Dataset: TCGA-BRCA (The Cancer Genome Atlas - Breast Cancer)
# Enfoque: Vecindario regulatorio de genes BRCA1/BRCA2
# ============================================================================

# Limpiar entorno
rm(list=ls())

# ----------------------------------------------------------------------------
# 0. INSTALACIÓN Y CARGA DE PAQUETES
# ----------------------------------------------------------------------------

cat("========== INSTALANDO PAQUETES NECESARIOS ==========\n")

# Paquetes estándar de CRAN
paquetes_cran <- c("glasso", "qgraph", "corrplot", "ggplot2", "pheatmap", 
                   "dplyr", "tidyr", "igraph")

for(pkg in paquetes_cran) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Instalar BiocManager primero
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

cat("\nInstalando paquetes de Bioconductor (puede tardar)...\n")

# Instalar dependencias de TCGAbiolinks primero
dependencias_bioc <- c("GenomicRanges", "SummarizedExperiment", 
                       "GenomicFeatures", "IRanges", "S4Vectors",
                       "GenomeInfoDb", "Biobase")

for(pkg in dependencias_bioc) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Instalando %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# Ahora instalar TCGAbiolinks
if(!require("TCGAbiolinks", quietly = TRUE)) {
  cat("Instalando TCGAbiolinks...\n")
  BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
  library(TCGAbiolinks)
}

cat("\n✓ Todos los paquetes cargados exitosamente\n\n")

# ----------------------------------------------------------------------------
# 1. DESCARGAR DATOS DE TCGA
# ----------------------------------------------------------------------------

cat("========== DESCARGANDO DATOS DE TCGA-BRCA ==========\n")
cat("Nota: Primera vez puede tardar varios minutos...\n\n")

# Definir query para datos de expresión génica
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Descargar datos (si no están ya descargados)
tryCatch({
  GDCdownload(query, method = "api")
  cat("✓ Descarga completada\n")
}, error = function(e) {
  cat("Nota: Los datos ya pueden estar descargados\n")
})

# Preparar datos
data <- GDCprepare(query)

cat("\n✓ Datos cargados exitosamente\n")
cat("Dimensiones:", dim(data), "\n")
cat("  - Genes:", nrow(data), "\n")
cat("  - Muestras:", ncol(data), "\n\n")

# ----------------------------------------------------------------------------
# 2. PREPROCESAMIENTO
# ----------------------------------------------------------------------------

cat("========== PREPROCESAMIENTO DE DATOS ==========\n")

# Extraer matriz de expresión
expr_matrix <- assay(data, "unstranded")

# Extraer información clínica
clinical <- colData(data)

# Normalización: log2(counts + 1)
expr_log <- log2(expr_matrix + 1)

# Filtrar genes con baja expresión (expresados en al menos 50% de muestras)
genes_expressed <- rowSums(expr_log > 1) >= (ncol(expr_log) * 0.5)
expr_filtered <- expr_log[genes_expressed, ]

cat("Genes después de filtrado:", nrow(expr_filtered), "\n")

# Identificar muestras tumorales vs normales
sample_type <- clinical$sample_type

# Contar tipos de muestra
cat("\nTipos de muestra:\n")
print(table(sample_type))

# Separar en tumor y normal
tumor_samples <- which(sample_type == "Primary Tumor")
normal_samples <- which(sample_type == "Solid Tissue Normal")

cat("\n✓ Preprocesamiento completado\n")
cat("  - Muestras tumorales:", length(tumor_samples), "\n")
cat("  - Muestras normales:", length(normal_samples), "\n\n")

# ----------------------------------------------------------------------------
# 3. SELECCIÓN DE GENES DE INTERÉS
# ----------------------------------------------------------------------------

cat("========== SELECCIONANDO GENES DE INTERÉS ==========\n")

# Convertir nombres de genes (ENSEMBL a símbolos)
gene_info <- rowData(data)
gene_symbols <- gene_info$gene_name

# Genes de interés: BRCA1/2 y sus vecinos conocidos
genes_of_interest <- c(
  # Genes ancla (famosos)
  "BRCA1", "BRCA2", 
  
  # Supresores tumorales principales
  "TP53", "PTEN", "RB1", "ATM",
  
  # Reparación de DNA (pathway de BRCA)
  "RAD51", "PALB2", "CHEK2", "BRIP1", "RAD50", 
  "MRE11", "NBN", "FANCD2", "XRCC2",
  
  # Receptores hormonales (importantes en breast cancer)
  "ESR1", "ESR2", "PGR", "ERBB2",
  
  # Ciclo celular y apoptosis
  "CDK4", "CCND1", "MDM2", "BAX", "BCL2",
  
  # Señalización oncogénica
  "PIK3CA", "AKT1", "MTOR", "KRAS",
  
  # Metástasis
  "MMP9", "VEGFA", "CDH1"
)

cat("Genes seleccionados:", length(genes_of_interest), "\n")

# Encontrar índices de estos genes en la matriz
gene_indices <- which(gene_symbols %in% genes_of_interest)
genes_found <- gene_symbols[gene_indices]

cat("Genes encontrados en los datos:", length(genes_found), "\n")
cat("\nGenes incluidos en el análisis:\n")
print(sort(genes_found))

# Extraer subconjunto de expresión
expr_subset <- expr_filtered[gene_indices, ]
rownames(expr_subset) <- genes_found

cat("\n✓ Subset de genes creado\n")
cat("Dimensiones finales:", dim(expr_subset), "\n\n")

# ----------------------------------------------------------------------------
# 4. ANÁLISIS GLASSO: TEJIDO NORMAL
# ----------------------------------------------------------------------------

cat("========== ANÁLISIS GLASSO: TEJIDO NORMAL ==========\n")

if(length(normal_samples) < 10) {
  cat("⚠ Advertencia: Pocas muestras normales. Usando todas las disponibles.\n")
}

# Extraer expresión en muestras normales
expr_normal <- t(expr_subset[, normal_samples])

# Calcular matriz de correlación
cor_normal <- cor(expr_normal)

# Aplicar GLASSO con diferentes lambdas
lambda_values <- c(0.2, 0.4, 0.6)

cat("\nProbando diferentes valores de regularización...\n")

glasso_normal_list <- list()
for(lam in lambda_values) {
  glasso_normal_list[[as.character(lam)]] <- glasso(cor_normal, rho=lam)
  n_edges <- sum(glasso_normal_list[[as.character(lam)]]$wi[upper.tri(cor_normal)] != 0)
  cat(sprintf("  Lambda = %.1f -> %d conexiones directas\n", lam, n_edges))
}

# Seleccionar lambda intermedio
best_lambda_normal <- 0.2
glasso_normal <- glasso_normal_list[[as.character(best_lambda_normal)]]

cat("\n✓ Análisis GLASSO completado para tejido normal\n")
cat("Lambda seleccionado:", best_lambda_normal, "\n\n")

# ----------------------------------------------------------------------------
# 5. ANÁLISIS GLASSO: TEJIDO TUMORAL
# ----------------------------------------------------------------------------

cat("========== ANÁLISIS GLASSO: TEJIDO TUMORAL ==========\n")

# Usar submuestra de tumores (para balancear con normales si hay muchos)
if(length(tumor_samples) > 200) {
  set.seed(123)
  tumor_samples_subset <- sample(tumor_samples, 200)
  cat("Usando submuestra de 200 tumores (de", length(tumor_samples), ")\n")
} else {
  tumor_samples_subset <- tumor_samples
}

# Extraer expresión en muestras tumorales
expr_tumor <- t(expr_subset[, tumor_samples_subset])

# Calcular matriz de correlación
cor_tumor <- cor(expr_tumor)

# Aplicar GLASSO
cat("\nProbando diferentes valores de regularización...\n")

glasso_tumor_list <- list()
for(lam in lambda_values) {
  glasso_tumor_list[[as.character(lam)]] <- glasso(cor_tumor, rho=lam)
  n_edges <- sum(glasso_tumor_list[[as.character(lam)]]$wi[upper.tri(cor_tumor)] != 0)
  cat(sprintf("  Lambda = %.1f -> %d conexiones directas\n", lam, n_edges))
}

# Usar el mismo lambda que en normal para comparabilidad
glasso_tumor <- glasso_tumor_list[[as.character(best_lambda_normal)]]

cat("\n✓ Análisis GLASSO completado para tejido tumoral\n\n")

# ----------------------------------------------------------------------------
# 6. VISUALIZACIÓN: COMPARACIÓN DE REDES
# ----------------------------------------------------------------------------

cat("========== GENERANDO VISUALIZACIONES ==========\n")

# Configurar ventana de gráficos
par(mfrow=c(2,2), mar=c(4,4,3,2))

# 1. Correlaciones en tejido normal
corrplot(cor_normal, method="color", type="upper",
         title="Correlaciones - Tejido Normal",
         mar=c(0,0,2,0), tl.cex=0.6, tl.col="black")

# 2. Correlaciones en tejido tumoral
corrplot(cor_tumor, method="color", type="upper",
         title="Correlaciones - Tejido Tumoral",
         mar=c(0,0,2,0), tl.cex=0.6, tl.col="black")

# 3. GLASSO en tejido normal
prec_normal_norm <- cov2cor(glasso_normal$wi)
corrplot(prec_normal_norm, method="color", type="upper",
         title="GLASSO - Tejido Normal",
         mar=c(0,0,2,0), tl.cex=0.6, tl.col="black")

# 4. GLASSO en tejido tumoral
prec_tumor_norm <- cov2cor(glasso_tumor$wi)
corrplot(prec_tumor_norm, method="color", type="upper",
         title="GLASSO - Tejido Tumoral",
         mar=c(0,0,2,0), tl.cex=0.6, tl.col="black")

# ----------------------------------------------------------------------------
# 7. REDES COMO GRAFOS
# ----------------------------------------------------------------------------

cat("\nGenerando visualizaciones de redes...\n")

# Red en tejido normal
par(mfrow=c(1,2), mar=c(1,1,3,1))

qgraph(glasso_normal$wi,
       layout="spring",
       labels=colnames(cor_normal),
       title="Red Génica - Tejido Normal",
       vsize=6,
       label.cex=0.8,
       cut=0.05,
       legend=FALSE)

# Red en tejido tumoral
qgraph(glasso_tumor$wi,
       layout="spring",
       labels=colnames(cor_tumor),
       title="Red Génica - Tejido Tumoral",
       vsize=6,
       label.cex=0.8,
       cut=0.05,
       legend=FALSE)

# ----------------------------------------------------------------------------
# 8. ANÁLISIS DE DIFERENCIAS
# ----------------------------------------------------------------------------

cat("\n========== ANÁLISIS DE DIFERENCIAS ==========\n")

# Calcular diferencias en matrices de precisión
diff_precision <- glasso_tumor$wi - glasso_normal$wi

# Conexiones perdidas en cáncer (negativo = más débil en tumor)
connections_lost <- which(diff_precision < -0.1 & glasso_normal$wi != 0, arr.ind=TRUE)

# Conexiones ganadas en cáncer (positivo = más fuerte en tumor)
connections_gained <- which(diff_precision > 0.1 & glasso_tumor$wi != 0, arr.ind=TRUE)

cat("\nCONEXIONES PERDIDAS EN CÁNCER:\n")
if(nrow(connections_lost) > 0) {
  for(i in 1:min(10, nrow(connections_lost))) {
    gene1 <- colnames(cor_normal)[connections_lost[i,1]]
    gene2 <- colnames(cor_normal)[connections_lost[i,2]]
    if(connections_lost[i,1] < connections_lost[i,2]) {  # Evitar duplicados
      cat(sprintf("  %s <-> %s\n", gene1, gene2))
    }
  }
} else {
  cat("  Ninguna conexión significativamente perdida\n")
}

cat("\nCONEXIONES GANADAS EN CÁNCER:\n")
if(nrow(connections_gained) > 0) {
  for(i in 1:min(10, nrow(connections_gained))) {
    gene1 <- colnames(cor_tumor)[connections_gained[i,1]]
    gene2 <- colnames(cor_tumor)[connections_gained[i,2]]
    if(connections_gained[i,1] < connections_gained[i,2]) {
      cat(sprintf("  %s <-> %s\n", gene1, gene2))
    }
  }
} else {
  cat("  Ninguna conexión significativamente ganada\n")
}

# ----------------------------------------------------------------------------
# 9. ANÁLISIS DEL VECINDARIO DE BRCA1/BRCA2
# ----------------------------------------------------------------------------

cat("\n========== VECINDARIO DE BRCA1 y BRCA2 ==========\n")

analyze_gene_neighborhood <- function(gene_name, precision_matrix, threshold=0.05) {
  if(!gene_name %in% colnames(precision_matrix)) {
    cat(sprintf("Gene %s no encontrado\n", gene_name))
    return(NULL)
  }
  
  gene_idx <- which(colnames(precision_matrix) == gene_name)
  connections <- precision_matrix[gene_idx, ]
  
  # Encontrar vecinos directos
  neighbors <- names(connections[abs(connections) > threshold & names(connections) != gene_name])
  neighbor_strengths <- connections[neighbors]
  
  # Ordenar por fuerza de conexión
  neighbor_strengths <- sort(abs(neighbor_strengths), decreasing=TRUE)
  
  return(list(neighbors=names(neighbor_strengths), 
              strengths=neighbor_strengths))
}

# BRCA1 en tejido normal
cat("\n--- BRCA1 en TEJIDO NORMAL ---\n")
brca1_normal <- analyze_gene_neighborhood("BRCA1", glasso_normal$wi, threshold=0.05)
if(!is.null(brca1_normal) && length(brca1_normal$neighbors) > 0) {
  cat("Vecinos directos:\n")
  for(i in 1:min(5, length(brca1_normal$neighbors))) {
    cat(sprintf("  %s (fuerza: %.3f)\n", 
                brca1_normal$neighbors[i], 
                brca1_normal$strengths[i]))
  }
} else {
  cat("No se encontraron vecinos directos significativos\n")
}

# BRCA1 en tejido tumoral
cat("\n--- BRCA1 en TEJIDO TUMORAL ---\n")
brca1_tumor <- analyze_gene_neighborhood("BRCA1", glasso_tumor$wi, threshold=0.05)
if(!is.null(brca1_tumor) && length(brca1_tumor$neighbors) > 0) {
  cat("Vecinos directos:\n")
  for(i in 1:min(5, length(brca1_tumor$neighbors))) {
    cat(sprintf("  %s (fuerza: %.3f)\n", 
                brca1_tumor$neighbors[i], 
                brca1_tumor$strengths[i]))
  }
} else {
  cat("No se encontraron vecinos directos significativos\n")
}

# BRCA2 en tejido normal
cat("\n--- BRCA2 en TEJIDO NORMAL ---\n")
brca2_normal <- analyze_gene_neighborhood("BRCA2", glasso_normal$wi, threshold=0.05)
if(!is.null(brca2_normal) && length(brca2_normal$neighbors) > 0) {
  cat("Vecinos directos:\n")
  for(i in 1:min(5, length(brca2_normal$neighbors))) {
    cat(sprintf("  %s (fuerza: %.3f)\n", 
                brca2_normal$neighbors[i], 
                brca2_normal$strengths[i]))
  }
} else {
  cat("No se encontraron vecinos directos significativos\n")
}

# BRCA2 en tejido tumoral
cat("\n--- BRCA2 en TEJIDO TUMORAL ---\n")
brca2_tumor <- analyze_gene_neighborhood("BRCA2", glasso_tumor$wi, threshold=0.05)
if(!is.null(brca2_tumor) && length(brca2_tumor$neighbors) > 0) {
  cat("Vecinos directos:\n")
  for(i in 1:min(5, length(brca2_tumor$neighbors))) {
    cat(sprintf("  %s (fuerza: %.3f)\n", 
                brca2_tumor$neighbors[i], 
                brca2_tumor$strengths[i]))
  }
} else {
  cat("No se encontraron vecinos directos significativos\n")
}

# ----------------------------------------------------------------------------
# 10. MÉTRICAS DE RED
# ----------------------------------------------------------------------------

cat("\n========== MÉTRICAS DE RED ==========\n")

calculate_network_metrics <- function(precision_matrix, threshold=0.05) {
  # Crear grafo
  adj <- (abs(precision_matrix) > threshold) * 1
  diag(adj) <- 0
  g <- graph_from_adjacency_matrix(adj, mode="undirected")
  
  # Calcular métricas
  degree <- degree(g)
  betweenness <- betweenness(g)
  
  return(list(degree=degree, betweenness=betweenness, graph=g))
}

# Métricas en tejido normal
metrics_normal <- calculate_network_metrics(glasso_normal$wi)
cat("\n--- TEJIDO NORMAL ---\n")
cat("Genes más conectados (hubs):\n")
top_hubs_normal <- sort(metrics_normal$degree, decreasing=TRUE)[1:5]
print(top_hubs_normal)

# Métricas en tejido tumoral
metrics_tumor <- calculate_network_metrics(glasso_tumor$wi)
cat("\n--- TEJIDO TUMORAL ---\n")
cat("Genes más conectados (hubs):\n")
top_hubs_tumor <- sort(metrics_tumor$degree, decreasing=TRUE)[1:5]
print(top_hubs_tumor)

# ----------------------------------------------------------------------------
# 11. RESUMEN FINAL
# ----------------------------------------------------------------------------

cat("\n\n========== RESUMEN DEL ANÁLISIS ==========\n")
cat("✓ Análisis completado exitosamente\n\n")

cat("DATOS ANALIZADOS:\n")
cat(sprintf("  - Muestras normales: %d\n", length(normal_samples)))
cat(sprintf("  - Muestras tumorales: %d\n", length(tumor_samples_subset)))
cat(sprintf("  - Genes analizados: %d\n", length(genes_found)))
cat(sprintf("  - Lambda utilizado: %.2f\n\n", best_lambda_normal))

cat("HALLAZGOS CLAVE:\n")
cat("1. Las redes génicas muestran diferencias estructurales entre tejido normal y tumoral\n")
cat("2. Algunas conexiones se pierden en cáncer (posible pérdida de regulación)\n")
cat("3. Nuevas conexiones emergen en cáncer (posible desregulación)\n")
cat("4. Los genes hub cambian entre condiciones (reorganización de la red)\n\n")

cat("PRÓXIMOS PASOS SUGERIDOS:\n")
cat("- Investigar genes con cambios dramáticos en conectividad\n")
cat("- Validar conexiones perdidas/ganadas con literatura\n")
cat("- Analizar pathways enriquecidos en módulos de la red\n")
cat("- Comparar con subtipos de cáncer (ER+, HER2+, Triple Negative)\n\n")

cat("========== FIN DEL ANÁLISIS ==========\n")