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

# Genes de interés: BRCA1/2 y su top50 de vecinos a dos pasos
genes_of_interest <- c(
  "AASDH", "ABRAXAS1", "ABRAXAS2", "AGO2",  "ANAPC5",  "APEX1",  "ATM",  "ATP5F1D",  "AURKA",
  "BABAM1",  "BABAM2",  "BARD1",  "BCCIP",  "BRCA1",  "BRCA2",  "BRCC3",  "BRIP1",  "BRIX1",
  "BUB1B",  "CBX2",  "CBX4",  "CBX7",  "CBX8",  "CCNB1",  "CCND1",  "CCNT1",  "CCNT2",
  "CD4",  "CDC6",  "CDK13",  "CDK2",  "CDK7",  "CDK9",  "CDKN1A",  "CDKN1B",  "CIAO3",
  "CSTF1",  "DHX15",  "DHX32",  "DMC1",  "EIF3E",  "EIF6",  "EMSY",  "ERCC1",  "ERCC3",
  "FANCC",  "FANCD2",  "FANCE",  "FANCG",  "FANCI"
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

# Probar múltiples lambdas para encontrar el óptimo
lambda_values <- c(0.2, 0.3, 0.4, 0.5, 0.6)

cat("\n=== EXPLORANDO DIFERENTES VALORES DE LAMBDA ===\n")
cat("(Buscamos densidad de red ~20-30%)\n\n")

glasso_normal_list <- list()
for(lam in lambda_values) {
  glasso_normal_list[[as.character(lam)]] <- glasso(cor_normal, rho=lam)
  n_edges <- sum(glasso_normal_list[[as.character(lam)]]$wi[upper.tri(cor_normal)] != 0)
  n_genes <- ncol(cor_normal)
  densidad <- n_edges / (n_genes*(n_genes-1)/2)
  cat(sprintf("Lambda = %.1f -> %3d conexiones (densidad = %.1f%%)\n", 
              lam, n_edges, densidad*100))
}

# Seleccionar lambda óptimo (apuntamos a ~25-30% densidad)
best_lambda_normal <- 0.4
glasso_normal <- glasso_normal_list[[as.character(best_lambda_normal)]]

n_genes <- ncol(cor_normal)
n_edges_normal <- sum(glasso_normal$wi[upper.tri(cor_normal)] != 0)
densidad_normal <- n_edges_normal / (n_genes*(n_genes-1)/2)

cat("\n✓ Lambda seleccionado:", best_lambda_normal, "\n")
cat("  Genes analizados:", n_genes, "\n")
cat("  Conexiones encontradas:", n_edges_normal, "\n")
cat("  Densidad de red:", sprintf("%.1f%%", densidad_normal*100), "\n\n")

# Visualizar comparación de lambdas
par(mfrow=c(2,3), mar=c(2,2,3,1))
for(lam in lambda_values) {
  prec_norm <- cov2cor(glasso_normal_list[[as.character(lam)]]$wi)
  n_conn <- sum(glasso_normal_list[[as.character(lam)]]$wi[upper.tri(cor_normal)] != 0)
  qgraph(prec_norm,
         layout="spring",
         labels=colnames(cor_normal),
         title=sprintf("Lambda=%.1f (%d conex.)", lam, n_conn),
         vsize=4,
         label.cex=0.6,
         cut=0.01,
         legend=FALSE)
}

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

# Aplicar GLASSO con el mismo lambda para comparabilidad
glasso_tumor <- glasso(cor_tumor, rho=best_lambda_normal)

n_edges_tumor <- sum(glasso_tumor$wi[upper.tri(cor_tumor)] != 0)
densidad_tumor <- n_edges_tumor / (n_genes*(n_genes-1)/2)

cat("\n✓ Lambda utilizado:", best_lambda_normal, "\n")
cat("  Genes analizados:", n_genes, "\n")
cat("  Conexiones encontradas:", n_edges_tumor, "\n")
cat("  Densidad de red:", sprintf("%.1f%%", densidad_tumor*100), "\n\n")

# ----------------------------------------------------------------------------
# 6. DEFINIR GRUPOS FUNCIONALES PARA COLOREAR
# ----------------------------------------------------------------------------

# Definir grupos de genes
gene_groups <- list(
  ancla = c("BRCA1", "BRCA2"),
  supresores = c("TP53", "PTEN", "RB1", "ATM"),
  reparacion = c("RAD51", "PALB2", "CHEK2", "BRIP1", "RAD50", "MRE11", "NBN", "FANCD2", "XRCC2"),
  receptores = c("ESR1", "ESR2", "PGR", "ERBB2"),
  ciclo = c("CDK4", "CCND1", "MDM2", "BAX", "BCL2"),
  senalizacion = c("PIK3CA", "AKT1", "MTOR", "KRAS"),
  metastasis = c("MMP9", "VEGFA", "CDH1")
)

# Asignar colores
gene_colors <- rep("gray80", n_genes)
names(gene_colors) <- colnames(cor_normal)

color_map <- c(
  ancla = "#E41A1C",        # Rojo intenso
  supresores = "#377EB8",   # Azul
  reparacion = "#4DAF4A",   # Verde
  receptores = "#FF7F00",   # Naranja
  ciclo = "#984EA3",        # Púrpura
  senalizacion = "#F781BF", # Rosa
  metastasis = "#A65628"    # Marrón
)

for(grupo in names(gene_groups)) {
  genes_in_group <- intersect(gene_groups[[grupo]], colnames(cor_normal))
  gene_colors[genes_in_group] <- color_map[grupo]
}

# ----------------------------------------------------------------------------
# 7. VISUALIZACIÓN: COMPARACIÓN DE REDES CON COLORES
# ----------------------------------------------------------------------------

cat("========== GENERANDO VISUALIZACIONES ==========\n")

# Visualización mejorada: Correlaciones vs GLASSO
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
# 8. REDES COMO GRAFOS CON COLORES POR GRUPO FUNCIONAL
# ----------------------------------------------------------------------------

cat("\nGenerando visualizaciones de redes con grupos funcionales...\n")

# Red en tejido normal
par(mfrow=c(1,2), mar=c(1,1,3,1))

qgraph(glasso_normal$wi,
       layout="spring",
       labels=colnames(cor_normal),
       color=gene_colors,
       title=sprintf("Red Génica - Tejido Normal\n(%d conexiones)", n_edges_normal),
       vsize=7,
       label.cex=0.7,
       cut=0.01,
       legend=FALSE)

# Red en tejido tumoral
qgraph(glasso_tumor$wi,
       layout="spring",
       labels=colnames(cor_tumor),
       color=gene_colors,
       title=sprintf("Red Génica - Tejido Tumoral\n(%d conexiones)", n_edges_tumor),
       vsize=7,
       label.cex=0.7,
       cut=0.01,
       legend=FALSE)

# Agregar leyenda de colores
par(mfrow=c(1,1), mar=c(1,1,3,1))
plot.new()
legend("center", 
       legend=c("BRCA1/2 (Ancla)", "Supresores Tumorales", "Reparación DNA", 
                "Receptores Hormonales", "Ciclo Celular", "Señalización", "Metástasis"),
       fill=c(color_map["ancla"], color_map["supresores"], color_map["reparacion"],
              color_map["receptores"], color_map["ciclo"], color_map["senalizacion"],
              color_map["metastasis"]),
       title="Grupos Funcionales",
       cex=1.2,
       bty="n")

# ----------------------------------------------------------------------------
# 9. ANÁLISIS DE DIFERENCIAS
# ----------------------------------------------------------------------------

cat("\n========== ANÁLISIS DE DIFERENCIAS ENTRE NORMAL Y TUMOR ==========\n")

# Calcular diferencias en matrices de precisión
diff_precision <- glasso_tumor$wi - glasso_normal$wi

# Usar threshold más bajo (0.03 en lugar de 0.1)
threshold_diff <- 0.03

# Conexiones perdidas en cáncer (negativo = más débil en tumor)
connections_lost <- which(diff_precision < -threshold_diff & 
                            glasso_normal$wi != 0, arr.ind=TRUE)

# Conexiones ganadas en cáncer (positivo = más fuerte en tumor)
connections_gained <- which(diff_precision > threshold_diff & 
                              glasso_tumor$wi != 0, arr.ind=TRUE)

# Remover duplicados (matriz simétrica)
connections_lost <- connections_lost[connections_lost[,1] < connections_lost[,2], , drop=FALSE]
connections_gained <- connections_gained[connections_gained[,1] < connections_gained[,2], , drop=FALSE]

cat("\nCONEXIONES PERDIDAS EN CÁNCER (", nrow(connections_lost), " encontradas):\n", sep="")
if(nrow(connections_lost) > 0) {
  lost_df <- data.frame(
    Gen1 = colnames(cor_normal)[connections_lost[,1]],
    Gen2 = colnames(cor_normal)[connections_lost[,2]],
    Cambio = diff_precision[connections_lost]
  )
  lost_df <- lost_df[order(lost_df$Cambio), ]
  
  for(i in 1:min(10, nrow(lost_df))) {
    cat(sprintf("  %s <-> %s (cambio: %.3f)\n", 
                lost_df$Gen1[i], lost_df$Gen2[i], lost_df$Cambio[i]))
  }
  if(nrow(lost_df) > 10) cat(sprintf("  ... y %d más\n", nrow(lost_df)-10))
} else {
  cat("  Ninguna conexión significativamente perdida\n")
}

cat("\nCONEXIONES GANADAS EN CÁNCER (", nrow(connections_gained), " encontradas):\n", sep="")
if(nrow(connections_gained) > 0) {
  gained_df <- data.frame(
    Gen1 = colnames(cor_tumor)[connections_gained[,1]],
    Gen2 = colnames(cor_tumor)[connections_gained[,2]],
    Cambio = diff_precision[connections_gained]
  )
  gained_df <- gained_df[order(-gained_df$Cambio), ]
  
  for(i in 1:min(10, nrow(gained_df))) {
    cat(sprintf("  %s <-> %s (cambio: +%.3f)\n", 
                gained_df$Gen1[i], gained_df$Gen2[i], gained_df$Cambio[i]))
  }
  if(nrow(gained_df) > 10) cat(sprintf("  ... y %d más\n", nrow(gained_df)-10))
} else {
  cat("  Ninguna conexión significativamente ganada\n")
}

# Distribución de cambios
cat("\nDistribución de cambios en conexiones:\n")
changes <- diff_precision[upper.tri(diff_precision)]
print(summary(changes))

# ----------------------------------------------------------------------------
# 10. ANÁLISIS DEL VECINDARIO DE BRCA1/BRCA2 (MEJORADO)
# ----------------------------------------------------------------------------

cat("\n========== VECINDARIO DE BRCA1 y BRCA2 ==========\n")

# Función mejorada: muestra top N vecinos sin importar threshold
analyze_gene_neighborhood_v2 <- function(gene_name, precision_matrix, top_n=5) {
  if(!gene_name %in% colnames(precision_matrix)) {
    cat(sprintf("  Gen %s no encontrado en el análisis\n", gene_name))
    return(NULL)
  }
  
  gene_idx <- which(colnames(precision_matrix) == gene_name)
  connections <- precision_matrix[gene_idx, ]
  
  # Remover el gen mismo
  connections <- connections[names(connections) != gene_name]
  
  # Ordenar por valor absoluto
  connections_sorted <- sort(abs(connections), decreasing=TRUE)
  
  # Tomar top N
  top_connections <- connections_sorted[1:min(top_n, length(connections_sorted))]
  
  # Obtener signos originales
  original_signs <- sign(connections[names(top_connections)])
  
  return(list(
    neighbors = names(top_connections), 
    strengths = top_connections,
    signs = original_signs
  ))
}

# Análisis de BRCA1
cat("\n--- BRCA1 ---\n")

cat("\nEn TEJIDO NORMAL:\n")
brca1_normal <- analyze_gene_neighborhood_v2("BRCA1", glasso_normal$wi, top_n=5)
if(!is.null(brca1_normal) && length(brca1_normal$neighbors) > 0) {
  cat("Top 5 vecinos más cercanos:\n")
  for(i in 1:length(brca1_normal$neighbors)) {
    sign_str <- ifelse(brca1_normal$signs[i] > 0, "+", "-")
    cat(sprintf("  %d. %s (fuerza: %.4f, dirección: %s)\n", 
                i,
                brca1_normal$neighbors[i], 
                brca1_normal$strengths[i],
                sign_str))
  }
} else {
  cat("  Sin vecinos significativos (todas las conexiones fueron eliminadas por GLASSO)\n")
}

cat("\nEn TEJIDO TUMORAL:\n")
brca1_tumor <- analyze_gene_neighborhood_v2("BRCA1", glasso_tumor$wi, top_n=5)
if(!is.null(brca1_tumor) && length(brca1_tumor$neighbors) > 0) {
  cat("Top 5 vecinos más cercanos:\n")
  for(i in 1:length(brca1_tumor$neighbors)) {
    sign_str <- ifelse(brca1_tumor$signs[i] > 0, "+", "-")
    cat(sprintf("  %d. %s (fuerza: %.4f, dirección: %s)\n", 
                i,
                brca1_tumor$neighbors[i], 
                brca1_tumor$strengths[i],
                sign_str))
  }
} else {
  cat("  Sin vecinos significativos\n")
}

# Análisis de BRCA2
cat("\n--- BRCA2 ---\n")

cat("\nEn TEJIDO NORMAL:\n")
brca2_normal <- analyze_gene_neighborhood_v2("BRCA2", glasso_normal$wi, top_n=5)
if(!is.null(brca2_normal) && length(brca2_normal$neighbors) > 0) {
  cat("Top 5 vecinos más cercanos:\n")
  for(i in 1:length(brca2_normal$neighbors)) {
    sign_str <- ifelse(brca2_normal$signs[i] > 0, "+", "-")
    cat(sprintf("  %d. %s (fuerza: %.4f, dirección: %s)\n", 
                i,
                brca2_normal$neighbors[i], 
                brca2_normal$strengths[i],
                sign_str))
  }
} else {
  cat("  Sin vecinos significativos\n")
}

cat("\nEn TEJIDO TUMORAL:\n")
brca2_tumor <- analyze_gene_neighborhood_v2("BRCA2", glasso_tumor$wi, top_n=5)
if(!is.null(brca2_tumor) && length(brca2_tumor$neighbors) > 0) {
  cat("Top 5 vecinos más cercanos:\n")
  for(i in 1:length(brca2_tumor$neighbors)) {
    sign_str <- ifelse(brca2_tumor$signs[i] > 0, "+", "-")
    cat(sprintf("  %d. %s (fuerza: %.4f, dirección: %s)\n", 
                i,
                brca2_tumor$neighbors[i], 
                brca2_tumor$strengths[i],
                sign_str))
  }
} else {
  cat("  Sin vecinos significativos\n")
}

# ----------------------------------------------------------------------------
# 11. MÉTRICAS DE RED
# ----------------------------------------------------------------------------

cat("\n========== MÉTRICAS DE RED ==========\n")

calculate_network_metrics <- function(precision_matrix, threshold=0.01) {
  # Crear grafo
  adj <- (abs(precision_matrix) > threshold) * 1
  diag(adj) <- 0
  g <- graph_from_adjacency_matrix(adj, mode="undirected")
  
  # Calcular métricas
  deg <- degree(g)
  betw <- betweenness(g)
  
  return(list(degree=deg, betweenness=betw, graph=g))
}

# Métricas en tejido normal
metrics_normal <- calculate_network_metrics(glasso_normal$wi, threshold=0.01)
cat("\n--- TEJIDO NORMAL ---\n")
cat("Genes más conectados (hubs):\n")
top_hubs_normal <- sort(metrics_normal$degree, decreasing=TRUE)[1:5]
for(i in 1:length(top_hubs_normal)) {
  cat(sprintf("  %d. %s (%d conexiones)\n", i, names(top_hubs_normal)[i], top_hubs_normal[i]))
}

cat("\nGenes con mayor intermediación (brokers):\n")
top_betw_normal <- sort(metrics_normal$betweenness, decreasing=TRUE)[1:5]
for(i in 1:length(top_betw_normal)) {
  cat(sprintf("  %d. %s (betweenness: %.1f)\n", i, names(top_betw_normal)[i], top_betw_normal[i]))
}

# Métricas en tejido tumoral
metrics_tumor <- calculate_network_metrics(glasso_tumor$wi, threshold=0.01)
cat("\n--- TEJIDO TUMORAL ---\n")
cat("Genes más conectados (hubs):\n")
top_hubs_tumor <- sort(metrics_tumor$degree, decreasing=TRUE)[1:5]
for(i in 1:length(top_hubs_tumor)) {
  cat(sprintf("  %d. %s (%d conexiones)\n", i, names(top_hubs_tumor)[i], top_hubs_tumor[i]))
}

cat("\nGenes con mayor intermediación (brokers):\n")
top_betw_tumor <- sort(metrics_tumor$betweenness, decreasing=TRUE)[1:5]
for(i in 1:length(top_betw_tumor)) {
  cat(sprintf("  %d. %s (betweenness: %.1f)\n", i, names(top_betw_tumor)[i], top_betw_tumor[i]))
}

# ----------------------------------------------------------------------------
# 12. RESUMEN FINAL (DINÁMICO)
# ----------------------------------------------------------------------------

cat("\n\n========== RESUMEN DEL ANÁLISIS ==========\n")
cat("✓ Análisis completado exitosamente\n\n")

cat("DATOS ANALIZADOS:\n")
cat(sprintf("  - Muestras normales: %d\n", nrow(expr_normal)))
cat(sprintf("  - Muestras tumorales: %d\n", nrow(expr_tumor)))
cat(sprintf("  - Genes analizados: %d\n", n_genes))
cat(sprintf("  - Lambda utilizado: %.2f\n\n", best_lambda_normal))

cat("TOPOLOGÍA DE REDES:\n")
cat(sprintf("  - Conexiones en tejido normal: %d (densidad: %.1f%%)\n", 
            n_edges_normal, densidad_normal*100))
cat(sprintf("  - Conexiones en tejido tumoral: %d (densidad: %.1f%%)\n", 
            n_edges_tumor, densidad_tumor*100))
cat(sprintf("  - Cambio neto: %+d conexiones\n\n", n_edges_tumor - n_edges_normal))

cat("CAMBIOS ESTRUCTURALES:\n")
cat(sprintf("  - Conexiones perdidas en cáncer: %d\n", nrow(connections_lost)))
cat(sprintf("  - Conexiones ganadas en cáncer: %d\n", nrow(connections_gained)))

if(nrow(connections_lost) > 0 | nrow(connections_gained) > 0) {
  cat("\nINTERPRETACIÓN:\n")
  cat("Las redes génicas muestran reorganización estructural en cáncer:\n")
  if(nrow(connections_lost) > 0) {
    cat(sprintf("  - %d conexiones se debilitan/pierden (posible pérdida de regulación)\n", 
                nrow(connections_lost)))
  }
  if(nrow(connections_gained) > 0) {
    cat(sprintf("  - %d nuevas conexiones emergen (posible desregulación)\n", 
                nrow(connections_gained)))
  }
}

cat("\nGENES HUB (más conectados):\n")
cat("  Normal: ", paste(names(top_hubs_normal)[1:3], collapse=", "), "\n")
cat("  Tumor:  ", paste(names(top_hubs_tumor)[1:3], collapse=", "), "\n")

cat("\n\nPRÓXIMOS PASOS SUGERIDOS:\n")
cat("- Investigar genes con cambios dramáticos en conectividad\n")
cat("- Validar conexiones perdidas/ganadas con literatura\n")
cat("- Analizar pathways enriquecidos en módulos de la red\n")
cat("- Comparar con subtipos moleculares de cáncer de mama\n")
cat("- Evaluar asociación con supervivencia de pacientes\n\n")

cat("========== FIN DEL ANÁLISIS ==========\n")


# DEBUG

# Ejecutá esto para ver las correlaciones en tus datos
cat("\n=== TOP 15 CORRELACIONES CON BRCA1 (NORMAL) ===\n")
if("BRCA1" %in% colnames(cor_normal)) {
  brca1_cors <- cor_normal["BRCA1", ]
  brca1_cors <- brca1_cors[names(brca1_cors) != "BRCA1"]
  top15 <- sort(abs(brca1_cors), decreasing=TRUE)[1:15]
  for(i in 1:length(top15)) {
    cat(sprintf("%2d. %s (r=%.3f)\n", i, names(top15)[i], top15[i]))
  }
} else {
  cat("BRCA1 no encontrado\n")
}

cat("\n=== TOP 15 CORRELACIONES CON BRCA1 (TUMOR) ===\n")
if("BRCA1" %in% colnames(cor_tumor)) {
  brca1_cors_t <- cor_tumor["BRCA1", ]
  brca1_cors_t <- brca1_cors_t[names(brca1_cors_t) != "BRCA1"]
  top15_t <- sort(abs(brca1_cors_t), decreasing=TRUE)[1:15]
  for(i in 1:length(top15_t)) {
    cat(sprintf("%2d. %s (r=%.3f)\n", i, names(top15_t)[i], top15_t[i]))
  }
} else {
  cat("BRCA1 no encontrado\n")
}

# Ver cuáles de los genes de STRING están en tu dataset actual
cat("\n=== GENES DE STRING (top score) EN TU DATASET ===\n")
string_genes <- c("ATM", "TP53", "PALB2", "FANCD2", "MRE11", "BARD1", 
                  "BRIP1", "TOPBP1", "BABAM1", "BABAM2")
genes_presentes <- intersect(string_genes, colnames(cor_normal))
genes_ausentes <- setdiff(string_genes, colnames(cor_normal))

cat("Presentes en tu análisis:\n")
print(genes_presentes)
cat("\nAUSENTES de tu análisis:\n")
print(genes_ausentes)