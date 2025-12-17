# ============================================================================
# ANÁLISIS GLASSO: RED DE REPARACIÓN DNA EN CÁNCER DE MAMA
# Dataset: TCGA-BRCA (solo muestras tumorales, solo Seq RNA)
# Enfoque: Análisis profundo con exploración exhaustiva de lambda
# ============================================================================

rm(list=ls())

# ----------------------------------------------------------------------------
# 0. INSTALACIÓN Y CARGA DE PAQUETES
# ----------------------------------------------------------------------------

cat("========== INSTALANDO PAQUETES ==========\n")

paquetes_cran <- c("glasso", "qgraph", "corrplot", "ggplot2", "igraph")
for(pkg in paquetes_cran) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

dependencias_bioc <- c("GenomicRanges", "SummarizedExperiment", 
                       "GenomicFeatures", "IRanges", "S4Vectors",
                       "GenomeInfoDb", "Biobase")

for(pkg in dependencias_bioc) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

if(!require("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
  library(TCGAbiolinks)
}

cat("\n✓ Paquetes cargados\n\n")

# ----------------------------------------------------------------------------
# 1. DESCARGAR DATOS DE TCGA
# ----------------------------------------------------------------------------

cat("========== DESCARGANDO DATOS DE TCGA-BRCA ==========\n")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

tryCatch({
  GDCdownload(query, method = "api")
  cat("✓ Descarga completada\n")
}, error = function(e) {
  cat("Nota: Los datos pueden estar descargados\n")
})

data <- GDCprepare(query)

cat("\n✓ Datos cargados\n")
cat("Dimensiones:", dim(data), "\n\n")

# ----------------------------------------------------------------------------
# 2. PREPROCESAMIENTO (SOLO TUMORES)
# ----------------------------------------------------------------------------

cat("========== PREPROCESAMIENTO ==========\n")

expr_matrix <- assay(data, "unstranded")
clinical <- colData(data)

# Normalización log2
expr_log <- log2(expr_matrix + 1)

# Filtrar genes de baja expresión
genes_expressed <- rowSums(expr_log > 1) >= (ncol(expr_log) * 0.5)
expr_filtered <- expr_log[genes_expressed, ]

cat("Genes después de filtrado:", nrow(expr_filtered), "\n")

# Seleccionar SOLO muestras tumorales
sample_type <- clinical$sample_type
tumor_samples <- which(sample_type == "Primary Tumor")

cat("Muestras tumorales:", length(tumor_samples), "\n")

expr_tumor <- expr_filtered[, tumor_samples]

cat("\n✓ Preprocesamiento completado\n")
cat("Matriz final:", dim(expr_tumor), "\n\n")

# ----------------------------------------------------------------------------
# 3. SELECCIÓN DE GENES EGO NETWORK STRING
# ----------------------------------------------------------------------------

cat("========== SELECCIONANDO GENES DE INTERÉS ==========\n")

# Nombre del archivo (como están en el mismo directorio, basta con el nombre)
archivo_genes <- "brca_ego_network_genes_simple.txt"

# Verificamos que exista para evitar errores
if(file.exists(archivo_genes)) {
  
  # Leemos el archivo línea por línea creando un vector de caracteres
  genes_of_interest <- readLines(archivo_genes)
  
  # (Opcional) Limpieza: Quitar espacios vacíos si hubieran quedado
  genes_of_interest <- trimws(genes_of_interest)
  genes_of_interest <- genes_of_interest[genes_of_interest != ""]
  
  cat(sprintf("✓ Se cargaron %d genes desde '%s'\n", length(genes_of_interest), archivo_genes))
  cat("Primeros 5 genes:", paste(head(genes_of_interest, 5), collapse=", "), "...\n")
  
} else {
  stop("ERROR CRÍTICO: No se encuentra el archivo 'brca_ego_network_genes_simple.txt'. \nAsegúrate de que el script de R y el txt estén en la misma carpeta o revisa tu working directory con getwd().")
}

cat("Genes solicitados:", length(genes_of_interest), "\n")

gene_info <- rowData(data)

genes_filtered_symbols <- gene_info$gene_name[genes_expressed]

genes_found <- intersect(genes_of_interest, genes_filtered_symbols)
cat("Genes encontrados:", length(genes_found), "\n")

genes_not_found <- setdiff(genes_of_interest, genes_found)
if(length(genes_not_found) > 0) {
  cat("\n⚠ Genes NO encontrados:", length(genes_not_found), "\n")
  if(length(genes_not_found) <= 10) {
    print(genes_not_found)
  }
}

if(length(genes_found) < 10) {
  stop("ERROR: Muy pocos genes encontrados")
}

# Extraer subset
expr_subset <- expr_tumor[genes_filtered_symbols %in% genes_found, ]
rownames(expr_subset) <- genes_filtered_symbols[genes_filtered_symbols %in% genes_found]

cat("\n✓ Genes en análisis:", nrow(expr_subset), "\n")
cat("Muestras:", ncol(expr_subset), "\n\n")

# Transponer para GLASSO (muestras en filas, genes en columnas)
expr_t <- t(expr_subset)

# Calcular correlación
cor_tumor <- cor(expr_t)

n_genes <- ncol(expr_t)
cat("Matriz de correlación:", dim(cor_tumor), "\n\n")

# ----------------------------------------------------------------------------
# 4.1 EXPLORACIÓN EXHAUSTIVA DE LAMBDA
# ----------------------------------------------------------------------------

cat("========== EXPLORACIÓN DE LAMBDA (CRITERIO BIOLÓGICO) ==========\n")


lambda_grid <- seq(0.05, 0.5, by=0.01) 
genes_interes <- c("BRCA1", "BRCA2")

cat("Probando", length(lambda_grid), "valores de lambda (", min(lambda_grid), "-", max(lambda_grid), ")\n\n")

glasso_results <- list()
metrics <- data.frame(
  lambda = lambda_grid,
  n_edges = numeric(length(lambda_grid)),
  density = numeric(length(lambda_grid)),
  brca1_conn = integer(length(lambda_grid)), 
  brca2_conn = integer(length(lambda_grid)) 
)

for(i in 1:length(lambda_grid)) {
  lam <- lambda_grid[i]
  
  # penalize.diagonal=FALSE es importante para no alterar la varianza propia
  res <- glasso(cor_tumor, rho=lam, penalize.diagonal=FALSE)
  glasso_results[[i]] <- res
  
  # Procesar matriz de precisión para métricas
  wi <- res$wi

  rownames(wi) <- rownames(cor_tumor)
  colnames(wi) <- colnames(cor_tumor)
  
  # Anular diagonal para conteos
  diag(wi) <- 0
  
  # Métricas Globales
  n_edges <- sum(wi[upper.tri(wi)] != 0)
  density <- n_edges / (n_genes * (n_genes - 1) / 2)
  
  # Métricas Específicas (Supervivencia de Genes)
  if("BRCA1" %in% rownames(wi)) {
    metrics$brca1_conn[i] <- sum(wi["BRCA1", ] != 0)
  }
  if("BRCA2" %in% rownames(wi)) {
    metrics$brca2_conn[i] <- sum(wi["BRCA2", ] != 0)
  }
  
  metrics$n_edges[i] <- n_edges
  metrics$density[i] <- density
  
  if(i %% 5 == 0) {
    cat(sprintf("Lambda=%.2f | Densidad: %.1f%% | BRCA1: %d enlaces | BRCA2: %d enlaces\n", 
                lam, density*100, metrics$brca1_conn[i], metrics$brca2_conn[i]))
  }
}

cat("\n✓ Grid de lambdas completado\n\n")

# ----------------------------------------------------------------------------
# 4.2 SELECCIÓN DE LAMBDA (Lógica de Supervivencia)
# ----------------------------------------------------------------------------

# Buscamos lambdas donde AMBOS genes tengan al menos 11 conexiones
indices_supervivencia <- which(metrics$brca1_conn > 11 & metrics$brca2_conn > 11)

if(length(indices_supervivencia) > 0) {
  best_idx <- max(indices_supervivencia)
  criterion_msg <- "Criterio: Máximo lambda con BRCA1 y BRCA2 conectados"
} else {
  best_idx <- max(indices_supervivencia)
  criterion_msg <- "ALERTA: No se encontró lambda que conecte ambos genes. Usando mínimo."
}

best_lambda <- metrics$lambda[best_idx]

cat("========== LAMBDA SELECCIONADO ==========\n")
cat(criterion_msg, "\n")
cat(sprintf("Lambda óptimo: %.2f\n", best_lambda))
cat(sprintf("Conexiones BRCA1: %d\n", metrics$brca1_conn[best_idx]))
cat(sprintf("Conexiones BRCA2: %d\n", metrics$brca2_conn[best_idx]))
cat(sprintf("Densidad global: %.1f%%\n\n", metrics$density[best_idx]*100))

# ----------------------------------------------------------------------------
# 5.1 VISUALIZACIÓN DE MÉTRICAS ACTUALIZADA
# ----------------------------------------------------------------------------
par(mfrow=c(2,2), mar=c(4,4,3,2))

# Plot 1: Conexiones de BRCA1 y BRCA2 vs Lambda (El gráfico clave)
plot(metrics$lambda, metrics$brca1_conn, type="b", pch=19, col="firebrick3",
     ylim=c(0, max(c(metrics$brca1_conn, metrics$brca2_conn)) + 2),
     xlab="Lambda", ylab="Grado (Nº Conexiones)",
     main="Supervivencia de Genes de Interés")
points(metrics$lambda, metrics$brca2_conn, type="b", pch=19, col="dodgerblue3")
abline(v=best_lambda, lty=2, lwd=2, col="black") # Marca el elegido
legend("topright", legend=c("BRCA1", "BRCA2"), col=c("firebrick3", "dodgerblue3"), 
       pch=19, bty="n", cex=0.8)

# Plot 2: Densidad Global vs Lambda
plot(metrics$lambda, metrics$density, type="b", pch=19, col="gray40",
     xlab="Lambda", ylab="Densidad Global",
     main="Densidad de la Red")
points(metrics$lambda[best_idx], metrics$density[best_idx], pch=19, col="red", cex=1.5)
text(best_lambda, metrics$density[best_idx], paste0("Sel: ", best_lambda), pos=3, col="red")

# Plot 3: Sparsity (Visualizando la limpieza)
plot(metrics$lambda, 1-metrics$density, type="b", pch=19, col="purple",
     xlab="Lambda", ylab="Sparsity", main="Sparsity vs Lambda")

# Plot 4: Resumen Textual Gráfico
plot.new()
text(0.5, 0.8, "RESULTADO FINAL", font=2, cex=1.5)
text(0.5, 0.6, paste("Lambda:", best_lambda), cex=1.2)
text(0.5, 0.5, paste("Enlaces BRCA1:", metrics$brca1_conn[best_idx]), cex=1.2, col="firebrick3")
text(0.5, 0.4, paste("Enlaces BRCA2:", metrics$brca2_conn[best_idx]), cex=1.2, col="dodgerblue3")
text(0.5, 0.2, "(Este valor equilibra limpieza\ncon retención biológica)", cex=0.8, font=3)


# ----------------------------------------------------------------------------
# 5.2. VISUALIZACIÓN DE REDES
# ----------------------------------------------------------------------------
# Muestro 3 estados: Uno laxo, el ÓPTIMO, y uno muy estricto
par(mfrow=c(1,3), mar=c(2,2,3,1))
lambdas_demo <- unique(sort(c(0.1, best_lambda, min(0.5, max(lambda_grid)))))

if(length(lambdas_demo) < 3) lambdas_demo <- c(lambdas_demo[1], best_lambda, lambdas_demo[2]) 

for(lam in lambdas_demo) {
  idx <- which(metrics$lambda == lam)[1]
  prec_temp <- glasso_results[[idx]]$wi
  colnames(prec_temp) <- colnames(cor_tumor)
  
  titulo <- ifelse(lam == best_lambda, paste("★ ÓPTIMO λ=", lam), paste("λ=", lam))
  color_borde <- ifelse(lam == best_lambda, "red", "black")
  
  qgraph(cov2cor(prec_temp),
         layout="spring",
         labels=colnames(cor_tumor),
         title=titulo,
         vsize=4, label.cex=0.5, cut=0.01, legend=FALSE)
  box(col=color_borde, lwd=ifelse(lam == best_lambda, 3, 1))
}

# ----------------------------------------------------------------------------
# 5.3 ASIGNACIÓN FINAL
# ----------------------------------------------------------------------------
glasso_final <- glasso_results[[best_idx]]
n_edges_final <- metrics$n_edges[best_idx]

# ----------------------------------------------------------------------------
# 6. DETECCIÓN DE COMUNIDADES
# ----------------------------------------------------------------------------

cat("\n--- Detectando comunidades ---\n")

adj <- (abs(glasso_final$wi) > 0.01) * 1
diag(adj) <- 0
g <- graph_from_adjacency_matrix(adj, mode="undirected")

communities <- cluster_louvain(g)
n_communities <- length(unique(membership(communities)))
modularity_val <- modularity(communities)

cat(sprintf("Comunidades detectadas: %d\n", n_communities))
cat(sprintf("Modularidad: %.3f\n\n", modularity_val))

# Mostrar composición
for(i in 1:n_communities) {
  genes_in_comm <- colnames(cor_tumor)[membership(communities) == i]
  cat(sprintf("Comunidad %d (%d genes): %s\n", 
              i, length(genes_in_comm),
              paste(genes_in_comm, collapse=", ")))
}

# Colores por comunidad
if(n_communities <= 8) {
  palette_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                      "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
} else {
  palette_colors <- rainbow(n_communities, alpha=0.8)
}

node_colors <- palette_colors[membership(communities)]
names(node_colors) <- colnames(cor_tumor)

# ----------------------------------------------------------------------------
# 7. VISUALIZACIÓN DE RED FINAL
# ----------------------------------------------------------------------------

cat("\n--- Generando visualización ---\n")

par(mfrow=c(1,2), mar=c(1,1,3,1))

prec_norm <- cov2cor(glasso_final$wi)
diag(prec_norm) <- 0 

# Red con comunidades
qgraph(glasso_final$wi,
       layout="spring",
       labels=colnames(cor_tumor),
       color=node_colors,
       title=sprintf("Red GLASSO (λ=%.2f)\n%d comunidades, %d conexiones", 
                     best_lambda, n_communities, n_edges_final),
       vsize=7,
       label.cex=0.7,
       cut=0.01,
       legend=FALSE,
       groups=as.factor(membership(communities)),
       diag = FALSE)

# 1. Librerías necesarias
if (!require("ggcorrplot")) install.packages("ggcorrplot")
if (!require("gridExtra")) install.packages("gridExtra") # Para ponerlos lado a lado
library(ggcorrplot)
library(gridExtra)
library(ggplot2)

# --- 2. PREPARACIÓN DE DATOS ---

# A. Matriz de Correlación (Datos Observados)
matriz_cor <- cor_tumor 

# B. Matriz de Precisión (Glasso -> Correlación Parcial)
# Invertimos signo, normalizamos y recuperamos nombres
matriz_prec <- -cov2cor(glasso_final$wi)
diag(matriz_prec) <- 1 
rownames(matriz_prec) <- rownames(cor_tumor)
colnames(matriz_prec) <- colnames(cor_tumor)

# --- 3. GENERACIÓN DE GRÁFICOS ---

library(ggcorrplot)
library(ggplot2)

# --- ESTRATEGIA PARA LOS LABELS ---
# Obtenemos los nombres de las variables (asumiendo que están en las columnas)
nombres_vars <- colnames(matriz_cor) 

# Definimos un "salto" para no mostrar todos. 
# Por ejemplo, 5 significa que mostrará el nombre 1, el 6, el 11, etc.
# AJUSTA ESTE NÚMERO según la cantidad de variables que tengas.
salto <- 5 
breaks_visibles <- nombres_vars[seq(1, length(nombres_vars), by = salto)]

# --- GRÁFICO 1: Correlación Estándar ---
plot_cor <- ggcorrplot(matriz_cor, 
                       method = "square", 
                       type = "lower", 
                       lab = FALSE,           
                       show.diag = FALSE,     
                       colors = c("navy", "white", "firebrick3"),
                       title = "Correlación (Pearson)",
                       tl.cex = 9
) + 
  # Aquí aplicamos el filtro de nombres a los ejes X e Y
  scale_x_discrete(breaks = breaks_visibles) +
  scale_y_discrete(breaks = breaks_visibles) +
  theme(legend.position = "bottom",
        # Opcional: rotar texto del eje X si aun así se tocan
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- GRÁFICO 2: Precisión Glasso ---
plot_prec <- ggcorrplot(matriz_prec, 
                        method = "square", 
                        type = "lower", 
                        lab = FALSE, 
                        show.diag = FALSE,
                        title = "Correlación Parcial (Glasso)",
                        tl.cex = 9
) + 
  scale_fill_gradient2(low = "navy", high = "firebrick3", mid = "white", 
                       midpoint = 0, limit = c(-0.25, 0.25), 
                       space = "Lab", name = "Corr",
                       oob = scales::squish) +
  # Aplicamos el mismo filtro de nombres al segundo gráfico
  scale_x_discrete(breaks = breaks_visibles) +
  scale_y_discrete(breaks = breaks_visibles) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(plot_cor)
print(plot_prec)

# ----------------------------------------------------------------------------
# 8. ANÁLISIS DE VECINDARIOS DE BRCA1/2
# ----------------------------------------------------------------------------

cat("\n========== VECINDARIOS DE BRCA1/2 ==========\n")

matriz_precision_final <- glasso_final$wi
rownames(matriz_precision_final) <- rownames(cor_tumor)
colnames(matriz_precision_final) <- colnames(cor_tumor)

analyze_neighborhood <- function(gene_name, precision_matrix, top_n=10) {
  if(!gene_name %in% colnames(precision_matrix)) {
    cat(sprintf("  %s no encontrado\n", gene_name))
    return(NULL)
  }
  
  gene_idx <- which(colnames(precision_matrix) == gene_name)
  connections <- precision_matrix[gene_idx, ]
  connections <- connections[names(connections) != gene_name]
  
  connections_sorted <- sort(abs(connections), decreasing=TRUE)
  top_connections <- connections_sorted[1:min(top_n, length(connections_sorted))]
  original_signs <- sign(connections[names(top_connections)])
  
  return(list(
    neighbors = names(top_connections),
    strengths = top_connections,
    signs = original_signs
  ))
}

# BRCA1
cat("\n--- BRCA1 ---\n")
brca1_neighbors <- analyze_neighborhood("BRCA1", matriz_precision_final, top_n=10)

if(!is.null(brca1_neighbors) && length(brca1_neighbors$neighbors) > 0) {
  cat("Top 10 vecinos directos:\n")
  for(i in 1:length(brca1_neighbors$neighbors)) {
    sign_str <- ifelse(brca1_neighbors$signs[i] > 0, "+", "-")
    cat(sprintf("  %2d. %-12s (fuerza: %.4f, dirección: %s)\n", 
                i, brca1_neighbors$neighbors[i], 
                brca1_neighbors$strengths[i], sign_str))
  }
} else {
  cat("  Sin vecinos significativos\n")
}

# BRCA2
cat("\n--- BRCA2 ---\n")
brca2_neighbors <- analyze_neighborhood("BRCA2", matriz_precision_final, top_n=10)

if(!is.null(brca2_neighbors) && length(brca2_neighbors$neighbors) > 0) {
  cat("Top 10 vecinos directos:\n")
  for(i in 1:length(brca2_neighbors$neighbors)) {
    sign_str <- ifelse(brca2_neighbors$signs[i] > 0, "+", "-")
    cat(sprintf("  %2d. %-12s (fuerza: %.4f, dirección: %s)\n", 
                i, brca2_neighbors$neighbors[i], 
                brca2_neighbors$strengths[i], sign_str))
  }
} else {
  cat("  Sin vecinos significativos\n")
}

# Overlap entre BRCA1 y BRCA2
if(!is.null(brca1_neighbors) && !is.null(brca2_neighbors)) {
  common_neighbors <- intersect(brca1_neighbors$neighbors, brca2_neighbors$neighbors)
  if(length(common_neighbors) > 0) {
    cat("\nVecinos compartidos por BRCA1 y BRCA2:\n")
    cat(" ", paste(common_neighbors, collapse=", "), "\n")
  }
}

# ----------------------------------------------------------------------------
# 9. ANÁLISIS DE ESTABILIDAD DE VECINDARIOS [No incluido en el análisis]
# ----------------------------------------------------------------------------

cat("\n========== ESTABILIDAD DE VECINDARIOS ==========\n")
cat("(Analizando cómo cambian con lambda)\n\n")

# Para BRCA1, ver vecinos en diferentes lambdas
lambda_stability <- c(0.15, 0.25, 0.35, best_lambda)

cat("--- Estabilidad de vecinos de BRCA1 ---\n")
brca1_stability <- list()

for(lam in lambda_stability) {
  idx <- which.min(abs(lambda_grid - lam))
  neighbors <- analyze_neighborhood("BRCA1", matriz_precision_final, top_n=5)
  
  if(!is.null(neighbors)) {
    brca1_stability[[as.character(lam)]] <- neighbors$neighbors
    cat(sprintf("\nλ=%.2f: %s\n", lam, paste(neighbors$neighbors, collapse=", ")))
  }
}

# Genes que aparecen en todos los lambdas (estables)
if(length(brca1_stability) > 0) {
  stable_neighbors <- Reduce(intersect, brca1_stability)
  cat("\nVecinos ESTABLES de BRCA1 (en todos los lambdas):\n")
  if(length(stable_neighbors) > 0) {
    cat(" ", paste(stable_neighbors, collapse=", "), "\n")
  } else {
    cat("  Ninguno (los vecindarios varían con lambda)\n")
  }
}

# ----------------------------------------------------------------------------
# 10. MÉTRICAS DE RED
# ----------------------------------------------------------------------------

cat("\n========== MÉTRICAS DE RED (TOPOLOGÍA GLASSO) ==========\n")

# 1. Crear objeto grafo desde tu matriz limpia
# Usamos 'matriz_precision_final' (la que tiene nombres)
# Convertimos a adyacencia: 1 si hay conexión, 0 si no (binario)
adj_matrix <- matriz_precision_final != 0
diag(adj_matrix) <- 0 # Sin auto-bucles
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# 2. Calcular Métricas Clave
grados <- degree(g)
betweenness <- betweenness(g, normalized = TRUE)

# 3. Crear Dataframe de Importancia
df_metrics <- data.frame(
  Gen = names(grados),
  Grado = grados,
  Intermediacion = betweenness
)

# Ordenar por Grado (Hubs)
top_hubs <- df_metrics[order(df_metrics$Grado, decreasing = TRUE), ][1:10, ]

# Ordenar por Intermediación (Puentes)
top_bridges <- df_metrics[order(df_metrics$Intermediacion, decreasing = TRUE), ][1:10, ]

cat("--- Top 10 Hubs (Mayor conectividad directa) ---\n")
print(top_hubs[, c("Gen", "Grado")])

cat("\n--- Top 10 Puentes (Mayor flujo de información) ---\n")
print(top_bridges[, c("Gen", "Intermediacion")])

# 4. Ver específicamente BRCA1
if("BRCA1" %in% df_metrics$Gen) {
  cat("\n--- Métricas BRCA1 ---\n")
  print(df_metrics[df_metrics$Gen == "BRCA1", ])
}

# 5. Ver específicamente BRCA2
if("BRCA2" %in% df_metrics$Gen) {
  cat("\n--- Métricas BRCA2 ---\n")
  print(df_metrics[df_metrics$Gen == "BRCA2", ])
}

# ----------------------------------------------------------------------------
# 11. VALIDACIÓN CON STRING
# ----------------------------------------------------------------------------

cat("========== VALIDACIÓN INTEGRAL: STRING vs TCGA vs GLASSO ==========\n")

# ----------------------------------------------------------------------------
# 11.1. CARGAR DATOS
# ----------------------------------------------------------------------------

file_string <- "brca_ego_network_data.csv"
if(!file.exists(file_string)) stop("Falta el archivo CSV de STRING.")

df_string <- read.csv(file_string, stringsAsFactors = FALSE)

library(dplyr)

df_edges <- df_string %>%
  filter(origin_gene != "SEED" & origin_gene != "From_Step1") %>% 
  select(GenA = origin_gene, GenB = gene_symbol, STRING_Score = score)

# ID único ordenado
df_edges$PairID <- apply(df_edges[, c("GenA", "GenB")], 1, function(x) paste(sort(x), collapse = "_"))

cat(sprintf("Conexiones de STRING cargadas: %d\n", nrow(df_edges)))


# ----------------------------------------------------------------------------
# 11.2. CRUZAR CON TUS DATOS (TCGA & GLASSO)
# ----------------------------------------------------------------------------

# A. Mapear Pearson (TCGA)
# ------------------------
get_pearson <- function(g1, g2, mat) {
  if(g1 %in% rownames(mat) && g2 %in% rownames(mat)) return(mat[g1, g2])
  else return(NA)
}
df_edges$Pearson_TCGA <- mapply(get_pearson, df_edges$GenA, df_edges$GenB, MoreArgs = list(mat = cor_tumor))

# B. Mapear Glasso (Lambda Óptimo)
# ------------------------
# Aseguramos tener la matriz parcial lista (Lambda 0.35 según tu elección)
lambda_best <- best_lambda
if(!exists("partial_corr")) {
  cat(sprintf("Recalculando Glasso (Lambda %.2f)...\n", lambda_best))
  glasso_tmp <- glasso(cor_tumor, rho = lambda_best, penalize.diagonal = FALSE)
  prec_tmp <- glasso_tmp$wi
  colnames(prec_tmp) <- colnames(cor_tumor); rownames(prec_tmp) <- rownames(cor_tumor)
  partial_corr <- -cov2cor(prec_tmp) # Invertimos signo
  diag(partial_corr) <- 1
}

get_glasso <- function(g1, g2, mat) {
  if(g1 %in% rownames(mat) && g2 %in% rownames(mat)) return(mat[g1, g2])
  else return(NA)
}
df_edges$Glasso_Partial <- mapply(get_glasso, df_edges$GenA, df_edges$GenB, MoreArgs = list(mat = partial_corr))


# ----------------------------------------------------------------------------
# 11.3. CREAR DF_FINAL Y TRUTH_SET
# ----------------------------------------------------------------------------

# Eliminamos NAs (pares que no están en tu matriz TCGA)
df_final <- na.omit(df_edges)

# Definimos detecciones
df_final$Detected_Glasso <- abs(df_final$Glasso_Partial) > 0.0001
df_final$Strong_Pearson <- abs(df_final$Pearson_TCGA) > 0.3

# --- AQUÍ CREAMOS TRUTH_SET ---
# Solo nos interesan las que STRING dice que son muy seguras (>700)
truth_set <- df_final %>% filter(STRING_Score >= 700)
total_truth <- nrow(truth_set)

cat(sprintf("Pares 'Verdad' (STRING Alto Score y presentes en TCGA): %d\n\n", total_truth))

# ----------------------------------------------------------------------------
# 11.4. ANÁLISIS DE EMBUDO (PÉRDIDAS Y GANANCIAS)
# ----------------------------------------------------------------------------

cat("========== DESGLOSE DE PÉRDIDAS ==========\n")

# PASO 1: Pearson > 0.3
loss_pearson <- truth_set %>% 
  filter(abs(Pearson_TCGA) < 0.3) %>% 
  arrange(desc(STRING_Score))

cat(sprintf("1. PÉRDIDA POR DATOS (Pearson < 0.3): %d (%.1f%%)\n", 
            nrow(loss_pearson), (nrow(loss_pearson)/total_truth)*100))
if(nrow(loss_pearson)>0) {
  cat("   Ejemplo: Proteínas que interactúan pero no se co-expresan.\n")
  print(head(loss_pearson[, c("GenA","GenB","STRING_Score","Pearson_TCGA")], 10))
}

# PASO 2: GLASSO
pearson_thresh <- 0.3
loss_glasso <- truth_set %>% 
  filter(abs(Pearson_TCGA) >= pearson_thresh & !Detected_Glasso) %>% 
  arrange(desc(abs(Pearson_TCGA)))

cat(sprintf("\n2. PÉRDIDA POR MODELO (Pearson > 0.3 pero Glasso borró): %d (%.1f%%)\n", 
            nrow(loss_glasso), (nrow(loss_glasso)/total_truth)*100))
if(nrow(loss_glasso)>0) {
  cat("   Correlaciones fuertes que Glasso consideró indirectas.\n")
  print(head(loss_glasso[, c("GenA","GenB","Pearson_TCGA","Glasso_Partial")], 10))
}

# PASO 3: Sobrevivientes
survivors <- truth_set %>% filter(Detected_Glasso & abs(Pearson_TCGA) >= 0.3)
cat(sprintf("\n3. RECUPERADOS (ÉXITO): %d (%.1f%%)\n", 
            nrow(survivors), (nrow(survivors)/total_truth)*100))


cat("\n========== RESCATES (GANANCIAS) ==========\n")
# PASO EXTRA: ¿Glasso vio algo que Pearson no?
gains <- truth_set %>% filter(abs(Pearson_TCGA) < 0.3 & Detected_Glasso)

if(nrow(gains) > 0) {
  cat(sprintf("¡ATENCIÓN! Glasso rescató %d conexiones con Pearson bajo:\n", nrow(gains)))
  print(gains[, c("GenA","GenB","STRING_Score","Pearson_TCGA","Glasso_Partial")])
} else {
  cat("No hubo rescates (Glasso no agregó conexiones donde Pearson era débil).\n")
}

# ----------------------------------------------------------------------------
# RESUMEN FINAL DEL EMBUDO
# ----------------------------------------------------------------------------
survivors <- truth_set %>% filter(Detected_Glasso & abs(Pearson_TCGA) >= pearson_thresh)

cat("========== RESUMEN DEL EMBUDO ==========\n")
cat(sprintf("1. Input STRING (High Conf): %d\n", total_truth))
cat(sprintf("2. Pasan filtro de Datos (Pearson): %d (perdimos %d)\n", total_truth - nrow(loss_pearson), nrow(loss_pearson)))
cat(sprintf("3. Pasan filtro de Modelo (Glasso): %d (perdimos %d en este paso)\n", nrow(survivors), nrow(loss_glasso)))
cat(sprintf("4. Tasa de recuperación final: %.1f%%\n", (nrow(survivors)/total_truth)*100))

# ----------------------------------------------------------------------------
# 11.5. LISTA DE SOBREVIVIENTES
# ----------------------------------------------------------------------------

cat("\n========== DETALLE DE LAS CONEXIONES VALIDADAS ==========\n")

if(nrow(survivors) > 0) {
  
  # Ordenamos por fuerza de interacción DIRECTA (Glasso Absoluto)
  survivors_sorted <- survivors %>% 
    arrange(desc(abs(Glasso_Partial))) %>%
    dplyr::select(GenA, GenB, STRING_Score, Pearson_TCGA, Glasso_Partial)
  
  # Mostramos el TOP 30 en consola para inspección rápida
  cat("Top 30 conexiones validadas (Ordenadas por fuerza en Glasso):\n")
  print(head(survivors_sorted, 30), row.names = FALSE)
  
  # Guardamos todo en un CSV para que lo puedas poner en un Anexo de tu tesis
  archivo_salida <- "validcion_biologica_survivors.csv"
  write.csv(survivors_sorted, archivo_salida, row.names = FALSE)
  
  cat(sprintf("\n✓ Se guardó la lista completa de %d conexiones en '%s'\n", 
              nrow(survivors), archivo_salida))
  
  # --- ANÁLISIS DE SENTIDO BIOLÓGICO RÁPIDO ---
  # Chequeo automático de parejas famosas
  parejas_famosas <- c("BRCA1_BARD1", "BRCA2_PALB2", "BRCA1_BRIP1", "BRCA1_MYC", "TP53_MDM2")
  
  cat("\n--- Chequeo de 'Celebridades' ---\n")
  for(par in parejas_famosas) {
    genes <- unlist(strsplit(par, "_"))
    # Buscamos en ambos sentidos A-B o B-A
    match <- survivors_sorted %>% 
      filter((GenA == genes[1] & GenB == genes[2]) | (GenA == genes[2] & GenB == genes[1]))
    
    if(nrow(match) > 0) {
      cat(sprintf("✅ %s DETECTADO (Glasso: %.4f | STRING: %d)\n", 
                  par, match$Glasso_Partial, match$STRING_Score))
    } else {
      cat(sprintf("❌ %s NO DETECTADO en el set final\n", par))
    }
  }
  
} else {
  cat("⚠️ No hubo sobrevivientes con los umbrales actuales.\n")
}


cat("========== ANÁLISIS DE NOVEDAD: LO QUE GLASSO VIO Y STRING NO ==========\n")

# 1. APLANAR LA MATRIZ GLASSO COMPLETA
# ----------------------------------------------------
# Convertimos la matriz n x n en una lista de pares
genes <- rownames(partial_corr)
pares_idx <- which(upper.tri(partial_corr), arr.ind = TRUE)

df_glasso_all <- data.frame(
  GenA = genes[pares_idx[,1]],
  GenB = genes[pares_idx[,2]],
  Glasso_Partial = partial_corr[pares_idx]
)

# Filtramos solo lo que Glasso detectó (distinto de 0)
# Usamos un umbral mínimo técnico, o si quieres ser estricto, > 0.01
df_glasso_detected <- df_glasso_all %>% 
  filter(abs(Glasso_Partial) > 0.0001) %>%
  mutate(PairID = apply(cbind(GenA, GenB), 1, function(x) paste(sort(x), collapse="_")))

cat(sprintf("Total de conexiones en la red Glasso (Lambda=0.35): %d\n", nrow(df_glasso_detected)))


# 2. CARGAR LO CONOCIDO (STRING)
# ----------------------------------------------------
# Usamos df_edges cargado anteriormente 
known_interactions <- df_edges %>% 
  filter(STRING_Score >= 700) %>%
  dplyr::select(PairID, STRING_Score)


# 3. ENCONTRAR LAS "NUEVAS" (ANTI-JOIN)
# ----------------------------------------------------
# Buscamos pares en Glasso que NO estén en 'known_interactions'
novel_connections <- df_glasso_detected %>%
  anti_join(known_interactions, by = "PairID") %>%
  arrange(desc(abs(Glasso_Partial)))

cat(sprintf("Conexiones conocidas (STRING High Conf): %d\n", nrow(df_glasso_detected) - nrow(novel_connections)))
cat(sprintf("Conexiones 'Novedosas' (Solo Glasso): %d\n\n", nrow(novel_connections)))


# 4. ANÁLISIS DE LAS NOVEDADES
# ----------------------------------------------------

cat("--- Top 150 Conexiones 'Novedosas' (Más fuertes en Glasso) ---\n")
cat("Interpretación: Pares que se co-regulan fuertemente en el tumor, \npero no son interacciones proteicas directas clásicas.\n\n")

print(head(novel_connections[, c("GenA", "GenB", "Glasso_Partial")], 150))


# 5. ¿SON REALMENTE NUEVAS O SOLO STRING DE BAJA CONFIANZA?
# ----------------------------------------------------
# Cruzamos con la base STRING completa (incluso score bajo) para ver si
# aparecen con score 400, 200, etc.

df_string_full <- df_string %>%
  mutate(PairID = apply(cbind(origin_gene, gene_symbol), 1, function(x) paste(sort(x), collapse="_"))) %>%
  dplyr::select(PairID, STRING_Full_Score = score)

novel_checked <- novel_connections %>%
  left_join(df_string_full, by = "PairID") %>%
  mutate(Status = ifelse(is.na(STRING_Full_Score), "Desconocido en STRING", 
                         paste("STRING Bajo Score:", STRING_Full_Score)))

cat("\n--- Detalle de las Novedades (Check de Score Bajo) ---\n")
print(head(novel_checked[, c("GenA", "GenB", "Glasso_Partial", "Status")], 15))

# 6. GUARDAR 
# ----------------------------------------------------
write.csv(novel_checked, "glasso_novel_connections.csv", row.names = FALSE)
cat("\n✓ Lista de conexiones novedosas guardada en 'glasso_novel_connections.csv'\n")


# ----------------------------------------------------------------------------
# 12. RESUMEN FINAL
# ----------------------------------------------------------------------------

cat("\n\n========== RESUMEN FINAL ==========\n")

cat("\nDATOS:\n")
cat(sprintf("  - Muestras tumorales: %d\n", nrow(expr_t)))
cat(sprintf("  - Genes analizados: %d\n", n_genes))

cat("\nRED GLASSO:\n")
cat(sprintf("  - Lambda óptimo: %.2f\n", best_lambda))
cat(sprintf("  - Conexiones: %d\n", n_edges_final))
cat(sprintf("  - Densidad: %.1f%%\n", metrics$density[best_idx]*100))
cat(sprintf("  - Comunidades: %d\n", n_communities))
cat(sprintf("  - Modularidad: %.3f\n", modularity_val))

cat("\nGENES CLAVE:\n")
cat(sprintf("  - Top hub: %s (%d conexiones)\n", top_hubs$Gen[1], top_hubs$Grado[1]))
cat(sprintf("  - Top broker: %s (betweenness: %.1f)\n", top_bridges$Gen[1], top_bridges$Intermediacion[1]))


cat("\n========== ANÁLISIS COMPLETADO ==========\n")