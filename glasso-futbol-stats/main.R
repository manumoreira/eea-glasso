# ============================================================================
# GLASSO APLICADO A ANÁLISIS DE REDES SOCIALES (SIN ATRIBUTOS EXTERNOS)
# ============================================================================

# Limpiar entorno
rm(list=ls())

# Instalar paquetes si no están instalados
paquetes <- c("igraph", "glasso", "qgraph", "corrplot")
for(pkg in paquetes) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Cargar librerías explícitamente
library(igraph)
library(glasso)
library(qgraph)
library(corrplot)

# ----------------------------------------------------------------------------
# 1. CARGAR Y PREPARAR LOS DATOS
# ----------------------------------------------------------------------------

# Ver directorio actual
# cat("Directorio de trabajo actual:", getwd(), "\n\n")

# OPCIÓN A: Si el archivo está en el directorio de trabajo
edges <- read.csv("redinventada.txt", sep="\t", header=TRUE)

# OPCIÓN B: Seleccionar archivo manualmente (más fácil)
#cat("Seleccioná el archivo redinventada.txt en la ventana que se abre...\n")
#archivo <- file.choose()
#edges <- read.csv(archivo, sep="\t", header=TRUE)

# OPCIÓN C: Si preferís especificar la ruta completa, descomentá y modificá:
# edges <- read.csv("C:/Users/TuUsuario/Documents/redinventada.txt", sep="\t", header=TRUE)

# Crear el grafo
g <- graph_from_data_frame(edges, directed=FALSE)

# Ver información básica
cat("========== INFORMACIÓN DE LA RED ==========\n")
cat("Número de nodos:", vcount(g), "\n")
cat("Número de aristas:", ecount(g), "\n")
cat("\nTipos de relaciones:\n")
print(table(edges$label))

# ----------------------------------------------------------------------------
# 2. ENFOQUE 1: GLASSO SOBRE PATRONES DE VECINDAD
# ----------------------------------------------------------------------------
# Pregunta: ¿Qué nodos tienen patrones estructurales directamente similares?

cat("\n\n========== ENFOQUE 1: VECINDARIOS DE NODOS ==========\n")

# Seleccionar nodos con suficiente conectividad
min_degree <- 3  # Nodos con al menos 3 conexiones
connected_nodes <- V(g)[degree(g) >= min_degree]

# Si hay muchos nodos, tomar los más conectados
if(length(connected_nodes) > 50) {
  top_n <- 50
  connected_nodes <- V(g)[order(degree(g), decreasing=TRUE)[1:top_n]]
  cat("Analizando los", top_n, "nodos más conectados\n")
} else {
  cat("Analizando", length(connected_nodes), "nodos con degree >=", min_degree, "\n")
}

subg <- induced_subgraph(g, connected_nodes)

# Crear matriz de adyacencia
adj_matrix <- as_adjacency_matrix(subg, sparse=FALSE)

# Correlación entre patrones de conexión de nodos
# Cada nodo es una "variable", sus conexiones son "observaciones"
S_nodes <- cor(adj_matrix)

# Aplicar GLASSO con diferentes valores de lambda
cat("\nProbando diferentes valores de regularización...\n")
lambda_values <- c(0.05, 0.1, 0.2, 0.3)

par(mfrow=c(2,2))
for(lambda in lambda_values) {
  glasso_result <- glasso(S_nodes, rho=lambda)
  n_edges <- sum(glasso_result$wi[upper.tri(glasso_result$wi)] != 0)
  
  qgraph(glasso_result$wi,
         layout="spring",
         title=paste0("Lambda = ", lambda, " (", n_edges, " conexiones)"),
         labels=V(subg)$name,
         vsize=4,
         label.cex=0.8,
         cut=0.01)
  
  cat(sprintf("  Lambda=%.2f -> %d conexiones directas\n", lambda, n_edges))
}

# Elegir un lambda intermedio para análisis detallado
best_lambda <- 0.1
glasso_nodes <- glasso(S_nodes, rho=best_lambda)

cat("\n--- Interpretación (Lambda =", best_lambda, ") ---\n")
cat("Este grafo muestra nodos con ROLES ESTRUCTURALES similares.\n")
cat("Dos nodos están conectados si tienen patrones de conexión directamente relacionados,\n")
cat("no solo porque estén en la misma comunidad o cerca en el grafo.\n")

# Identificar nodos centrales en el grafo GLASSO
prec_matrix <- glasso_nodes$wi
glasso_degree <- colSums(abs(prec_matrix) > 0.01) - 1  # -1 para diagonal
top_glasso <- order(glasso_degree, decreasing=TRUE)[1:min(5, length(glasso_degree))]

cat("\nNodos con más dependencias estructurales directas:\n")
for(i in top_glasso) {
  cat(sprintf("  Nodo %s: %d conexiones directas (degree original: %d)\n",
              V(subg)$name[i], 
              glasso_degree[i],
              degree(subg)[i]))
}

# ----------------------------------------------------------------------------
# 3. ENFOQUE 2: GLASSO SOBRE TIPOS DE RELACIÓN
# ----------------------------------------------------------------------------
# Pregunta: ¿Qué tipos de relaciones co-ocurren directamente?

cat("\n\n========== ENFOQUE 2: TIPOS DE RELACIÓN ==========\n")

# Crear matriz: nodo x tipo_de_relación
relaciones <- unique(edges$label[edges$label != ""])
cat("Tipos de relación encontrados:", length(relaciones), "\n")

relation_matrix <- matrix(0, 
                          nrow=vcount(g), 
                          ncol=length(relaciones),
                          dimnames=list(V(g)$name, relaciones))

# Llenar la matriz
for(i in 1:nrow(edges)) {
  if(edges$label[i] != "") {
    src <- as.character(edges$source[i])
    tgt <- as.character(edges$target[i])
    rel <- edges$label[i]
    
    if(src %in% V(g)$name) relation_matrix[src, rel] <- 1
    if(tgt %in% V(g)$name) relation_matrix[tgt, rel] <- 1
  }
}

# Ver distribución
cat("\nFrecuencia de cada tipo de relación:\n")
rel_counts <- sort(colSums(relation_matrix), decreasing=TRUE)
print(rel_counts)

# Filtrar relaciones frecuentes (aparecen en al menos 5 nodos)
frequent_rels <- names(rel_counts[rel_counts >= 5])

if(length(frequent_rels) >= 3) {
  cat("\nAnalizando", length(frequent_rels), "tipos de relación frecuentes\n")
  
  rel_subset <- relation_matrix[, frequent_rels]
  
  # Correlación entre tipos de relación
  S_rel <- cor(rel_subset)
  
  # Visualizar correlaciones
  par(mfrow=c(1,2))
  corrplot(S_rel, method="color", type="upper",
           title="Correlaciones entre tipos de relación",
           mar=c(0,0,2,0),
           tl.cex=0.7)
  
  # Aplicar GLASSO
  glasso_rel <- glasso(S_rel, rho=0.1)
  
  # Normalizar matriz de precisión para visualización
  prec_rel_norm <- cov2cor(glasso_rel$wi)
  corrplot(prec_rel_norm, method="color", type="upper",
           title="Dependencias parciales (GLASSO)",
           mar=c(0,0,2,0),
           tl.cex=0.7)
  
  # Visualizar como red
  par(mfrow=c(1,1))
  qgraph(glasso_rel$wi,
         layout="spring",
         labels=colnames(rel_subset),
         title="Red de Dependencias entre Tipos de Relación",
         vsize=10,
         label.cex=1.2,
         esize=5,
         cut=0.05)
  
  cat("\n--- Interpretación ---\n")
  cat("Este grafo muestra qué tipos de relaciones co-ocurren DIRECTAMENTE.\n")
  cat("Por ejemplo, si 'familiar' y 'escuchas' están conectados, significa que\n")
  cat("personas con vínculos familiares tienden a estar en escuchas,\n")
  cat("más allá de cualquier otro factor intermedio.\n")
  
  # Analizar conexiones específicas
  prec_rel <- glasso_rel$wi
  strong_connections <- which(abs(prec_rel[upper.tri(prec_rel)]) > 0.1)
  
  if(length(strong_connections) > 0) {
    cat("\nConexiones directas fuertes entre tipos de relación:\n")
    indices <- which(upper.tri(prec_rel), arr.ind=TRUE)
    for(i in strong_connections) {
      rel1 <- colnames(rel_subset)[indices[i,1]]
      rel2 <- colnames(rel_subset)[indices[i,2]]
      strength <- prec_rel[indices[i,1], indices[i,2]]
      cat(sprintf("  %s <-> %s: %.3f\n", rel1, rel2, strength))
    }
  }
  
} else {
  cat("\nNo hay suficientes tipos de relación frecuentes para análisis GLASSO\n")
  cat("(Se necesitan al menos 3 tipos con más de 5 ocurrencias)\n")
}

# ----------------------------------------------------------------------------
# 4. RESUMEN Y CONCLUSIONES
# ----------------------------------------------------------------------------

cat("\n\n========== RESUMEN ==========\n")
cat("GLASSO te permite descubrir:\n\n")
cat("1. ENFOQUE DE VECINDARIOS:\n")
cat("   - Nodos con roles estructurales similares\n")
cat("   - Puentes ocultos entre comunidades\n")
cat("   - Equivalencias estructurales más allá de proximidad\n\n")
cat("2. ENFOQUE DE TIPOS DE RELACIÓN:\n")
cat("   - Qué tipos de vínculos van juntos directamente\n")
cat("   - Relaciones fundamentales vs. derivadas\n")
cat("   - Patrones de co-ocurrencia no obvios\n\n")
cat("Las conexiones en GLASSO son DIRECTAS, no transitivas ni espurias.\n")