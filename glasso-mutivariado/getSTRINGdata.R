# ============================================================================
# EGO-NETWORK DE 2 PASOS DESDE STRING
# Input: Genes semilla (BRCA1, BRCA2)
# Output: Lista expandida de genes (vecinos + vecinos de vecinos)
# ============================================================================

# Limpiar entorno
rm(list=ls())

# ----------------------------------------------------------------------------
# 0. INSTALACIÓN Y CARGA DE PAQUETES
# ----------------------------------------------------------------------------

cat("========== INSTALANDO PAQUETES ==========\n")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("STRINGdb", quietly = TRUE)) {
  cat("Instalando STRINGdb...\n")
  BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
}

library(STRINGdb)

cat("✓ Paquetes cargados\n\n")

# ----------------------------------------------------------------------------
# 1. CONFIGURACIÓN
# ----------------------------------------------------------------------------

# Genes semilla (punto de partida)
seed_genes <- c("BRCA1", "BRCA2")

# Threshold de score en STRING (400 = 0.4, 700 = 0.7, 900 = 0.9)
score_threshold <- 700  # Ajustá según qué tan estricto querés ser

# Número máximo de vecinos por nodo en cada paso
max_neighbors_step1 <- 30  # Vecinos directos de BRCA1/2
max_neighbors_step2 <- 20  # Vecinos de cada vecino

cat("========== CONFIGURACIÓN ==========\n")
cat("Genes semilla:", paste(seed_genes, collapse=", "), "\n")
cat("Score threshold:", score_threshold/1000, "\n")
cat("Max vecinos paso 1:", max_neighbors_step1, "\n")
cat("Max vecinos paso 2:", max_neighbors_step2, "\n\n")

# ----------------------------------------------------------------------------
# 2. CONECTAR A STRING
# ----------------------------------------------------------------------------

cat("========== CONECTANDO A STRING ==========\n")
cat("(Esto puede tardar un momento la primera vez...)\n")

string_db <- STRINGdb$new(version="12.0", 
                          species=9606,  # Homo sapiens
                          score_threshold=score_threshold,
                          input_directory="")

cat("✓ Conexión establecida\n\n")

# ----------------------------------------------------------------------------
# 3. PASO 1: VECINOS DIRECTOS DE GENES SEMILLA
# ----------------------------------------------------------------------------

cat("========== PASO 1: VECINOS DIRECTOS ==========\n")

# Mapear genes semilla a IDs de STRING
seed_mapped <- string_db$map(data.frame(gene=seed_genes), "gene")

if(nrow(seed_mapped) == 0) {
  stop("ERROR: No se pudieron mapear los genes semilla a STRING")
}

cat("Genes semilla mapeados:\n")
print(seed_mapped)

# Obtener vecinos de cada gen semilla
step1_neighbors <- c()

for(i in 1:nrow(seed_mapped)) {
  gene_name <- seed_mapped$gene[i]
  gene_id <- seed_mapped$STRING_id[i]
  
  cat(sprintf("\nBuscando vecinos de %s...\n", gene_name))
  
  # Obtener vecinos directos
  neighbors <- string_db$get_neighbors(gene_id)
  
  cat(sprintf("  Encontrados: %d vecinos\n", length(neighbors)))
  
  # Limitar cantidad si es necesario
  if(length(neighbors) > max_neighbors_step1) {
    # Obtener scores y tomar los top
    interactions <- string_db$get_interactions(c(gene_id, neighbors))
    interactions <- interactions[interactions$from == gene_id | interactions$to == gene_id, ]
    interactions <- interactions[order(-interactions$combined_score), ]
    neighbors <- unique(c(interactions$from[1:max_neighbors_step1], 
                          interactions$to[1:max_neighbors_step1]))
    neighbors <- neighbors[neighbors != gene_id]
    neighbors <- head(neighbors, max_neighbors_step1)
    
    cat(sprintf("  Limitado a top %d\n", length(neighbors)))
  }
  
  step1_neighbors <- c(step1_neighbors, neighbors)
}

# Eliminar duplicados y genes semilla
step1_neighbors <- unique(step1_neighbors)
step1_neighbors <- step1_neighbors[!step1_neighbors %in% seed_mapped$STRING_id]

cat(sprintf("\n✓ Total vecinos únicos en Paso 1: %d\n", length(step1_neighbors)))

# ----------------------------------------------------------------------------
# 4. PASO 2: VECINOS DE VECINOS
# ----------------------------------------------------------------------------

cat("\n========== PASO 2: VECINOS DE VECINOS ==========\n")
cat("(Esto puede tardar varios minutos...)\n")

step2_neighbors <- c()

# Tomar subset de vecinos más importantes del paso 1 para expandir
# (expandir desde todos sería demasiado)
n_to_expand <- min(10, length(step1_neighbors))

cat(sprintf("Expandiendo desde los top %d vecinos del Paso 1...\n", n_to_expand))

# Obtener scores de paso 1 para priorizar
all_seed_ids <- seed_mapped$STRING_id
interactions_step1 <- string_db$get_interactions(c(all_seed_ids, step1_neighbors))

# Calcular score promedio de cada vecino con genes semilla
neighbor_scores <- sapply(step1_neighbors, function(neighbor_id) {
  relevant <- interactions_step1[
    (interactions_step1$from == neighbor_id & interactions_step1$to %in% all_seed_ids) |
      (interactions_step1$to == neighbor_id & interactions_step1$from %in% all_seed_ids), ]
  if(nrow(relevant) > 0) {
    return(mean(relevant$combined_score))
  } else {
    return(0)
  }
})

# Tomar los n_to_expand con mayor score
top_to_expand <- names(sort(neighbor_scores, decreasing=TRUE))[1:n_to_expand]

for(i in 1:length(top_to_expand)) {
  neighbor_id <- top_to_expand[i]
  
  if(i %% 3 == 0) {  # Progress update cada 3
    cat(sprintf("  Procesando %d/%d...\n", i, length(top_to_expand)))
  }
  
  # Obtener vecinos de este vecino
  neighbors_of_neighbor <- string_db$get_neighbors(neighbor_id)
  
  # Limitar cantidad
  if(length(neighbors_of_neighbor) > max_neighbors_step2) {
    neighbors_of_neighbor <- head(neighbors_of_neighbor, max_neighbors_step2)
  }
  
  step2_neighbors <- c(step2_neighbors, neighbors_of_neighbor)
}

# Eliminar duplicados
step2_neighbors <- unique(step2_neighbors)

# Eliminar genes que ya están en paso 1 o semilla
step2_neighbors <- step2_neighbors[!step2_neighbors %in% c(step1_neighbors, seed_mapped$STRING_id)]

cat(sprintf("\n✓ Nuevos vecinos en Paso 2: %d\n", length(step2_neighbors)))

# ----------------------------------------------------------------------------
# 5. COMBINAR Y MAPEAR A SÍMBOLOS DE GENES
# ----------------------------------------------------------------------------

cat("\n========== GENERANDO LISTA FINAL ==========\n")

# Combinar todos los IDs
all_ids <- unique(c(seed_mapped$STRING_id, step1_neighbors, step2_neighbors))

cat(sprintf("Total de proteínas en la red: %d\n", length(all_ids)))

# Mapear de vuelta a símbolos de genes
all_proteins <- string_db$get_proteins()
final_genes <- all_proteins[all_proteins$protein_external_id %in% all_ids, ]

# Usar preferred_name (símbolo del gen)
gene_list <- unique(final_genes$preferred_name)

# Ordenar alfabéticamente
gene_list <- sort(gene_list)

cat(sprintf("\n✓ Lista final: %d genes\n\n", length(gene_list)))

# ----------------------------------------------------------------------------
# 6. RESULTADOS
# ----------------------------------------------------------------------------

cat("========== LISTA DE GENES ==========\n")
cat("(Genes ordenados alfabéticamente)\n\n")

# Imprimir en columnas de 5 para legibilidad
genes_per_row <- 5
for(i in seq(1, length(gene_list), by=genes_per_row)) {
  end_idx <- min(i + genes_per_row - 1, length(gene_list))
  cat(paste(sprintf("%-12s", gene_list[i:end_idx]), collapse=""), "\n")
}

# ----------------------------------------------------------------------------
# 7. ESTADÍSTICAS
# ----------------------------------------------------------------------------

cat("\n========== ESTADÍSTICAS ==========\n")
cat(sprintf("Genes semilla: %d\n", length(seed_genes)))
cat(sprintf("Vecinos directos (Paso 1): %d\n", length(step1_neighbors)))
cat(sprintf("Vecinos de vecinos (Paso 2): %d\n", length(step2_neighbors)))
cat(sprintf("TOTAL: %d genes\n", length(gene_list)))

# Calcular qué % de crecimiento hubo
growth_pct <- ((length(gene_list) - length(seed_genes)) / length(seed_genes)) * 100
cat(sprintf("Expansión: %.0f%%\n\n", growth_pct))

# ----------------------------------------------------------------------------
# 8. EXPORTAR A ARCHIVO
# ----------------------------------------------------------------------------

cat("========== EXPORTANDO ==========\n")

# Crear data frame con información adicional
output_df <- data.frame(
  gene_symbol = gene_list,
  step = ifelse(gene_list %in% seed_genes, "seed",
                ifelse(final_genes$protein_external_id[match(gene_list, final_genes$preferred_name)] %in% step1_neighbors,
                       "step1", "step2"))
)

# Guardar como CSV
output_file <- "brca_ego_network_genes.csv"
write.csv(output_df, output_file, row.names=FALSE)

cat(sprintf("✓ Lista guardada en: %s\n", output_file))

# También guardar solo los símbolos (para copiar/pegar fácilmente)
output_file_simple <- "brca_ego_network_genes_simple.txt"
writeLines(gene_list, output_file_simple)

cat(sprintf("✓ Lista simple guardada en: %s\n\n", output_file_simple))

# ----------------------------------------------------------------------------
# 9. CREAR VECTOR PARA R
# ----------------------------------------------------------------------------

cat("========== VECTOR PARA USAR EN R ==========\n")
cat("Copiá y pegá esto en tu script principal:\n\n")

# Generar código R
cat("genes_of_interest <- c(\n")
for(i in 1:length(gene_list)) {
  if(i < length(gene_list)) {
    cat(sprintf('  "%s",\n', gene_list[i]))
  } else {
    cat(sprintf('  "%s"\n', gene_list[i]))
  }
  
  # No imprimir toda la lista si es muy larga
  if(i == 50 && length(gene_list) > 50) {
    cat(sprintf('  # ... y %d genes más (ver archivo CSV)\n', length(gene_list) - 10))
    break
  }
}
cat(")\n\n")

# ----------------------------------------------------------------------------
# 10. ANÁLISIS ADICIONAL (OPCIONAL)
# ----------------------------------------------------------------------------

cat("========== ANÁLISIS ADICIONAL ==========\n")

# Ver qué genes del paso 1 son más conectados
if(length(step1_neighbors) > 0) {
  cat("\nTop 10 vecinos directos más conectados:\n")
  
  step1_genes <- final_genes[final_genes$protein_external_id %in% step1_neighbors, ]
  
  # Contar cuántos vecinos tiene cada uno
  step1_connectivity <- sapply(step1_neighbors, function(id) {
    length(string_db$get_neighbors(id))
  })
  
  top_connected <- sort(step1_connectivity, decreasing=TRUE)[1:min(10, length(step1_connectivity))]
  
  for(i in 1:length(top_connected)) {
    gene_id <- names(top_connected)[i]
    gene_name <- step1_genes$preferred_name[step1_genes$protein_external_id == gene_id]
    cat(sprintf("%2d. %s (%d conexiones)\n", i, gene_name, top_connected[i]))
  }
}

cat("\n========== PROCESO COMPLETADO ==========\n")
cat("Archivos generados:\n")
cat(sprintf("  1. %s (con información de pasos)\n", output_file))
cat(sprintf("  2. %s (solo nombres de genes)\n", output_file_simple))
cat("\nUsá estos genes en tu análisis GLASSO principal.\n")