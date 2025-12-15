# ============================================================================
# EGO-NETWORK DE 2 PASOS DESDE STRING (CON SCORES Y ORDENAMIENTO)
# Input: Genes semilla (BRCA1, BRCA2)
# Output: Dataframe con score de interacción y jerarquía
# ============================================================================

rm(list=ls())

# ----------------------------------------------------------------------------
# 1. SETUP
# ----------------------------------------------------------------------------
if(!require("STRINGdb", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("STRINGdb", update = FALSE, ask = FALSE)
}
library(STRINGdb)

# PARÁMETROS
seed_genes <- c("BRCA1", "BRCA2")
score_threshold <- 700       # Confianza alta
max_neighbors_step1 <- 30    # Top N vecinos directos
max_neighbors_step2 <- 20    # Top N vecinos de segunda linea

# Inicializar DB
cat("Conectando a STRING 12.0...\n")
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=score_threshold, input_directory="")

# ----------------------------------------------------------------------------
# 2. PROCESAMIENTO
# ----------------------------------------------------------------------------

# --- A. MAPEAR SEMILLAS ---
seed_mapped <- string_db$map(data.frame(gene=seed_genes), "gene")
# Inicializamos el dataframe maestro con las semillas
master_df <- data.frame(
  STRING_id = seed_mapped$STRING_id,
  gene = seed_mapped$gene,
  step = 0,
  origin_gene = "SEED",
  score = 1000 # Score máximo para las semillas
)

cat(sprintf("\nGenes semilla mapeados: %s\n", paste(seed_mapped$gene, collapse=", ")))

# --- B. PASO 1: VECINOS DIRECTOS ---
cat("\n>>> PASO 1: Buscando vecinos directos de BRCA1/2...\n")

# Obtenemos TODOS los vecinos que cumplan el threshold
# Usamos get_interactions para obtener el score inmediatamente
candidates_step1 <- list()

for(i in 1:nrow(seed_mapped)) {
  seed_id <- seed_mapped$STRING_id[i]
  seed_name <- seed_mapped$gene[i]
  
  # Traer vecinos
  neighbors <- string_db$get_neighbors(seed_id)
  
  if(length(neighbors) > 0) {
    # Traer interacciones para ver scores
    interactions <- string_db$get_interactions(c(seed_id, neighbors))
    
    # Filtrar solo las conexiones directas con la semilla actual
    links <- interactions[(interactions$from == seed_id & interactions$to %in% neighbors) |
                            (interactions$to == seed_id & interactions$from %in% neighbors), ]
    
    # Ordenar por score
    links <- links[order(-links$combined_score), ]
    
    # Quedarnos con el top N
    links <- head(links, max_neighbors_step1)
    
    # Identificar al vecino (el que NO es la semilla)
    neighbor_ids <- ifelse(links$from == seed_id, links$to, links$from)
    
    # Guardar info temporal
    df_temp <- data.frame(
      STRING_id = neighbor_ids,
      step = 1,
      origin_gene = seed_name, # Guardamos quién lo trajo
      score = links$combined_score
    )
    candidates_step1[[seed_name]] <- df_temp
  }
}

# Unir candidatos del paso 1
df_step1 <- do.call(rbind, candidates_step1)

# Limpieza: Si un gen es vecino de BRCA1 y BRCA2, nos quedamos con la conexión más fuerte
df_step1 <- df_step1[order(-df_step1$score), ] # Ordenar por score
df_step1 <- df_step1[!duplicated(df_step1$STRING_id), ] # Eliminar duplicados (queda el mejor)

# Eliminar si ya estaba en semillas (auto-bucles)
df_step1 <- df_step1[!df_step1$STRING_id %in% master_df$STRING_id, ]

cat(sprintf("✓ Paso 1 completado: %d vecinos únicos agregados.\n", nrow(df_step1)))

# --- C. PASO 2: VECINOS DE VECINOS ---
cat("\n>>> PASO 2: Expandiendo desde los mejores vecinos...\n")

# Seleccionamos los top 10 genes del paso 1 para expandir (para no explotar la red)
top_expanders <- head(df_step1$STRING_id, 10)
candidates_step2 <- list()

counter <- 0
for(expander_id in top_expanders) {
  counter <- counter + 1
  # Recuperar nombre del gen origen para referencia
  expander_name <- df_step1$origin_gene[df_step1$STRING_id == expander_id] # Es aproximado, aun no tenemos mapped name del expander
  
  neighbors <- string_db$get_neighbors(expander_id)
  
  if(length(neighbors) > 0) {
    interactions <- string_db$get_interactions(c(expander_id, neighbors))
    links <- interactions[(interactions$from == expander_id & interactions$to %in% neighbors) |
                            (interactions$to == expander_id & interactions$from %in% neighbors), ]
    
    links <- links[order(-links$combined_score), ]
    links <- head(links, max_neighbors_step2)
    
    neighbor_ids <- ifelse(links$from == expander_id, links$to, links$from)
    
    df_temp <- data.frame(
      STRING_id = neighbor_ids,
      step = 2,
      origin_id_temp = expander_id, # Usamos ID temporalmente
      score = links$combined_score
    )
    candidates_step2[[counter]] <- df_temp
  }
}

df_step2 <- do.call(rbind, candidates_step2)

# Limpieza Paso 2
if(!is.null(df_step2)) {
  df_step2 <- df_step2[order(-df_step2$score), ]
  df_step2 <- df_step2[!duplicated(df_step2$STRING_id), ]
  
  # IMPORTANTE: No incluir genes que ya están en Semillas o Paso 1
  ids_existentes <- c(master_df$STRING_id, df_step1$STRING_id)
  df_step2 <- df_step2[!df_step2$STRING_id %in% ids_existentes, ]
  
  # Arreglar columna de origen (mapear ID -> Nombre genérico "NeighborLvl1")
  # Para hacerlo simple en el output, ponemos "Neighbor_Lvl1"
  df_step2$origin_gene <- "From_Step1"
  df_step2$origin_id_temp <- NULL
  
  cat(sprintf("✓ Paso 2 completado: %d vecinos lejanos agregados.\n", nrow(df_step2)))
}

# ----------------------------------------------------------------------------
# 3. CONSOLIDACIÓN Y MAPEO
# ----------------------------------------------------------------------------

cat("\n>>> Generando tabla final...\n")

# Unir todo sin las columnas de nombre aun (porque df_step1 y 2 no lo tienen)
cols_common <- c("STRING_id", "step", "origin_gene", "score")

# Preparar master para merge
part1 <- master_df[, cols_common]
part2 <- df_step1[, cols_common]
part3 <- if(!is.null(df_step2)) df_step2[, cols_common] else NULL

final_network <- rbind(part1, part2, part3)

# Mapear IDs a Nombres Reales
all_proteins <- string_db$get_proteins()
indices <- match(final_network$STRING_id, all_proteins$protein_external_id)
final_network$gene_symbol <- all_proteins$preferred_name[indices]

# Reordenar columnas
final_network <- final_network[, c("gene_symbol", "step", "score", "origin_gene", "STRING_id")]

# ORDENAMIENTO FINAL PARA EL TXT (Tu requerimiento)
# 1. Por Step (0, 1, 2)
# 2. Por Score (Descendiente)
final_network <- final_network[order(final_network$step, -final_network$score), ]

# ----------------------------------------------------------------------------
# 4. EXPORTACIÓN
# ----------------------------------------------------------------------------

# Archivo 1: CSV Completo con Scores
file_csv <- "brca_ego_network_data.csv"
write.csv(final_network, file_csv, row.names = FALSE)
cat(sprintf("\n✓ CSV guardado: %s (incluye scores)\n", file_csv))

# Archivo 2: TXT Simple Ordenado
file_txt <- "brca_ego_network_genes_simple.txt"
# Solo guardamos la columna de genes
writeLines(final_network$gene_symbol, file_txt)
cat(sprintf("✓ TXT guardado: %s (Ordenado: Seed -> Step1(Score) -> Step2(Score))\n", file_txt))

# ----------------------------------------------------------------------------
# 5. VISUALIZACIÓN DE MUESTRA
# ----------------------------------------------------------------------------
cat("\nTop 15 genes de tu red (ordenados por importancia):\n")
print(head(final_network[, c("gene_symbol", "step", "score", "origin_gene")], 15))

cat("\nTotal genes para analizar:", nrow(final_network), "\n")
