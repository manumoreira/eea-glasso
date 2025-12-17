#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


library(shiny)
library(tidyverse)
library(data.table)
library(glasso)
library(corrplot)
library(plotly)
library(stringr)

#========================================================
# CARGA DE DATOS + LIMPIEZA
#========================================================
file_path <- "player_stats.csv"

df <- read.csv(file_path)

df <- df %>%
  mutate(
    value_clean = str_remove_all(value, "\\$"),
    value_clean = str_replace(value_clean, "\\.00", "000"),
    value_clean = str_remove_all(value_clean, "\\."),
    value_num   = as.numeric(value_clean)
  ) %>%
  select(-value_clean)

df_num <- df[, sapply(df, is.numeric)]
vars <- colnames(df_num)
m_cor <- cor(df_num, use = "complete.obs")

#========================================================
# FUNCIONES AUXILIARES
#========================================================

build_edge_list <- function(PN, vars) {
  edges <- data.table()
  
  for (i in 1:(nrow(PN)-1)) {
    for (j in (i+1):ncol(PN)) {
      if (abs(PN[i,j]) > 1e-6) {
        edges <- rbind(edges, data.table(
          from = vars[i], 
          to   = vars[j],
          weight = PN[i,j]
        ))
      }
    }
  }
  
  edges
}

create_3d_graph <- function(PN, vars, selected_node = NULL,
                            min_node_distance = 1.5,
                            min_line_width = 1,
                            max_line_width = 5) {
  # Crear lista de edges
  edges <- build_edge_list(PN, vars)
  
  if (nrow(edges) == 0) {
    return(NULL)
  }
  
  # Distancias (1 - |rho_parcial|)
  distancia <- as.dist(1 - abs(PN))
  mds_3d <- cmdscale(distancia, k = 3)
  
  posiciones <- data.table(
    var = vars,
    x = mds_3d[,1],
    y = mds_3d[,2],
    z = mds_3d[,3]
  )
  
  # ALGORITMO DE SEPARACIÓN MÍNIMA
  max_iter <- 100
  for (iter in 1:max_iter) {
    moved <- FALSE
    
    for (i in 1:nrow(posiciones)) {
      for (j in 1:nrow(posiciones)) {
        if (i != j) {
          dx <- posiciones$x[i] - posiciones$x[j]
          dy <- posiciones$y[i] - posiciones$y[j]
          dz <- posiciones$z[i] - posiciones$z[j]
          dist <- sqrt(dx^2 + dy^2 + dz^2)
          
          if (dist < min_node_distance && dist > 0) {
            factor <- (min_node_distance - dist) / dist * 0.1
            
            posiciones$x[i] <- posiciones$x[i] + dx * factor
            posiciones$y[i] <- posiciones$y[i] + dy * factor
            posiciones$z[i] <- posiciones$z[i] + dz * factor
            
            moved <- TRUE
          }
        }
      }
    }
    
    if (!moved) break
  }
  
  # Escalar
  posiciones[, x := x * 2]
  posiciones[, y := y * 2]
  posiciones[, z := z * 2]
  
  # Agregar coords a edges
  edges <- merge(edges, posiciones, by.x = "from", by.y = "var")
  setnames(edges, c("x", "y", "z"), c("x1", "y1", "z1"))
  
  edges <- merge(edges, posiciones, by.x = "to", by.y = "var")
  setnames(edges, c("x", "y", "z"), c("x2", "y2", "z2"))
  
  # Normalizar pesos para ancho de línea
  edges[, abs_weight := abs(weight)]
  max_weight <- max(edges$abs_weight)
  edges[, line_width := min_line_width + (abs_weight / max_weight) * (max_line_width - min_line_width)]
  
  # Vecinos del nodo seleccionado (si existe)
  neighbors <- character(0)
  if (!is.null(selected_node)) {
    neighbors <- unique(c(
      edges[from == selected_node, to],
      edges[to   == selected_node, from]
    ))
  }
  
  # Crear plot con fuente identificable
  fig <- plot_ly(source = "graph3D")
  
  # Agregar edges (líneas)
  for (i in 1:nrow(edges)) {
    is_incident <- is.null(selected_node) ||
      edges$from[i] == selected_node ||
      edges$to[i]   == selected_node
    
    edge_opacity <- if (is.null(selected_node)) 0.6 else if (is_incident) 0.9 else 0.05
    
    fig <- fig %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = c(edges$x1[i], edges$x2[i]),
        y = c(edges$y1[i], edges$y2[i]),
        z = c(edges$z1[i], edges$z2[i]),
        line = list(
          color = ifelse(edges$weight[i] > 0, "rgba(255,0,0,1)", "rgba(0,0,255,1)"),
          width = edges$line_width[i]
        ),
        opacity = edge_opacity,
        hoverinfo = "text",
        text = paste0(
          edges$from[i], " - ", edges$to[i], 
          ": ", round(edges$weight[i], 3)
        ),
        showlegend = FALSE
      )
  }
  
  # Preparar estilos de nodos
  posiciones[, is_selected := (!is.null(selected_node) & var == selected_node)]
  posiciones[, is_neighbor := var %in% neighbors]
  
  if (is.null(selected_node)) {
    posiciones[, marker_size := 12]
    posiciones[, marker_opacity := 1]
  } else {
    posiciones[, marker_size := ifelse(
      is_selected, 14,
      ifelse(is_neighbor, 10, 6)
    )]
    
    posiciones[, marker_opacity := ifelse(
      is_selected | is_neighbor, 0.9,
      0.15
    )]
  }
  
  # Nodos
  fig <- fig %>%
    add_trace(
      data = posiciones,
      type = "scatter3d",
      mode = "markers+text",
      x = ~x, y = ~y, z = ~z,
      marker = list(
        size = posiciones$marker_size,
        color = "steelblue",
        opacity = posiciones$marker_opacity,
        line = list(color = "white", width = 2)
      ),
      text = ~var,
      textposition = "top center",
      textfont = list(size = 8, color = "black"),
      hoverinfo = "text",
      showlegend = FALSE
    )

  axis_off <- list(
    title = "",
    showgrid = FALSE,
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    ticks = "",
    showspikes = FALSE
  )
  
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = axis_off,
        yaxis = axis_off,
        zaxis = axis_off,
        camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
      ),
      showlegend = FALSE,
      hovermode = "closest"
    )
  
  
  return(fig)
}

#========================================================
# UI
#========================================================
ui <- fluidPage(
  
  titlePanel("Glasso Explorer — Matrices, Red 3D y Conexiones Perdidas"),
  
  fluidRow(
    
    #-----------------------------------
    # COLUMNA IZQUIERDA (CONTROLES)
    #-----------------------------------
    column(
      width = 3,
      wellPanel(
        sliderInput("lambda",
                    "Lambda (rho) para Glasso:",
                    min = 0.001,
                    max = 1.0,
                    value = 0.3,
                    step = 0.01),
        br(),
        actionButton("reset_node", "Mostrar red completa"),
        br(), br(),
        h5("Información:"),
        textOutput("info_edges"),
        textOutput("info_selected")
      )
    ),
    
    #-----------------------------------
    # COLUMNA DERECHA (CONTENIDO)
    #-----------------------------------
    column(
      width = 9,
      
      tabsetPanel(
        
        #====================================================
        # DATASET
        #====================================================
        tabPanel("Dataset",
                 br(),
                 h4("Información del Dataset"),
                 textOutput("dataset_info"),
                 br(),
                 h5("Primeras 15 filas:"),
                 div(
                   style = "overflow-x: auto;",
                   tableOutput("dataset_preview")
                 )
        ),
        
        #====================================================
        # MATRICES
        #====================================================
        tabPanel("Matrices",
                 
                 fluidRow(
                   column(
                     width = 6,
                     h4("Matriz de Correlaciones"),
                     plotOutput("heat_cor", height = "650px", width = "650px")
                   ),
                   
                   column(
                     width = 6,
                     h4("Matriz de Precisión Glasso (Partial Corr)"),
                     plotOutput("heat_prec", height = "650px", width = "650px")
                   )
                 )
        ),
        
        #====================================================
        # RED 3D
        #====================================================
        tabPanel("Red Glasso 3D",
                 plotlyOutput("graph3D", height = "700px")
        ),
        
        #====================================================
        # TABLA DE CONEXIONES PERDIDAS
        #====================================================
        tabPanel("Conexiones Perdidas",
                 h4("Top 20 correlaciones perdidas y sus mediadores"),
                 p("Las variables en 'mediadores' son las que explican la correlación indirecta entre var1 y var2."),
                 tableOutput("tabla_perdidas")
        )
      )
    )
  )
)

#========================================================
# SERVER
#========================================================
server <- function(input, output, session) {
  
  # Nodo seleccionado en el grafo (NULL = red completa)
  selected_node <- reactiveVal(NULL)
  
  observeEvent(input$reset_node, {
    selected_node(NULL)
  })
  
  # GLASSO
  gl_fit <- reactive({
    glasso(m_cor, rho = input$lambda)
  })
  
  # PRECISION NORMALIZADA (PARTIAL CORR)
  precision_norm <- reactive({
    W <- gl_fit()$wi
    PN <- W / sqrt(outer(diag(W), diag(W)))
    rownames(PN) <- vars
    colnames(PN) <- vars
    PN
  })
  
  #--------------------------------------------------------
  # INFO DEL DATASET
  #--------------------------------------------------------
  output$dataset_info <- renderText({
    paste0("El dataset contiene ", nrow(df), " filas y ", ncol(df), " columnas",
           " (", length(vars), " variables numéricas)")
  })
  
  output$dataset_preview <- renderTable({
    head(df, 15)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  #--------------------------------------------------------
  # INFO DE EDGES
  #--------------------------------------------------------
  output$info_edges <- renderText({
    PN <- precision_norm()
    edges <- build_edge_list(PN, vars)
    paste("Conexiones activas totales:", nrow(edges))
  })
  
  output$info_selected <- renderText({
    sn <- selected_node()
    if (is.null(sn)) {
      "Nodo seleccionado: ninguno (se muestra toda la red)"
    } else {
      PN <- precision_norm()
      edges <- build_edge_list(PN, vars)
      deg <- nrow(edges[from == sn | to == sn])
      paste0("Nodo seleccionado: ", sn, " (", deg, " conexiones)")
    }
  })
  
  #--------------------------------------------------------
  # HEATMAP DE CORRELACIONES
  #--------------------------------------------------------
  output$heat_cor <- renderPlot({
    corrplot(m_cor,
             method="circle",
             type="upper",
             diag=FALSE,
             tl.cex=1,
             tl.col="black")
  }, height = 650, width = 650)
  
  #--------------------------------------------------------
  # HEATMAP PRECISIÓN GLASSO (PARCIAL) - MEJORADO
  #--------------------------------------------------------
  output$heat_prec <- renderPlot({
    PN <- precision_norm()
    
    # Forzar diagonal a NA para mejor visualización
    diag(PN) <- NA
    
    corrplot(PN,
             method = "color",        # Cambiar a "color" en lugar de "circle"
             type = "upper",
             diag = FALSE,
             tl.cex = 1,
             tl.col = "black",
             na.label = " ",
             col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),  # Colores más contrastados
             cl.lim = c(-1, 1))       # Rango fijo
  }, height = 650, width = 650)
  #--------------------------------------------------------
  # GRAFO 3D CON PLOTLY (FILTRABLE POR NODO)
  #--------------------------------------------------------
  output$graph3D <- renderPlotly({
    PN <- precision_norm()
    fig <- create_3d_graph(PN, vars, selected_node = selected_node())
    
    if (is.null(fig)) {
      plot_ly() %>%
        layout(
          title = "No hay conexiones con este lambda. Reduce el valor.",
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE)
        )
    } else {
      fig
    }
  })
  
  # Capturar clics en nodos
  observeEvent(event_data("plotly_click", source = "graph3D"), {
    click <- event_data("plotly_click", source = "graph3D")
    if (is.null(click)) return()
    
    PN <- precision_norm()
    edges <- build_edge_list(PN, vars)
    n_edges <- nrow(edges)
    
    if (!is.null(click$curveNumber) && click$curveNumber == n_edges) {
      idx <- click$pointNumber + 1
      if (idx >= 1 && idx <= length(vars)) {
        selected_node(vars[idx])
      }
    }
  })
  
  #--------------------------------------------------------
  # TABLA CONEXIONES PERDIDAS CON MEDIADORES
  #--------------------------------------------------------
  output$tabla_perdidas <- renderTable({
    
    P <- precision_norm()
    
    # Construir lista de conexiones ACTIVAS (las que quedaron)
    edges_activas <- build_edge_list(P, vars)
    
    # Crear grafo de adyacencia para buscar caminos
    adj_list <- list()
    for (v in vars) {
      vecinos <- unique(c(
        edges_activas[from == v, to],
        edges_activas[to == v, from]
      ))
      adj_list[[v]] <- vecinos
    }
    
    # Función para encontrar mediadores (caminos de longitud 2)
    encontrar_mediadores <- function(var1, var2, adj_list) {
      if (is.null(adj_list[[var1]]) || is.null(adj_list[[var2]])) {
        return(character(0))
      }
      
      # Buscar nodos que conecten var1 y var2
      vecinos1 <- adj_list[[var1]]
      vecinos2 <- adj_list[[var2]]
      
      # Mediadores son los vecinos comunes
      mediadores <- intersect(vecinos1, vecinos2)
      
      return(mediadores)
    }
    
    # Crear tabla de correlaciones perdidas
    tabla <- data.table()
    
    for (i in 1:(nrow(P)-1)) {
      for (j in (i+1):ncol(P)) {
        perdida <- (abs(P[i,j]) < 1e-6)
        
        if (perdida) {
          var1 <- vars[i]
          var2 <- vars[j]
          
          # Buscar mediadores
          mediadores <- encontrar_mediadores(var1, var2, adj_list)
          
          # Si hay mediadores, construir string
          mediadores_str <- if (length(mediadores) > 0) {
            if (length(mediadores) <= 3) {
              paste(mediadores, collapse = ", ")
            } else {
              paste0(paste(head(mediadores, 3), collapse = ", "), 
                     " (y ", length(mediadores) - 3, " más)")
            }
          } else {
            "sin mediadores directos"
          }
          
          tabla <- rbind(tabla,
                         data.table(
                           var1 = var1,
                           var2 = var2,
                           correlacion = m_cor[i,j],
                           mediadores = mediadores_str,
                           n_mediadores = length(mediadores)
                         ))
        }
      }
    }
    
    # Ordenar por correlación absoluta
    tabla[, abs_cor := abs(correlacion)]
    tabla <- tabla[order(-abs_cor)]
    
    # Seleccionar top 20
    resultado <- head(tabla[, .(
      var1, 
      var2, 
      correlacion = round(correlacion, 4),
      mediadores,
      n_med = n_mediadores
    )], 20)
    
    return(resultado)
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
}

shinyApp(ui = ui, server = server)

