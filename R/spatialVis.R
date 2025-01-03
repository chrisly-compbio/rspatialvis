
#' A shiny visualizer for spatial objects
#' @export

spatialVis <- function() {
  ui <- bslib::page_sidebar(
    title = "rSpatialVis",
    sidebar = bslib::sidebar(
      bslib::accordion(
        bslib::accordion_panel(
          "Refresh Objects",
          shiny::actionButton("refreshSeuratObj",
                              "Refresh"),
          shiny::actionButton("done",
                              "Exit")
        ),
        bslib::accordion_panel(
          "Spot Properties",
          shiny::sliderInput("alpha",
                             "Alpha",
                             min = 0,
                             max = 1,
                             value = 1),
          shiny::sliderInput("spotSize",
                             "Spot size",
                             min = 1,
                             max = 30,
                             value = 6 )
        ),
        bslib::accordion_panel(
          "Object Selection",
          shiny::selectInput("seuratObjSelection",
                             "Select Spatial Object",
                             choices = "Please Refresh"),
          shiny::selectInput("seuratFOVSelection",
                             "Select Spatial FOV",
                             choices = "Please Refresh"),
          # assay/layer selection
          shiny::selectInput("seuratLayerSelection",
                              "Select Layer",
                              choices = "Please Refresh"),
          shiny::selectizeInput("seuratFeatureSelection",
                                "Select Feature or Metadata",
                                choices = "Please Refresh")
        )
      )
    ),
  # main plot
  bslib::card(
    full_screen = TRUE,
    bslib::card_body(
      plotly::plotlyOutput("plot",
                           width="100%",
                           height="100%")
    )
  )
  )

  server <- function(input, output, session) {
    # reactive vars
    seuratObj <- shiny::reactiveVal(NULL)
    spatialImage <- shiny::reactiveVal(NULL)
    spatialObjVarsList <- shiny::reactiveVal(NULL)

    alpha <- shiny::reactive({
      input$alpha
    })
    spotSize <- shiny::reactive({
      input$spotSize
    })

    # retrieve environment objects
    getSpatialObjs <- function() {
      supportedObjs <- c("Seurat")

      objs <- ls(envir = .GlobalEnv)
      spatialObjs <- data.frame(
        objName = objs,
        objClass = sapply(objs,
                          function(x) class(get(x, envir = .GlobalEnv))[1]),
        stringsAsFactors = FALSE
      )
      # add names to display object name and class in selectInput()
      spatialObjsVector <- spatialObjs$objName
      names(spatialObjsVector) <- paste0(spatialObjs$objName,
                                        " (",
                                        spatialObjs$objClass,
                                        ")")
      idx <- spatialObjs$objClass %in% supportedObjs
      spatialObjsVector <- spatialObjsVector[idx]
      return(spatialObjsVector)
    }

    # spatial object selection logic
    shiny::observeEvent(input$refreshSeuratObj, {
        shiny::updateSelectInput(session, "seuratObjSelection",
                          choices = getSpatialObjs())
    })

    # seurat object selected
    # alt pathway
    shiny::observeEvent(input$seuratObjSelection, {
      if(is.null(.GlobalEnv[[input$seuratObjSelection]])) {
        NULL #error cleanly?
      } else {
        seuratObj(.GlobalEnv[[input$seuratObjSelection]])
        spatialObjVars <- list()

        spatialObjVars[["FOVs"]] <- names(methods::slot(seuratObj(),
                                                          "images"))
        spatialObjVars[["Layers"]] <- SeuratObject::Layers(seuratObj())
        spatialObjVars[["Features"]] <- c(colnames(seuratObj()@meta.data),
                                          rownames(seuratObj()))
        spatialObjVarsList(spatialObjVars)

        # update inputs
        shiny::updateSelectInput(session, "seuratFOVSelection",
                                 choices = spatialObjVars[["FOVs"]])
        shiny::updateSelectInput(session, "seuratLayerSelection",
                                 choices = spatialObjVars[["Layers"]])
        shiny::updateSelectizeInput(session, "seuratFeatureSelection",
                                    choices = spatialObjVars[["Features"]],
                                    server = TRUE)
      }
    })

    shiny::observeEvent(input$seuratFOVSelection, {
      if(is.null(seuratObj())){
        NULL
      } else {
        spatialImage(magick::image_read(methods::slot(seuratObj(),
                                                      "images")[[input$seuratFOVSelection]]@image))
      }
    })

    # build reactive dataframe with neccessary and selected features
    spatialObjData <- shiny::reactive({
      genes <- rownames(seuratObj())
      meta.data <- colnames(seuratObj()@meta.data)
      cell_ids <- colnames(seuratObj())
      img_height <- (methods::slot(seuratObj(),
                                   "images")[[input$seuratFOVSelection]]@image |>
                       dim())[1]

      # get data
      seuratMetadata <- Seurat::GetTissueCoordinates(seuratObj(),
                                                     scale="lowres",
                                                     image = input$seuratFOVSelection)[,c("x","y")]
      # rotate and flip xy coordinates
      colnames(seuratMetadata) <- c("y","x")
      seuratMetadata["y"] <- abs(seuratMetadata["y"]-img_height)
      # add data
      seuratMetadata$cell_ids <- cell_ids
      if(input$seuratFeatureSelection %in% genes) {
        seuratMetadata[input$seuratFeatureSelection] <-
          SeuratObject::LayerData(seuratObj(),
                        assay = "Spatial",
                        layer = input$seuratLayerSelection)[input$seuratFeatureSelection,]
        } else if(input$seuratFeatureSelection %in% meta.data) {
          seuratMetadata[input$seuratFeatureSelection] <-
          seuratObj()@meta.data[,input$seuratFeatureSelection]
        }
      return(seuratMetadata)
    })

    image_txt <- shiny::reactive({
      plotly::raster2uri(spatialImage())
    })

    plot <- shiny::reactive({
      # get variables
      img_width <- (methods::slot(seuratObj(),
                                  "images")[[input$seuratFOVSelection]]@image |>
                      dim())[2]
      img_height <- (methods::slot(seuratObj(),
                                  "images")[[input$seuratFOVSelection]]@image |>
                      dim())[1]
      # config plot
      xconfig <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        range = c(0, img_width)
      )

      yconfig <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        range = c(0, img_height),
        scaleanchor="x"
      )

      plotly::plot_ly(data = spatialObjData(),
      x = ~x,
      y = ~y,
      type = "scatter",
      mode = "markers",
      alpha = input$alpha,
      color = ~get(input$seuratFeatureSelection),
      marker = list(size = spotSize()),
      text = ~paste0("Cell ID: ",
                     cell_ids,
                     "\n",
                     "Value: ",
                     get(input$seuratFeatureSelection))
      ) |>
      plotly::layout(
        xaxis = xconfig,
        yaxis = yconfig) |>
      plotly::layout(
        images = list(
          list(
            source = image_txt(),
            xref = "x",
            yref = "y",
            x = 0,
            sizex=img_width,
            y=img_height,
            sizey=img_height,
            layer="below",
            sizing="stretch"
          )
        ),
        dragmode = "pan"
      ) |>
      plotly::config(scrollZoom=TRUE)
    })

    # update legend based on feature plotted
    output$plot <- plotly::renderPlotly({
      if(is.numeric(unlist(spatialObjData()[input$seuratFeatureSelection]))) {
        plot() |>
          plotly::colorbar(title = input$seuratFeatureSelection)
      } else {
        plot() |>
          plotly::layout(
            legend = list(
              title = list(
                text = input$seuratFeatureSelection)
            )
          )
      }
    })

    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })

  }
  shiny::runGadget(ui, server, viewer = shiny::dialogViewer("rSpatialVis",
                                                            width = 1200,
                                                            height = 800))
}

