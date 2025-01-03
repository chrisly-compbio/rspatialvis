
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
                              "Done")
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
                             max = 20, #check
                             value = 6 )#check
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
    seuratObj <- shiny::reactiveVal(NULL) # is it better to pull from global instead?
    spatialImage <- shiny::reactiveVal(NULL)
    spatialObjVarsList <- shiny::reactiveVal(NULL)

    alpha <- reactive({
      input$alpha
    })
    spotSize <- reactive({
      input$spotSize
    })

    # retrieve environment objects
    getSpatialObjs <- function() {
      objs <- ls(envir = .GlobalEnv)
      seurat_objs <- data.frame(
        objName = objs,
        objClass = sapply(objs,
                          function(x) class(get(x, envir = .GlobalEnv))[1]),
        stringsAsFactors = FALSE
      )
      seurat_objs <- seurat_objs[seurat_objs$objClass == "Seurat",]
      seurat_objs <- seurat_objs$objName
      return(seurat_objs)
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

    spatialObjData <- reactive({
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
      # rotate and flip coordinates
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

    # observe event or reactive?
    image_txt <- reactive({
      plotly::raster2uri(spatialImage())
    })

    plot <- reactive({
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
      marker = list(size = spotSize(),
                    text = ~paste0("Cell ID: ",
                                   cell_ids,
                                   "\n",
                                   "Value: ",
                                   get(input$seuratFeatureSelection))
      )
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



#
#     # old
#     shiny::observeEvent(input$seuratObjSelection, {
#       # get spatial FOVs
#       # this could be made better
#       seuratObj_val <- .GlobalEnv[[input$seuratObjSelection]]
#       seuratObj(seuratObj_val)
#       if(is.null(seuratObj_val)) {
#         return()
#       } else {
#         shiny::updateSelectInput(session, "seuratFOVSelection",
#                                  choices = names(methods::slot(seuratObj_val,
#                                                                "images")))
#       }
#     })
#
#     # spatialFOV selected
#     shiny::observeEvent(input$seuratFOVSelection, {
#       if(is.null(seuratObj())){
#         return()
#       } else {
#         spatialImage(magick::image_read(methods::slot(seuratObj(),
#                                                       "images")[[input$seuratFOVSelection]]@image))
#         seuratFeatures <- c(colnames(seuratObj()@meta.data),
#                             rownames(seuratObj()))
#         shiny::updateSelectizeInput(session, "seuratFeatureSelection",
#                                     choices = seuratFeatures,
#                                     server = TRUE)
#       }
#     })
#
#     # feature selected
#     shiny::observeEvent(input$seuratFeatureSelection, {
#       if(is.null(seuratObj())){
#         return()
#       } else {
#         # get image vars
#         img_width <- (methods::slot(seuratObj(),
#                                     "images")[[input$seuratFOVSelection]]@image |>
#                         dim())[2]
#         img_height <- (methods::slot(seuratObj(),
#                                      "images")[[input$seuratFOVSelection]]@image |>
#                          dim())[1]
#         txt_img <- plotly::raster2uri(spatialImage())
#         genes <- rownames(seuratObj())
#         cell_ids <- Seurat::Cells(seuratObj())
#         meta.data <- colnames(seuratObj()@meta.data)
#
#         # get data
#         seuratMetadata <- Seurat::GetTissueCoordinates(seuratObj(), scale="lowres")[,c("x","y")]
#         colnames(seuratMetadata) <- c("y","x") # rotate?
#         seuratMetadata["y"] <- abs(seuratMetadata["y"]-img_height) #flip
#         seuratMetadata$cell_ids <- cell_ids
#
#         # feature data
#         # will have give layer selection choice...
#         if(input$seuratFeatureSelection %in% genes ||
#              input$seuratFeatureSelection %in% meta.data) {
#           if(input$seuratFeatureSelection %in% genes) {
#             seuratMetadata[input$seuratFeatureSelection] <-
#               SeuratObject::LayerData(seuratObj(),
#                             assay = "Spatial",
#                             layer = "counts")[input$seuratFeatureSelection,]
#           } else if (input$seuratFeatureSelection %in% meta.data) {
#             seuratMetadata[input$seuratFeatureSelection] <-
#             seuratObj()@meta.data[,input$seuratFeatureSelection]
#           }
#
#           # plotly config
#           fig <- plotly::plot_ly(data = seuratMetadata,
#                                  x= ~x,
#                                  y = ~y,
#                                  text = ~paste0("Cell ID: ", cell_ids,
#                                                 "\n",
#                                                 "Value: ",
#                                                 get(input$seuratFeatureSelection)),
#                                  color = ~get(input$seuratFeatureSelection),
#                                  type = "scatter",
#                                  mode = "markers",
#                                  alpha = alpha(),
#                                  marker = list(size = spotSize(),
#                                                sizemin = 1))
#
#           # axis config
#           # tick labels break
#           xconfig <- list(
#             title = "",
#             zeroline = FALSE,
#             showline = FALSE,
#             showticklabels = FALSE,
#             showgrid = FALSE,
#             range = c(0, img_width)
#           )
#
#           yconfig <- list(
#             title = "",
#             zeroline = FALSE,
#             showline = FALSE,
#             showticklabels = FALSE,
#             showgrid = FALSE,
#             range = c(0, img_height),
#             scaleanchor="x"
#           )
#           fig <- fig |> plotly::layout(xaxis = xconfig, yaxis = yconfig)
#
#           # add image
#           fig <- fig |>
#             plotly::layout(
#               images = list(
#                 list(
#                   source = txt_img,
#                   xref = "x",
#                   yref = "y",
#                   x = 0,
#                   sizex=img_width,
#                   y=img_height,
#                   sizey=img_height,
#                   layer="below",
#                   sizing="stretch"
#                 )
#               ),
#               dragmode = "pan"
#             ) |>
#             plotly::config(scrollZoom=TRUE)
#
#           # define legend
#           if(is.numeric(unlist(seuratMetadata[input$seuratFeatureSelection]))) {
#             fig <- fig |>
#               plotly::colorbar(title = input$seuratFeatureSelection)
#           } else {
#             fig <- fig |>
#               plotly::layout(
#                 legend = list(
#                   title = list(
#                     text = input$seuratFeatureSelection)
#                   )
#                 )
#           }
#
#           # output
#           output$plot <- plotly::renderPlotly(
#             fig#() for reactive fig creation?
#           )
#         } else {
#           # do nothing if selected feature not in
#         }
#       }
#     })

    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })

  }
  shiny::runGadget(ui, server, viewer = shiny::dialogViewer("rSpatialVis",
                                                            width = 1200,
                                                            height = 800))
}

