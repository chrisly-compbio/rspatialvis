
#' A shiny visualizer for spatial objects
#' @export

spatialVis <- function() {
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Visium Visualization"),
    shiny::fillRow(flex = c(1,4),
            # controls sidebar
            shiny::fillCol(flex=c(1,2,2,9),height="800px",
              # select and refresh button (too much work to put next to selection)
              miniUI::miniContentPanel(
                shiny::actionButton("refreshSeuratObj", "Refresh")
              ),
              # select Seurat object
              miniUI::miniContentPanel(
                shiny::selectInput("seuratObjSelection",
                                   "Select Spatial Object",
                                   choices = "Please Refresh")
                # select features
                # select metadata
              ),
              # select spatial FOV
              miniUI::miniContentPanel(
                shiny::selectInput("seuratFOVSelection",
                                   "Select Spatial FOV",
                                   choices = "Please Refresh")
              ),
              miniUI::miniContentPanel(
                shiny::selectizeInput("seuratFeatureSelection",
                                      "Select Feature or Metadata",
                                      choices = "Please Refresh")
              )

            ),
            # main plot
            miniUI::miniContentPanel(
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
    shiny::observeEvent(input$seuratObjSelection, {
      # get spatial FOVs
      # this could be made better
      seuratObj_val <- .GlobalEnv[[input$seuratObjSelection]]
      seuratObj(seuratObj_val)
      if(is.null(seuratObj_val)) {
        return()
      } else {
        shiny::updateSelectInput(session, "seuratFOVSelection",
                                 choices = names(methods::slot(seuratObj_val,
                                                               "images")))
      }
    })

    # spatialFOV selected
    shiny::observeEvent(input$seuratFOVSelection, {
      if(is.null(seuratObj())){
        return()
      } else {
        spatialImage(magick::image_read(methods::slot(seuratObj(),
                                                      "images")[[input$seuratFOVSelection]]@image))
        seuratFeatures <- c(colnames(seuratObj()@meta.data),
                            rownames(seuratObj()))
        shiny::updateSelectizeInput(session, "seuratFeatureSelection",
                                    choices = seuratFeatures,
                                    server = TRUE)
      }
    })

    # feature selected
    shiny::observeEvent(input$seuratFeatureSelection, {
      if(is.null(seuratObj())){
        return()
      } else {
        # get image vars
        img_width <- (methods::slot(seuratObj(),
                                    "images")[[input$seuratFOVSelection]]@image |>
                        dim())[2]
        img_height <- (methods::slot(seuratObj(),
                                     "images")[[input$seuratFOVSelection]]@image |>
                         dim())[1]
        txt_img <- plotly::raster2uri(spatialImage())
        genes <- rownames(seuratObj())
        cell_ids <- Seurat::Cells(seuratObj())
        meta.data <- colnames(seuratObj()@meta.data)
        # this doesn't work
        marker_size <-
          methods::slot(seuratObj(),
                        "images")[[input$seuratFOVSelection]]@scale.factors$lowres *
          methods::slot(seuratObj(),
                         "images")[[input$seuratFOVSelection]]@scale.factors$spot

        # get data
        seuratMetadata <- Seurat::GetTissueCoordinates(seuratObj(), scale="lowres")[,c("x","y")]
        colnames(seuratMetadata) <- c("y","x") # rotate?
        seuratMetadata["y"] <- abs(seuratMetadata["y"]-img_height) #flip
        seuratMetadata$cell_ids <- cell_ids

        # feature data
        # will have give layer selection choice...
        if(input$seuratFeatureSelection %in% genes ||
           input$seuratFeatureSelection %in% meta.data) {
        if(input$seuratFeatureSelection %in% genes) {
          seuratMetadata[input$seuratFeatureSelection] <-
            SeuratObject::LayerData(seuratObj(),
                          assay = "Spatial",
                          layer = "counts")[input$seuratFeatureSelection,]
        } else if (input$seuratFeatureSelection %in% meta.data) {
          seuratMetadata[input$seuratFeatureSelection] <-
            seuratObj()@meta.data[,input$seuratFeatureSelection]
        }

        # plotly config
        fig <- plotly::plot_ly(data = seuratMetadata,
                               x= ~x,
                               y = ~y,
                               text = ~paste0("Cell ID: ", cell_ids,
                                              "\n",
                                              "Value: ",
                                              get(input$seuratFeatureSelection)),
                               color = ~get(input$seuratFeatureSelection),
                               type = "scatter",
                               mode = "markers",
                               marker = list(size = marker_size,
                                             sizemin = marker_size*10))

        # axis config
        # tick labels break
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
        fig <- fig |> plotly::layout(xaxis = xconfig, yaxis = yconfig)

        # add image
        fig <- fig |>
          plotly::layout(
            images = list(
              list(
                source = txt_img,
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

        # define legend
        if(is.numeric(unlist(seuratMetadata[input$seuratFeatureSelection]))) {
          fig <- fig |>
            plotly::colorbar(title = input$seuratFeatureSelection)
        } else {
          fig <- fig |>
            plotly::layout(
              legend = list(
                title = list(
                  text = input$seuratFeatureSelection)
                )
              )
        }

        # output
        output$plot <- plotly::renderPlotly(
          fig
        )
        } else {
          # do nothing if selected feature not in
        }
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

