library(shiny)
library(shinydashboard)
library(jpeg)
library(dplyr)
library(magrittr)
library(ggplot2)
library(pheatmap)
library(grid)

#Load metadata
regions <- list(c("Primary", "A -Front"),
                c("Primary", "B -Back"),
                c("Primary", "C -Side"),
                c("R1", "A -Front"),
                c("R1", "B -Back"),
                c("R1", "C -Side"),
                c("R4", "A -Front"),
                c("R4", "B -Back"),
                c("R4", "C -Side"))

LCM_coordinates <- readRDS(paste0("LCM_coordinates.rds" ))
LCM_asCN_profile_samples <- readRDS("LCM_CN_available.rds")
asCN_mtx <- readRDS(paste0("LCM_CN_mtx.rds" ))
asCN_chr_probes <- readRDS(paste0("LCM_chr_probes.rds"))
both_haplo_count_combined_by_cell <- readRDS("LCM_haplotype_info.rds")

haplo_region_map <- c("Primary" = "Primary", "R1" = "T1.1", "R4" = "T4.1")

LCM_coordinates <- LCM_coordinates %>%
  mutate(spot_status = case_when(
    CN_profile == "Yes" ~ "LCM+WGA Successful",
    barcode != ""       ~ "LCM+WGA Unsuccessful",
    TRUE                ~ "Sequencing Not Performed"
  )) %>%
  mutate(spot_status = factor(spot_status, levels = c("LCM+WGA Successful", "LCM+WGA Unsuccessful", "Sequencing Not Performed")))

spot_status_colours <- c("LCM+WGA Successful" = "springgreen4",
                         "LCM+WGA Unsuccessful" = "deeppink1",
                         "Sequencing Not Performed" = "gray50")

#Image filename lookup (loaded on demand)
LCM_image_files <- c(
  "Primary_A -Front" = "Primary_12-24POSITIONS_WRITTEN.JPG",
  "Primary_B -Back"  = "Primary_B_OV20X_cut_positions.jpg",
  "Primary_C -Side"  = "Primary_C_side_OV20x_CUT_positions.jpg",
  "R1_A -Front"      = "T1.1_1.25X_CUT_24.JPG",
  "R1_B -Back"       = "T1.1_B_Overview_cut_positions.jpg",
  "R1_C -Side"       = "T1.1_side_OV_20x_CUT_positions.jpg",
  "R4_A -Front"      = "T4.1_Front_OVERVIEW_CUT_6.3x_positions.jpeg",
  "R4_B -Back"       = "T4.1_Back_OV_20x_cut.jpg",
  "R4_C -Side"       = "T4.1_side_OV_cut_20x.jpg"
)
LCM_image_cache <- new.env(parent = emptyenv())

get_LCM_image <- function(key) {
  if (is.null(LCM_image_cache[[key]])) {
    LCM_image_cache[[key]] <- readJPEG(LCM_image_files[[key]])
  }
  LCM_image_cache[[key]]
}

#Function for plotting copy number profile
CN_states <- c("0+0" = "royalblue3",
               "1+0" = "skyblue2",
               "2+0" = "grey80",
               "1+1" = "white",
               "3+0" = "lightgoldenrod1",
               "2+1" = "khaki1",
               "4+0" = "darkorange3",
               "3+1" = "darkorange1",
               "2+2" = "orange",
               "5+0" = "red4",
               "4+1" = "red",
               "3+2" = "orangered2",
               ">5" = "purple4")
as_CN_colour <- c("royalblue3", "skyblue2", "grey80", "white", "lightgoldenrod1", "khaki1", "darkorange3", "darkorange1", "orange", "red4", "red", "orangered2", "purple4")
as_CN_breaks = c(0,10,20,21,30,31,40,41,42,50,51,52,60,1000)-0.1

plot_asCN <- function(cell_id, ascn_mtx) {
  pheatmap(mat = ascn_mtx[cell_id,,drop = F], cluster_rows = F, cluster_cols = F,
                      show_rownames = F, show_colnames = F, color = as_CN_colour, breaks = as_CN_breaks, legend = F, fontsize = 14, main = "MPNST LCM Allele Specific CN Heatmap")
}

as_CN_legend <- t(as.matrix(unlist(lapply(2:nrow(asCN_chr_probes), function(c) {
  return(rep(abs(c%%2-1), asCN_chr_probes$total.probes[c]))
}))))

plot_chr_scale <- function(asCN_chr_probes, ascn_mtx) {
  pheatmap(mat = as_CN_legend, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, color = c("white", "grey"), legend = F)
  for (k in asCN_chr_probes$chr[-1]) {
    grid.text(asCN_chr_probes$chr[k+1], x=asCN_chr_probes$cum.probes[k+1]/ncol(ascn_mtx), y=(k%%2)*0.4+.25, gp=gpar(cex = 1.5), just = "right")
  }
}
  
samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("26787_8", "30095_2", "30363_3", "30177_3", "30177_4", "30177_2")

region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")
names(region_colours) <- samples

####################################################################################################################################
#Shiny App
####################################################################################################################################
ui <- navbarPage("MPNST LCM Explorer",
                 tabPanel("About",
                          fluidRow(
                            column(width = 8, offset = 2,
                                   br(),
                                   h2("Exploring Spatial Copy Number Heterogeneity through Laser Capture Microdissection",
                                      style = "text-align: center; font-weight: bold;"),
                                   br(),
                                   div(style = "text-align: center;",
                                       imageOutput("study_design_img", height = "auto", width = "100%")
                                   ),
                                   br(),
                                   h4("About this app"),
                                   p("This interactive application accompanies the study:",
                                     strong(em("Chromosomal instability drives spatial and temporal phenotypic diversity in Schwann cancer cells.")),
                                     style = "font-size: 16px;"),
                                   p("Paper: ", a("Link to bioRxiv preprint", href = "#", target = "_blank"),
                                     style = "font-size: 16px; color: grey;"),
                                   br(),
                                   p("We performed deep multi-omic profiling of a single MPNST (malignant peripheral nerve sheath tumour) case,
                                     combining single-cell and spatial genomic and transcriptomic data to characterise tumour evolution. As part of this study, we used
                                     laser capture microdissection (LCM) to isolate individual spots from tumour tissue sections in different orientations (front, side, back) across the
                                     primary tumour and two recurrences (R1 and R4). Each spot was profiled for allele-specific copy number
                                     and haplotype imbalance.",
                                     style = "font-size: 15px;"),
                                   p("This app allows you to interactively explore the LCM data:",
                                     style = "font-size: 15px;"),
                                   tags$ul(style = "font-size: 15px;",
                                     tags$li(strong("LCM Spot tab:"), " Select a tumour region and tissue section, then click on individual
                                             LCM spots to view their allele-specific copy number profile and haplotype balance."),
                                     tags$li(strong("LCM Region tab:"), " View copy number profiles of all successfully sequenced spots within
                                             a region, and click on a genomic segment to see how that segment's copy number state is
                                             distributed spatially across tissue sections.")
                                   ),
                                   br(),
                                   div(style = "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; padding: 12px; margin-bottom: 20px;",
                                       p(icon("info-circle"), strong("Note:"),
                                         "Plots may take a few seconds to render after changing regions or clicking on spots.
                                          Please wait for the current view to fully load before making another selection,
                                          as rapid repeated clicks may cause the app to become unresponsive.",
                                         style = "font-size: 14px; margin: 0;")
                                   ),
                                   br()
                            )
                          )
                 ),
                 tabPanel("LCM Spot",
                          fluidRow(
                            column(width = 6,
                                   fluidRow(
                                     column(width = 3, selectInput("spot_region", div(style = "font-size:20px;", strong("Region:")), choices = LCM_coordinates$region %>% unique())),
                                     column(width = 3, selectInput("spot_side", div(style = "font-size:20px;", strong("Side:")), choices = LCM_coordinates$side %>% unique()))
                                   ),
                                   p("Select a tumour region and tissue side above, then click on an LCM spot on the image below to view
                                     its copy number profile and haplotype information.",
                                     style = "font-size: 14px; color: #555; margin-bottom: 10px;"),
                                   fluidRow(
                                     column(width = 6,
                                            plotOutput("LCM_plot", height = 1000, width = 1000, click = "plot1_click")
                                     )
                                   )
                            ),
                            column(width = 6,
                                   fluidRow(
                                     column(width = 10,
                                            h4("Selected spot information"),
                                            p("Sample ID, barcode, and localisation details for the selected spot.",
                                              style = "font-size: 13px; color: #555;"),
                                            verbatimTextOutput("click_info"),
                                            h4("Allele-specific copy number profile"),
                                            p("Heatmap showing the allele-specific copy number state across the genome for the selected spot.
                                              Each column represents a genomic probe, coloured by its copy number state (major + minor allele).",
                                              style = "font-size: 13px; color: #555;"),
                                            plotOutput("spot_CN_plot", height = 100, width = 1000),
                                            plotOutput("chr_scale_spot", height = 50, width = 1000),
                                            plotOutput("CN_scale_spot", height = 100, width = 1000),
                                            h4("Haplotype balance"),
                                            p("Scatter plot of haplotype 1 vs haplotype 2 read counts for all spots on this tissue section.
                                              The selected spot is highlighted. The black diagonal lines delimit the expected region for
                                              normal diploid samples \u2014 spots falling outside these lines show evidence of haplotype imbalance.",
                                              style = "font-size: 13px; color: #555;"),
                                            plotOutput("genotype_plot", height = 800, width = 800)
                                     )
                                   )
                            )
                          )
                 ),
                 tabPanel("LCM Region",
                          fluidRow(
                            column(width = 3, selectInput("region_region", div(style = "font-size:20px;", strong("Region:")), choices = LCM_coordinates$region %>% unique())),
                          ),
                          fluidRow(
                            column(width = 6,
                                   h4("Allele-specific copy number profiles for all spots in this region"),
                                   p("The heatmap below shows the allele-specific CN profile for every successfully sequenced spot in the
                                     selected region. Click on any position in the heatmap to select a genomic segment \u2014 the tissue section
                                     images below will then display the copy number state at that segment for each spot, showing how it
                                     varies spatially across the tumour.",
                                     style = "font-size: 14px; color: #555; margin-bottom: 10px;"),
                                   plotOutput("slide_CN_plot", height = 200, width = 1000,
                                              click = "plot2_click"
                                   ),
                                   plotOutput("chr_scale_1", height = 50, width = 1000,
                                   ),
                                   plotOutput("CN_scale_region", height = 100, width = 1000)
                            )
                          ),
                          fluidRow(
                            column(width = 4, div(style = "font-size:20px;", strong("Front")),
                                   plotOutput("LCM_front", height = 800, width = 800),
                            ),
                            column(width = 4, div(style = "font-size:20px;", strong("Side")),
                                   plotOutput("LCM_side", height = 800, width = 800),
                            ),
                            column(width = 4, div(style = "font-size:20px;", strong("Back")),
                                   plotOutput("LCM_back", height = 800, width = 800),
                            )
                          )
                 )
)
  
server <- function(input, output) {
  ####################################################################################################################################
  #About tab
  ####################################################################################################################################
  output$study_design_img <- renderImage({
    list(src = "study_design.png", contentType = "image/png", width = "100%", alt = "Study design overview")
  }, deleteFile = FALSE)

  ####################################################################################################################################
  #Spot tab reactive values
  ####################################################################################################################################
  LCM_coordinates_reactive <- reactive({
    LCM_coordinates %>% filter(region == input$spot_region) %>% filter(side == input$spot_side)
  })
  
  samples_on_slide_reactive <- reactive({
    LCM_coordinates_reactive() %>% filter(CN_profile == "Yes") %>% pull(sample) %>% as.character()
  })
  
  sample_to_plot <- reactive({
    as.character(nearPoints(LCM_coordinates_reactive(), input$plot1_click, addDist = TRUE, maxpoints = 1, threshold = 10)$sample)
  })
  
  ####################################################################################################################################
  #Spot tab plots
  ####################################################################################################################################
  observeEvent(input$spot_region, {
    output$spot_CN_plot <- renderPlot({})
    output$genotype_plot <- renderPlot({})
    output$click_info <- renderPrint({})
  })

  observeEvent(input$spot_side, {
    output$spot_CN_plot <- renderPlot({})
    output$genotype_plot <- renderPlot({})
    output$click_info <- renderPrint({})
  })
  
  output$LCM_plot <- renderPlot({
    LCM_coordinates_reactive() %>%
      ggplot(aes(x = coord.x, y = coord.y, color = spot_status)) +
      annotation_raster(get_LCM_image(paste0(input$spot_region, "_", input$spot_side)),
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_point(size = 10) +
      scale_color_manual(values = spot_status_colours, drop = FALSE) +
      labs(color = "Spot Status") +
      coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
      theme(text=element_text(size=16), legend.position = c(0.85, 0.9),
            legend.background = element_rect(fill = alpha("white", 0.8)),
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  })
  
  observeEvent(input$plot1_click, {
    selected <- sample_to_plot()
    
    if (length(selected) == 0) {
      output$spot_CN_plot <- renderPlot({})
      output$genotype_plot <- renderPlot({})
      output$click_info <- renderPrint({ cat("No spot selected") })
      return()
    }
    
    output$spot_CN_plot<-renderPlot({
      if (selected %in% LCM_asCN_profile_samples$sample) {
        plot_asCN(selected, asCN_mtx)
      }
    })
    
    output$genotype_plot<-renderPlot({
      both_haplo_count_combined_by_cell %>% filter(region == haplo_region_map[input$spot_region], side == input$spot_side) %>%
        mutate(selected = Barcode == selected) %>%
        ggplot(aes(x=Haplo_1_Count, y=Haplo_2_Count, colour = selected)) + geom_point() + 
        scale_x_log10(limits = c(1,4000000)) + scale_y_log10(limits = c(1,4000000)) +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black"), labels = c("TRUE" = "Selected", "FALSE" = "Other"), name = "Spot") +
        xlab("Total number of counts of haplotype 1") + ylab("Total number of counts of haplotype 2") +
        geom_abline(intercept = .5, slope = 1) +
        geom_abline(intercept = -.5, slope = 1) +
        theme(plot.title = element_text(size = 30), text=element_text(size=16), aspect.ratio=1)
    })
    
    output$click_info <- renderPrint({
      nearPoints(LCM_coordinates_reactive(), input$plot1_click, addDist = TRUE, maxpoints = 1, threshold = 10)
    })
  })
  
  output$chr_scale_spot <- renderPlot({
    plot_chr_scale(asCN_chr_probes, asCN_mtx)
  })
  
  output$CN_scale_spot <- renderPlot({
    par(mar=c(1,1,1,1))
    text(barplot(rep(1,length(as_CN_colour)), col = as_CN_colour, axes = F),
         .2, c("0+0", "1+0", "2+0", "1+1", "3+0", "2+1", "4+0", "3+1", "2+2", "5+0", "4+1", "3+2", ">5"), 0, cex =1, pos=3)
  })
  
  ####################################################################################################################################
  #Region tab reactive values
  ####################################################################################################################################
  LCM_coordinates_region_reactive <- reactive({
    LCM_coordinates %>% filter(region == input$region_region)
  })
  
  samples_on_slide_region_reactive <- reactive({
    LCM_coordinates_region_reactive() %>% filter(CN_profile == "Yes") %>% pull(sample) %>% as.character()
  })
  
  probe_position_reactive <- reactive({
    round(input$plot2_click$coords_css$x/1000*ncol(asCN_mtx),digits = 0)
  })
  
  clicked_CN_reactive <- reactive({
    LCM_coordinates_region_reactive() %>% left_join(data.frame(sample = samples_on_slide_region_reactive(),
                                                        Segment_CN = asCN_mtx[samples_on_slide_region_reactive(), probe_position_reactive()]), by = "sample") %>%
      mutate(asCN = paste0(floor(Segment_CN/10)-(Segment_CN%%10), "+", Segment_CN%%10)) %>%
      mutate(asCN_plot  = ifelse(Segment_CN >= 60 ,">5",asCN))
  })
  
  ####################################################################################################################################
  #Region tab plots
  ####################################################################################################################################
  observeEvent(input$region_region, {
    output$LCM_front <- renderPlot({})
    output$LCM_side <- renderPlot({})
    output$LCM_back <- renderPlot({})
  })
  
  observeEvent(input$plot2_click, {
    output$LCM_front <- renderPlot({
      clicked_CN_reactive() %>% filter(side == "A -Front") %>%
        ggplot(aes(x = coord.x, y = coord.y, color = asCN_plot)) +
        annotation_raster(get_LCM_image(paste0(input$region_region, "_A -Front")),
                          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        geom_point(size = 10) + scale_color_manual(values = CN_states) +
        geom_text(aes(label = asCN_plot), size = 8, nudge_x = 3, nudge_y = 3) +
        coord_cartesian(xlim = c(0,100), ylim = c(0,100)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none")
    })
    
    output$LCM_side <- renderPlot({
      clicked_CN_reactive() %>% filter(side == "C -Side") %>%
        ggplot(aes(x = coord.x, y = coord.y, color = asCN_plot)) +
        annotation_raster(get_LCM_image(paste0(input$region_region, "_C -Side")),
                          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        geom_point(size = 10) + scale_color_manual(values = CN_states) +
        geom_text(aes(label = asCN_plot), size = 8, nudge_x = 3, nudge_y = 3) +
        coord_cartesian(xlim = c(0,100), ylim = c(0,100)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none")
    })
    
    output$LCM_back <- renderPlot({
      clicked_CN_reactive() %>% filter(side == "B -Back") %>%
        ggplot(aes(x = coord.x, y = coord.y, color = asCN_plot)) +
        annotation_raster(get_LCM_image(paste0(input$region_region, "_B -Back")),
                          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        geom_point(size = 10) + scale_color_manual(values = CN_states) +
        geom_text(aes(label = asCN_plot), size = 8, nudge_x = 3, nudge_y = 3) +
        coord_cartesian(xlim = c(0,100), ylim = c(0,100)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none")
    })
  })
  
  output$slide_CN_plot <- renderPlot({
    plot_asCN(samples_on_slide_region_reactive(), asCN_mtx)
  })
  
  output$chr_scale_1 <- renderPlot({
    plot_chr_scale(asCN_chr_probes, asCN_mtx)
  })

  output$CN_scale_region <- renderPlot({
    par(mar=c(1,1,1,1))
    text(barplot(rep(1,length(as_CN_colour)), col = as_CN_colour, axes = F),
         .2, c("0+0", "1+0", "2+0", "1+1", "3+0", "2+1", "4+0", "3+1", "2+2", "5+0", "4+1", "3+2", ">5"), 0, cex =1, pos=3)
  })
}

#Run Shiny app
graphics.off()
shinyApp(ui = ui, server = server)
