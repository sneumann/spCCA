# Shiny module for spCCA

COLOR_CHOICES <- c(
  "None"                  = "none",
  "NPC pathway"           = "npc",
  "ClassyFire superclass" = "cf"
)


mod_spcca_ui <- function(id) {
  ns <- shiny::NS(id)
  navbarPage(
    title = div(img(src='img/spcca_logo.png', alt = "spCCA", 
                    style="margin-top: -14px; padding-right:10px;padding-bottom:10px", 
                    height = 60)),
    windowTitle = "supervised penalised Canonical Correlation Analysis",
    header = tags$head(
      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = "css/ipb-styles.css"
      )
    ),
    tabPanel(
      "Data Import",
      
      shiny::sidebarPanel(
        width = 3,
        shiny::selectInput(ns("color_by"),   "Colour by",
                           choices = COLOR_CHOICES),
        shiny::tags$hr(),
      ),
      shiny::mainPanel(
        width = 9,
        shiny::tags$hr(),
        shiny::tags$h4("Here Data Import"),
      )
    ),
    tabPanel(
      "Data Analysis",
      
      shiny::sidebarPanel(
        width = 3,
        shiny::selectInput(ns("color_by"),   "Colour by",
                           choices = COLOR_CHOICES),
        shiny::tags$hr(),
      ),
      shiny::mainPanel(
        width = 9,
        shiny::tags$hr(),
        shiny::tags$h4("Here Data Analysis"),
      )
    ),
    tabPanel(
      "About",
      wellPanel(
        h4(HTML("<b>Session info</b>")),
        verbatimTextOutput(outputId = ns("rInfo"))
      )
    ),
    footer = fluidRow(
      column(
        12,
        HTML(readLines("www/ipbfooter.html"))
      )
    )
  )
}

mod_spcca_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    ## about page
    output$rInfo <- renderPrint({
      print(paste("init rInfo"))
      cat(
        R.Version()$version.string, 
        "\nspCCA App: ", spccaAppVersion, 
        "\nspCCA package: ", packageVersion,
        "\nHostname: ", system(command = "hostname", intern = TRUE)
      )
    })
  })
}
