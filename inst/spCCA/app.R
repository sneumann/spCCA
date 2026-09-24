# spCCA App -- minimal entry point. All logic lives in R/mod_spcca.R

suppressPackageStartupMessages({
  library(shiny)
})

source("version.R")
source("mod_spcca.R")


ui <- mod_spcca_ui("spcca")

server <- function(input, output, session) {
  mod_spcca_server("spcca")
}

shinyApp(ui, server)
