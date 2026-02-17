# Walters shiny app - Fraser Sockeye 

# 
library(shiny)
library(dplyr)
library(readxl)
library(ggplot2)

# ---- USER SETTINGS (edit once) ----
INPUT_XLSX <- "data/Walters_model_simplified.xlsx"
SHEET      <- "Sheet1"
LAG_YEARS  <- 4
FIT_YEARS  <- 1952:2019

# ---- UI ----
ui <- fluidPage(
  titlePanel("Retrospective Fraser Sockeye"),
  
  sidebarLayout(
    sidebarPanel(
      
      numericInput("retroU", "retroU (harvest rate)",
                   value = 0.3, min = 0, max = 1, step = 0.01),
      
      checkboxInput("useretro", "Use retrospective U?",
                    value = TRUE),
      
      numericInput("yrretro", "Retrospective start year",
                   value = 1990, min = 1900, max = 2100, step = 1)
    ),
    
    mainPanel(
      h4("Catch Summary"),
      verbatimTextOutput("catch_summary"),
      
      plotOutput("catch_plot", height = 300),
      plotOutput("adultreturn_plot", height = 300),
      plotOutput("retroS_plot", height = 300)
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  # ---- Load data once ----
  obs <- read_excel(
    path      = INPUT_XLSX,
    sheet     = SHEET,
    range     = "A3:I75",
    col_names = TRUE
  ) %>%
    rename(
      Year            = `Year`,
      AdultEscapement = `Sum of Adult Escapement`,
      JackEscapement  = `Sum of Jack Escapement`,
      TotalEscapement = `Sum of Total Escapement`,
      DBE             = `Sum of DBE`,
      BelowMissionC   = `Sum of Below Mission Catch (excluding Alaska catch)`,
      AboveMissionC   = `Sum of Above Mission Catch`,
      AlaskaCatch     = `Sum of Alaska Catch`,
      RunSize         = `Sum of Run Size`
    ) %>%
    mutate(Year = as.integer(Year)) %>%
    filter(!is.na(Year)) %>%
    arrange(Year)
  
  # ---- Reactive model ----
  model_out <- reactive({
    
    retroU   <- input$retroU
    useretro <- input$useretro
    yrretro  <- input$yrretro
    
    # ---- Derived observed columns ----
    obs2 <- obs %>%
      mutate(
        RunJacks = RunSize - JackEscapement,
        Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
        Ut_obs  = pmin(0.999, Catch / RunJacks),
        ENS     = pmin(1, pmax(0.0001,
                               AdultEscapement / RunJacks / (1 - Ut_obs))),
        migmort = 1 - ENS
      ) %>%
      mutate(
        AdultReturn = lead(RunJacks, n = LAG_YEARS),
        lnR_S       = log(AdultReturn / AdultEscapement)
      )
    
    # ---- Fit SR ----
    fit_df <- obs2 %>%
      filter(Year %in% FIT_YEARS) %>%
      filter(is.finite(lnR_S), is.finite(AdultEscapement))
    
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
    
    ra <- unname(coef(fit)[1])
    rb <- unname(-coef(fit)[2])
    
    obs3 <- obs2 %>%
      mutate(wt = lnR_S - (ra - rb * AdultEscapement))
    
    # ---- Retrospective vectors ----
    n <- nrow(obs3)
    
    retroR <- rep(NA_real_, n)
    retroU_vec <- if (useretro)
      ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
    else
      obs3$Ut_obs
    
    retroS <- rep(NA_real_, n)
    retroC <- rep(NA_real_, n)
    
    # Seed first lag years
    retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
    
    retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
      (1 - retroU_vec[1:LAG_YEARS]) *
      obs3$ENS[1:LAG_YEARS]
    
    retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
      retroU_vec[1:LAG_YEARS]
    
    # Forward loop
    for (i in (LAG_YEARS + 1):n) {
      j <- i - LAG_YEARS
      retroR[i] <- retroS[j] *
        exp(ra - rb * retroS[j] + obs3$wt[j])
      
      retroS[i] <- retroR[i] *
        (1 - retroU_vec[i]) *
        obs3$ENS[i]
      
      retroC[i] <- retroR[i] * retroU_vec[i]
    }
    
    out <- obs3 %>%
      mutate(
        retroR = retroR,
        retroU = retroU_vec,
        retroS = retroS,
        retroC = retroC
      )
    
    out
  })
  
  # ---- Outputs ----
  
  output$catch_summary <- renderPrint({
    out <- model_out()
    
    hist_catch  <- sum(out$Catch,  na.rm = TRUE)
    retro_catch <- sum(out$retroC, na.rm = TRUE)
    lost_Catch  <- retro_catch - hist_catch
    
    cat("Historical Catch:     ", format(round(hist_catch, 0), big.mark=","), "\n")
    cat("Retrospective Catch:  ", format(round(retro_catch, 0), big.mark=","), "\n")
    cat("Difference:    ", format(round(lost_Catch, 0), big.mark=","), "\n")
  })
  
  # catch plot 
  output$catch_plot <- renderPlot({
    out <- model_out()
    
    plot_df <- out %>%
      select(Year, Catch, retroC) %>%
      tidyr::pivot_longer(
        cols = c(Catch, retroC),
        names_to = "Type",
        values_to = "Value"
      )
    
    ggplot(plot_df, aes(x = Year, y = Value, color = Type, linetype = Type)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(
        values = c("Catch" = "black",
                   "retroC" = "firebrick"),
        labels = c("Catch" = "Observed Catch",
                   "retroC" = "Retrospective Catch")
      ) +
      labs(
        y = "Catch",
        color = "Series",
        title = "Observed vs Retrospective Catch"
      ) +
      theme_minimal() +
      theme(legend.position = "top")
  })
  
  # # lnrs plot 
  # output$lnrs_plot <- renderPlot({
  #   out <- model_out()
  #   
  #   ggplot(out, aes(x = Year, y = lnR_S)) +
  #     geom_line(linewidth = 1) +
  #     geom_point(size = 2) +
  #     labs(
  #       title = "ln(R/S) Over Time",
  #       y = "ln(R/S)"
  #     ) +
  #     theme_minimal()
  # })
  
  # adult return plot 
  output$adultreturn_plot <- renderPlot({
    out <- model_out()
    
    plot_df <- out %>%
      select(Year, AdultReturn, retroR) %>%
      tidyr::pivot_longer(
        cols = c(AdultReturn, retroR),
        names_to = "Type",
        values_to = "Value"
      )
    
    ggplot(plot_df, aes(x = Year, y = Value, color = Type)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(
        values = c("AdultReturn" = "black",
                   "retroR" = "steelblue"),
        labels = c("AdultReturn" = "Observed Adult Return",
                   "retroR" = "Retrospective Run")
      ) +
      labs(
        title = "Returns",
        y = "Number of Fish",
        color = "Series"
      ) +
      theme_minimal() +
      theme(legend.position = "top")
  })
  
  # Retro spawners plot 
  output$retroS_plot <- renderPlot({
    out <- model_out()
    
    plot_df <- out %>%
      select(Year, retroS, AdultEscapement) %>%
      tidyr::pivot_longer(
        cols = c(AdultEscapement, retroS),
        names_to = "Type",
        values_to = "Value"
      )
    
    ggplot(plot_df, aes(x = Year, y = Value, color = Type)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(
        values = c("AdultEscapement" = "black",
                   "retroS" = "steelblue"),
        labels = c("AdultEscapement" = "Observed Adult Spawners",
          "retroS" = "Retrospective Spawners")
      ) +
      labs(
        title = "Spawners",
        y = "Number of Fish",
        color = "Series"
      ) +
      theme_minimal() +
      theme(legend.position = "top")
  })
  
  output$preview <- renderTable({
    head(model_out())
  })
}

shinyApp(ui, server)

