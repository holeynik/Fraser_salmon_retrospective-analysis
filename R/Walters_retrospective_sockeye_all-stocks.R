# Sockeye retrospective model from Carl Walters excel model 
# January 2026
# Haley Oleynik 

#   * Observed exploitation (Ut) and ENS
#   * Fitting ln(R/S) = ra - rb * S
#   * Retrospective predictionwith 4-year lag using historical residuals

# load packages 
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})


# Configure  -----------------------------

# Years used in the Excel regression ranges:
FIT_YEARS  <- 1952:2019

# Lag in rows for R/S in the spreadsheet (X4 = O8, i.e., +4 rows)
LAG_YEARS  <- 4

# Read observed data -----------------------------

# new read observed data with all stocks 
obs <- read_csv("Data/Walters_model_all-stocks.csv")

obs <- obs %>%
  rename(
    Year            = `Year`,
    AdultEscapement = `Adult Escapement`,
    JackEscapement  = `Jack Escapement`,
    TotalEscapement = `Total Escapement`,
    DBE             = `DBE`,
    BelowMissionC   = `Below Mission Catch`,
    AboveMissionC   = `Above Mission Catch`,
    AlaskaCatch     = `Alaska Catch`,
    RunSize         = `Run Size`
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year)) %>%
  arrange(Stock,Year)

# Retrospective model function ------------------------------------
run_retro_model <- function(dat,
                            retroU = 0.3,
                            useretro = 1,
                            yrretro = 1990,
                            FIT_YEARS = 1952:2019,
                            LAG_YEARS = 4) {
  
  dat <- dat %>% arrange(Year)
  
  # ---- Observed derived columns ----
  obs2 <- dat %>%
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
  
  # ---- Fit Ricker ----
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS) %>%
    filter(is.finite(lnR_S), is.finite(AdultEscapement))
  
  fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
  
  ra <- coef(fit)[1]
  rb <- -coef(fit)[2]
  
  obs3 <- obs2 %>%
    mutate(wt = lnR_S - (ra - rb * AdultEscapement))
  
  n <- nrow(obs3)
  
  retroR <- rep(NA_real_, n)
  retroU_vec <- if (useretro)
    ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
  else
    obs3$Ut_obs
  
  retroS <- rep(NA_real_, n)
  retroC <- rep(NA_real_, n)
  
  # ---- Seed ----
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    (1 - retroU_vec[1:LAG_YEARS]) *
    obs3$ENS[1:LAG_YEARS]
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    retroU_vec[1:LAG_YEARS]
  
  # ---- Forward loop ----
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retroR[i] <- retroS[j] *
      exp(ra - rb * retroS[j] + obs3$wt[j])
    retroS[i] <- retroR[i] *
      (1 - retroU_vec[i]) *
      obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(
      retroR = retroR,
      retroU = retroU_vec,
      retroS = retroS,
      retroC = retroC,
      ra = ra,
      rb = rb
    )
}

# Run retrospective -----------------------------------------------------------
out_all <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 retroU = 0.3,
                                 useretro = 1,
                                 yrretro = 1990)) %>%
  ungroup()

# Get stock summaries ---------------------------------------------------------
stock_summary <- out_all %>%
  group_by(Stock) %>%
  summarise(
    hist_catch  = sum(Catch,  na.rm = TRUE),
    retro_catch = sum(retroC, na.rm = TRUE),
    lost_catch  = retro_catch - hist_catch,
    ra = first(ra),
    rb = first(rb)
  )

