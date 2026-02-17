# Chum retrospective model 
# from Carl Walters 
## January 2026

require(tidyverse)
require(readr)
require(ggplot2)

fp <- "Fraser chum Oleynik retrospective_v2.xlsx"

chum <- read_csv("data/walters_chum_retrospective_data.csv")

# ---- Pull the exact time series used by AO5:AO76 (1951–2022) ----
# Orig data row 4 is headers, data start row 5
dat <- read_excel(fp, sheet = "Orig data", skip = 3) %>%
  transmute(
    Year        = `Brood year`,
    total_stock = `total Stock, hr from Grant then Brittany`,
    U_hist      = `Historial Ut`,
    wt          = wt,
    p3          = rp3, # proportions 
    p4          = rp4,
    p5          = rp5
  ) %>%
  filter(!is.na(Year)) %>%
  arrange(Year)

# avg proportions at age -------------------- 
# # from Brittany 
# avg_prop_3 <- 0.3
# avg_prop_4 <- 0.6
# avg_prop_5 <- 0.1
# 
# check <- chum %>%
#   arrange(Year) %>%
#   mutate(age_3 = ifelse(Year < 2000, lead(total_stock, 3) * avg_prop_3,
#                         lead(total_stock,3) * lead(prop_3, 3)))
# 
# 
# # STOPPED HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# check <- chum %>%
#   arrange(Year) %>%
#   mutate(age_3 = ifelse(Year < 2000, lead(total_stock, 3) * avg_prop_3,0)) %>%
#   mutate(age_3 = case_when(Year > 2000 ~ lead(total_stock,3) * lead(prop_3, 3),
#          TRUE ~ age_3))
#          
# 
# #,
#          age_4 = lead(total_stock, 4) * avg_prop_4,
#          age_5 = lead(total_stock, 5) * avg_prop_5)


# calculate s-r params ------------------------
fit <- lm(lnrs ~ Spawners, data = chum)

aa <- unname(coef(fit)[[1]])
bb <- unname(-coef(fit)[[2]]) 

# residuals
chum <- chum %>%
  mutate(
    wt = lnrs - (aa - bb * Spawners)
  )

# ---- Simulator that RETURNS the full year-by-year table ----
simulate_policy <- function(dat, umethod = c("hist","constU","hcr"),
                            retroU = 0.5,
                            cslope = 0.7763355348064025,
                            nlim   = 487692.7340110076,
                            aa, bb,
                            util_exp = 0.6) {
  
  umethod <- match.arg(umethod)
  
  n <- nrow(dat)
  
  Run      <- numeric(n)
  U        <- numeric(n)
  Spawners <- numeric(n)
  Recruits <- numeric(n)
  Yield    <- numeric(n)
  
  for (i in seq_len(n)) {
    
    # AK: first 4 years are observed total stock, then lagged reconstruction
    if (i <= 4) {
      Run[i] <- dat$total_stock[i]
    } else {
      run_i <- 0
      if (i - 3 >= 1) run_i <- run_i + dat$p3[i - 3] * Recruits[i - 3]
      if (i - 4 >= 1) run_i <- run_i + dat$p4[i - 4] * Recruits[i - 4]
      if (i - 5 >= 1) run_i <- run_i + dat$p5[i - 5] * Recruits[i - 5]
      Run[i] <- run_i
    }
    
    # AL: policy rule for exploitation rate U
    if (umethod == "hist") {
      U[i] <- dat$U_hist[i]
    } else if (umethod == "constU") {
      U[i] <- retroU
    } else { # "hcr"
      U[i] <- max(0, cslope * (Run[i] - nlim) / Run[i])
    }
    
    # AM: spawners
    Spawners[i] <- Run[i] * (1 - U[i])
    
    # AN: recruits (Ricker with annual residual wt)
    Recruits[i] <- Spawners[i] * exp(aa - bb * Spawners[i] + dat$wt[i])
    
    # AO: yield
    Yield[i] <- Run[i] - Spawners[i]
  }
  
  tibble(
    Year = dat$Year,
    Run = Run,
    U = U,
    Spawners = Spawners,
    Recruits = Recruits,
    Yield = Yield,
    Utility_contrib = Yield^util_exp
  )
}

summarize_policy <- function(sim_tbl, util_exp = 0.6) {
  sim_tbl %>%
    summarise(
      total_yield = sum(Yield, na.rm = TRUE),
      utility     = sum(Yield^util_exp, na.rm = TRUE)
    )
}

# ---- Recreate policy compare (raw totals, not normalized) ----
util_exp <- 0.6

policies <- bind_rows(
  Historical = summarize_policy(simulate_policy(dat, "hist", aa = aa, bb = bb, util_exp = util_exp), util_exp),
  constU_0.5 = summarize_policy(simulate_policy(dat, "constU", retroU = 0.5, aa = aa, bb = bb, util_exp = util_exp), util_exp),
  MaxYield_HCR = summarize_policy(simulate_policy(dat, "hcr", cslope = 1.017, nlim = 841001, aa = aa, bb = bb, util_exp = util_exp), util_exp),
  MaxUtility_HCR = summarize_policy(simulate_policy(dat, "hcr", cslope = 0.7763355348064025, nlim = 487692.7340110076, aa = aa, bb = bb, util_exp = util_exp), util_exp),
  .id = "policy"
)

policies
