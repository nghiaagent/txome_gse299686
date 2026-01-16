here::i_am("r/00_define_contrasts.R")

# Define contrasts
dge_deseq2_contrasts <- list(
  ## Compare between days of treatment
  ### In non-responders
  d9_vs_d4_nonresponder = c(
    "treatment_day",
    "D9_NonResponder",
    "D4_NonResponder"
  ),
  d17_vs_d4_nonresponder = c(
    "treatment_day",
    "D17_NonResponder",
    "D4_NonResponder"
  ),
  d24_vs_d4_nonresponder = c(
    "treatment_day",
    "D24_NonResponder",
    "D4_NonResponder"
  ),
  d17_vs_d9_nonresponder = c(
    "treatment_day",
    "D17_NonResponder",
    "D9_NonResponder"
  ),
  d24_vs_d9_nonresponder = c(
    "treatment_day",
    "D24_NonResponder",
    "D9_NonResponder"
  ),
  d24_vs_d17_nonresponder = c(
    "treatment_day",
    "D24_NonResponder",
    "D17_NonResponder"
  ),

  ### In responders
  d9_vs_d4_responder = c(
    "treatment_day",
    "D9_Responder",
    "D4_Responder"
  ),
  d17_vs_d4_responder = c(
    "treatment_day",
    "D17_Responder",
    "D4_Responder"
  ),
  d24_vs_d4_responder = c(
    "treatment_day",
    "D24_Responder",
    "D4_Responder"
  ),
  d17_vs_d9_responder = c(
    "treatment_day",
    "D17_Responder",
    "D9_Responder"
  ),
  d24_vs_d9_responder = c(
    "treatment_day",
    "D24_Responder",
    "D9_Responder"
  ),
  d24_vs_d17_responder = c(
    "treatment_day",
    "D24_Responder",
    "D17_Responder"
  ),

  ## Compare between responders and non-responders
  response_d4 = c(
    "treatment_day",
    "D4_Responder",
    "D4_NonResponder"
  ),
  response_d9 = c(
    "treatment_day",
    "D9_Responder",
    "D9_NonResponder"
  ),
  response_d17 = c(
    "treatment_day",
    "D17_Responder",
    "D17_NonResponder"
  ),
  response_d24 = c(
    "treatment_day",
    "D24_Responder",
    "D24_NonResponder"
  )
)

message("Sourced ", here::here("r/00_define_contrasts.R"))
