# ============================================================================
# Generate Data Dictionary (CSV & Excel) for Open-Science Repository
# ============================================================================
library(dplyr)
library(labelled)
library(openxlsx)

# Load prepared datasets (run after main R data-preparation step)
# load("D:/dataset/SMBC/Myresult11r1/analysis_datasets.RData")

build_data_dictionary <- function(df, dataset_name, custom_notes = list()) {
  vars <- names(df)
  rows <- lapply(vars, function(v) {
    col <- df[[v]]
    lbl <- labelled::var_label(col)
    lbl_str <- if (!is.null(lbl) && nzchar(as.character(lbl))) as.character(lbl) else ""

    if (is.factor(col)) {
      var_class <- "Factor (categorical)"
      coding_desc <- paste(
        sprintf("%d = %s", seq_along(levels(col)), levels(col)),
        collapse = "; "
      )
    } else if (is.numeric(col) || is.integer(col)) {
      var_class <- if (is.integer(col)) "Integer" else "Numeric (continuous)"
      rng <- range(col, na.rm = TRUE)
      coding_desc <- if (all(is.finite(rng))) {
        sprintf("Range: %.2f to %.2f", rng[1], rng[2])
      } else {
        "Numeric values"
      }
    } else {
      var_class <- paste(class(col), collapse = ", ")
      coding_desc <- "Character / other"
    }

    note_str <- if (!is.null(custom_notes[[v]])) custom_notes[[v]] else ""
    n_missing <- sum(is.na(col))
    pct_missing <- round(100 * mean(is.na(col)), 2)

    data.frame(
      dataset = dataset_name,
      variable = v,
      data_type = var_class,
      label = lbl_str,
      levels_or_range = coding_desc,
      n_missing = n_missing,
      pct_missing = pct_missing,
      description_and_units = note_str,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

# Detailed definitions and units for key study variables
variable_notes <- list(
  code = "Unique de-identified participant identifier",
  fpg = "Maternal fasting plasma glucose (mmol/L), measured at 12-20 gestational weeks",
  h1pg2 = "Maternal 1-hour post-load plasma glucose (mmol/L), measured at 20-28 gestational weeks",
  fpg_gr3 = "Tertiles of maternal fasting plasma glucose (1 = T1 [ref], 2 = T2, 3 = T3)",
  pg1h_gr3 = "Tertiles of maternal 1-hr plasma glucose (1 = T1 [ref], 2 = T2, 3 = T3)",
  fpg_gr5 = "Quintiles of maternal fasting plasma glucose (1 = Q1 [ref] to 5 = Q5)",
  pg1h_gr5 = "Quintiles of maternal 1-hr plasma glucose (1 = Q1 [ref] to 5 = Q5)",
  pregifg = "Indicator for fasting glucose >= 5.1 mmol/L (1 = Yes, 0 = No)",
  h1pg78 = "Indicator for 1-hr glucose >= 7.8 mmol/L (1 = Yes, 0 = No)",
  weightz1 = "WHO weight-for-age z-score at birth",
  weightz2 = "WHO weight-for-age z-score at 6 months",
  weightz3 = "WHO weight-for-age z-score at 1 year",
  weightz4 = "WHO weight-for-age z-score at 2 years",
  weightz5 = "WHO weight-for-age z-score at 4 years",
  weightz6 = "WHO 2007 weight-for-age z-score at 6 years",
  weightz7 = "WHO 2007 weight-for-age z-score at 8 years",
  weightz8 = "WHO 2007 weight-for-age z-score at 10 years",
  heightz1 = "WHO length-for-age z-score at birth",
  heightz2 = "WHO length-for-age z-score at 6 months",
  heightz3 = "WHO length-for-age z-score at 1 year",
  heightz4 = "WHO length-for-age z-score at 2 years",
  heightz5 = "WHO length-for-age z-score at 4 years",
  heightz6 = "WHO 2007 height-for-age z-score at 6 years",
  heightz7 = "WHO 2007 height-for-age z-score at 8 years",
  heightz8 = "WHO 2007 height-for-age z-score at 10 years",
  bmiz1 = "WHO BMI-for-age z-score at birth",
  bmiz2 = "WHO BMI-for-age z-score at 6 months",
  bmiz3 = "WHO BMI-for-age z-score at 1 year",
  bmiz4 = "WHO BMI-for-age z-score at 2 years",
  bmiz5 = "WHO BMI-for-age z-score at 4 years",
  bmiz6 = "WHO 2007 BMI-for-age z-score at 6 years",
  bmiz7 = "WHO 2007 BMI-for-age z-score at 8 years",
  bmiz8 = "WHO 2007 BMI-for-age z-score at 10 years",
  wcz1 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 1",
  wcz2 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 2",
  wcz3 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 3",
  wcz4 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 4",
  wcz5 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 5",
  wcz6 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 6",
  wcz7 = "Internally standardized (sex- and age-specific) waist circumference z-score, visit 7",
  wtz_traj_model = "Group-based trajectory membership for weight-for-age z-score",
  htz_traj_model = "Group-based trajectory membership for height/length-for-age z-score",
  bmiz_traj_model = "Group-based trajectory membership for BMI-for-age z-score",
  wcz_traj_model = "Group-based trajectory membership for waist circumference z-score",
  weightz_class_fixed = "Fixed trajectory classification for weight z-score",
  heightz_class_fixed = "Fixed trajectory classification for height z-score",
  bmiz_class_fixed = "Fixed trajectory classification for BMI z-score",
  wcz_class_fixed = "Fixed trajectory classification for waist circumference z-score",
  owo = "Childhood overweight/obesity status (Yes vs. No based on WHO BMI z-score cut-offs)",
  owoage = "Age (years) at first attainment of overweight/obesity or last valid follow-up",
  obes = "Childhood obesity status (Yes vs. No based on WHO BMI z-score cut-offs)",
  obesage = "Age (years) at first attainment of obesity or last valid follow-up",
  mage = "Maternal age at delivery (years)",
  gweek = "Gestational age at delivery (weeks)",
  bweight = "Child birth weight (kg)",
  gwg_cat = "Gestational weight gain category per Chinese WS/T 801-2022 guideline",
  brefeeding = "Exclusive breastfeeding indicator (Yes vs. No)",
  in_exposure = "Inclusion indicator in maternal glucose exposure pathway (Include vs. Exclude)"
)

dd_main <- build_data_dictionary(maindata, "maindata", variable_notes)
dd_t1   <- build_data_dictionary(data_tables1, "data_tables1", variable_notes)
dd_all  <- bind_rows(dd_main, dd_t1)

write.csv(dd_all, file.path(output_dir, "data_dictionary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

wb_dd <- createWorkbook()
addWorksheet(wb_dd, "All_Variables")
writeData(wb_dd, "All_Variables", dd_all)
addWorksheet(wb_dd, "maindata")
writeData(wb_dd, "maindata", dd_main)
addWorksheet(wb_dd, "data_tables1")
writeData(wb_dd, "data_tables1", dd_t1)
setColWidths(wb_dd, "All_Variables", cols = 1:ncol(dd_all), widths = "auto")
setColWidths(wb_dd, "maindata", cols = 1:ncol(dd_main), widths = "auto")
setColWidths(wb_dd, "data_tables1", cols = 1:ncol(dd_t1), widths = "auto")
saveWorkbook(wb_dd, file.path(output_dir, "data_dictionary.xlsx"), overwrite = TRUE)

message("Data dictionary exported to: ", output_dir)