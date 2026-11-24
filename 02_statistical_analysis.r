# ============================================================================
# Main R Analysis Script:
#   - Baseline comparison tables (Table 1 & Table 2)
#   - Binary logistic & multinomial logistic regression (growth trajectories)
#   - Cox proportional hazards regression (OWO & Obesity) with Holm/BH-FDR
#   - Restricted cubic spline (RCS) dose-response curves (Cox HR & Logistic OR)
# ============================================================================
rm(list = ls(all.names = TRUE))

pkgs <- c(
  "dplyr", "tidyr", "tibble", "labelled", "tableone",
  "gtsummary", "flextable", "officer", "openxlsx",
  "nnet", "rms", "survival", "ggplot2", "patchwork"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------ 1. Directory configuration ------------------
input_dir  <- "D:/dataset/SMBC/Myresult11r1"
output_dir <- "D:/dataset/SMBC/Myresult11r1/output2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------ 2. Data import and standardization ------------------
data_tables1 <- read.csv(file.path(input_dir, "inandex1115.csv"), stringsAsFactors = FALSE)
data_raw     <- read.csv(file.path(input_dir, "tertile16.csv"), stringsAsFactors = FALSE)

data_tables1 <- data_tables1 %>%
  mutate(
    across(c("fsmoke", "fwine", "smoke1"), ~ factor(., levels = c(1, 2), labels = c("Yes", "No"))),
    parity      = factor(parity, levels = c(1, 2), labels = c("Nulliparous", "Parous")),
    across(c("medu", "fedu"), ~ factor(., levels = c(1, 2), labels = c("Under college", "College or above"))),
    income      = factor(income, levels = c(1, 2, 3), labels = c("<4000", "4000-7999", ">=8000")),
    mbmigr      = factor(mbmigr, levels = c(1, 3, 4), labels = c("<18.5", "18.5-23.9", ">=24.0")),
    delivery    = factor(delivery, levels = c(1, 2), labels = c("Vaginal", "Cesarean")),
    sex         = factor(sex, levels = c(1, 2), labels = c("Boys", "Girls")),
    in_exposure = factor(in_exposure, levels = c(0, 1), labels = c("Exclude", "Include")),
    gwg_cat     = factor(gwg_cat, levels = c(1, 2, 3), labels = c("Insufficient", "Optimal", "Excessive"))
  )

var_label(data_tables1) <- list(
  parity   = "Parity",
  mage     = "Maternal age (years)",
  gweek    = "Gestational age (weeks)",
  sex      = "Child sex",
  fedu     = "Paternal education",
  fsmoke   = "Paternal smoking",
  fwine    = "Paternal drinking",
  medu     = "Maternal education",
  smoke1   = "Maternal passive smoking before gestation",
  income   = "Household income per capita (CNY/month)",
  mbmigr   = "Maternal preconception BMI (kg/m2)",
  delivery = "Delivery method",
  gwg_cat  = "Gestational weight gain"
)

maindata <- data_raw %>%
  mutate(
    pregifg = ifelse(fpg >= 5.1, 1, 0),
    h1pg78  = factor(ifelse(h1pg2 >= 7.8, 1, 0)),
    across(c("fsmoke", "fwine", "smoke1"), ~ factor(., levels = c(1, 2), labels = c("Yes", "No"))),
    across(c("obes", "owo"), ~ factor(., levels = c(1, 0), labels = c("Yes", "No"))),
    parity     = factor(parity, levels = c(1, 2), labels = c("Nulliparous", "Parous")),
    across(c("medu", "fedu"), ~ factor(., levels = c(1, 2), labels = c("Under college", "College or above"))),
    income     = factor(income, levels = c(1, 2, 3), labels = c("<4000", "4000-7999", ">=8000")),
    mbmigr     = factor(mbmigr, levels = c(1, 3, 4), labels = c("<18.5", "18.5-23.9", ">=24.0")),
    delivery   = factor(delivery, levels = c(1, 2), labels = c("Vaginal delivery", "Cesarean section")),
    brefeeding = factor(brefeeding, levels = c(0, 1), labels = c("No", "Yes")),
    sex        = factor(sex, levels = c(1, 2), labels = c("Boys", "Girls")),
    gwg_cat    = factor(gwg_cat, levels = c(1, 2, 3), labels = c("Insufficient", "Optimal", "Excessive")),
    across(c("fpg_gr3", "fpg_gr5", "pg1h_gr3", "pg1h_gr5"), factor)
  )

data_renamed <- maindata %>%
  rename(
    "Parity"                          = "parity",
    "Maternal age (year)"             = "mage",
    "Gestational week"                = "gweek",
    "Child sex"                       = "sex",
    "Paternal education"              = "fedu",
    "Paternal smoking"                = "fsmoke",
    "Paternal drinking"               = "fwine",
    "Maternal education"              = "medu",
    "Maternal smoking"                = "smoke1",
    "Household income"                = "income",
    "Maternal BMI"                    = "mbmigr",
    "Delivery mode"                   = "delivery",
    "Birth weight (kg)"               = "bweight",
    "Breastfeeding"                   = "brefeeding",
    "Gestational weight gain"         = "gwg_cat",
    "Fasting plasma glucose (mmol/L)" = "fpg",
    "1h plasma glucose (mmol/L)"      = "h1pg2",
    "Childhood overweight/obesity"    = "owo",
    "Childhood obesity"               = "obes"
  )

save(data_tables1, maindata, data_renamed,
     file = file.path(input_dir, "analysis_datasets.RData"), compress = "xz")

# ------------------ 3. Table 1 & Table 2 export ------------------
var_t1    <- c("parity", "mage", "gweek", "fedu", "fsmoke", "fwine",
               "medu", "smoke1", "income", "mbmigr", "delivery", "sex", "gwg_cat")
catvar_t1 <- setdiff(var_t1, c("mage", "gweek"))

make_flex_tableone <- function(data, vars, factor_vars, strata_var, nonnormal_vars, footer_txt) {
  t1_obj <- CreateTableOne(vars = vars, factorVars = factor_vars, data = data,
                           strata = strata_var, addOverall = FALSE)
  t1_mat <- print(t1_obj, showAllLevels = TRUE, nonnormal = nonnormal_vars,
                  quote = FALSE, noSpaces = TRUE, varLabels = TRUE, printToggle = FALSE)
  df <- as_tibble(t1_mat, rownames = "Characteristics") %>%
    select(-test) %>%
    rename(Category = level, `P value` = p)
  df$Characteristics[1] <- "Total, n"
  for (col in names(df)) {
    df[[col]] <- as.character(df[[col]])
    Encoding(df[[col]]) <- "UTF-8"
  }
  flextable(df) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 9.5, part = "all") %>%
    bold(part = "header") %>%
    align(j = 1:2, align = "left", part = "all") %>%
    align(j = 3:ncol(df), align = "center", part = "all") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(color = "black", width = 1.5)) %>%
    hline_bottom(part = "header", border = fp_border(color = "black", width = 0.75)) %>%
    hline_bottom(part = "body", border = fp_border(color = "black", width = 1.5)) %>%
    add_footer_lines(values = footer_txt) %>%
    fontsize(size = 8.5, part = "footer") %>%
    italic(part = "footer") %>%
    autofit()
}

ft1 <- make_flex_tableone(
  data_tables1, var_t1, catvar_t1, "in_exposure", c("gweek"),
  "Values are mean (SD), median [IQR] for non-normal continuous variables, or n (%)."
)
save_as_docx("Table 1. Baseline characteristics by inclusion status." = ft1,
             path = file.path(output_dir, "Table1_Exclude_vs_Include.docx"))

df_t2 <- data_renamed %>%
  mutate(bmiz_traj_model = factor(bmiz_traj_model, levels = c(1, 2, 3),
                                  labels = c("Decreasing", "Moderate increase", "Rapid increase")))
var_t2    <- c("Parity", "Maternal age (year)", "Gestational week", "Child sex",
               "Paternal education", "Paternal smoking", "Paternal drinking",
               "Maternal education", "Maternal smoking", "Household income",
               "Maternal BMI", "Delivery mode", "Birth weight (kg)", "Breastfeeding",
               "Gestational weight gain", "Fasting plasma glucose (mmol/L)", "1h plasma glucose (mmol/L)")
catvar_t2 <- c("Parity", "Child sex", "Paternal education", "Paternal smoking",
               "Paternal drinking", "Maternal education", "Maternal smoking",
               "Household income", "Maternal BMI", "Delivery mode", "Breastfeeding",
               "Gestational weight gain")
ft2 <- make_flex_tableone(
  df_t2, var_t2, catvar_t2, "bmiz_traj_model", c("Gestational week"),
  "Values are mean (SD), median [IQR] for non-normal continuous variables, or n (%)."
)
save_as_docx("Table 2. Baseline characteristics across BMI z-score trajectories." = ft2,
             path = file.path(output_dir, "Table2_BMIz_Trajectories.docx"))

# ------------------ 4. Multi-model regression (trajectory outcomes) ------------------
format_or_ci <- function(est, se, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
          exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
}
format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
fmt_p_vec <- function(p) vapply(p, format_p, character(1))

adj_na <- function(p, method) {
  out <- rep(NA_real_, length(p))
  ok  <- !is.na(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = method)
  out
}
get_coef_any <- function(fit, term, cls = NULL) {
  if (inherits(fit, "error")) return(NULL)
  if (inherits(fit, "multinom")) {
    s <- summary(fit)
    if (!cls %in% rownames(s$coefficients) || !term %in% colnames(s$coefficients)) return(NULL)
    est <- s$coefficients[cls, term]; se <- s$standard.errors[cls, term]
    return(c(est = est, se = se, p = 2 * pnorm(-abs(est / se))))
  }
  cs <- summary(fit)$coefficients
  if (!term %in% rownames(cs)) return(NULL)
  c(est = cs[term, "Estimate"], se = cs[term, "Std. Error"], p = cs[term, "Pr(>|z|)"])
}

run_models_export <- function(dat, models_list, exposures, expo_labels,
                              bin_outcomes, multi_outcomes, outcome_labels,
                              out_cols, file_prefix, alpha = 0.05) {
  res_list <- list()
  for (m_name in names(models_list)) {
    covs <- models_list[[m_name]]
    for (x in exposures) {
      d_sub <- dat
      d_sub[[x]] <- factor(d_sub[[x]])
      ref_lev <- levels(d_sub[[x]])[1]
      nonref_levs <- levels(d_sub[[x]])[-1]
      d_sub$trend_var <- as.integer(d_sub[[x]])
      lev_prefix <- ifelse(grepl("_gr3$", x), "Tertile ", "Quintile ")
      
      add_rows <- function(y_code, fit_cat, fit_tr, cls = NULL) {
        rows <- list(
          data.frame(model = m_name, expo_code = x, expo_group = expo_labels[[x]],
                     row_no = 1L, level_text = "", outcome_code = y_code,
                     or_ci = NA_character_, p_str = NA_character_, p_num = NA_real_,
                     test_type = NA_character_, stringsAsFactors = FALSE),
          data.frame(model = m_name, expo_code = x, expo_group = expo_labels[[x]],
                     row_no = 2L, level_text = paste0(lev_prefix, ref_lev, " (ref)"),
                     outcome_code = y_code, or_ci = "1.00 (ref)", p_str = NA_character_,
                     p_num = NA_real_, test_type = NA_character_, stringsAsFactors = FALSE)
        )
        r_idx <- 3L
        for (lev in nonref_levs) {
          cf <- get_coef_any(fit_cat, paste0(x, lev), cls)
          rows[[r_idx]] <- data.frame(
            model = m_name, expo_code = x, expo_group = expo_labels[[x]],
            row_no = r_idx, level_text = paste0(lev_prefix, lev), outcome_code = y_code,
            or_ci = if (is.null(cf)) NA_character_ else format_or_ci(cf[["est"]], cf[["se"]]),
            p_str = if (is.null(cf)) NA_character_ else format_p(cf[["p"]]),
            p_num = if (is.null(cf)) NA_real_ else cf[["p"]],
            test_type = "level", stringsAsFactors = FALSE
          )
          r_idx <- r_idx + 1L
        }
        cf_tr <- get_coef_any(fit_tr, "trend_var", cls)
        rows[[r_idx]] <- data.frame(
          model = m_name, expo_code = x, expo_group = expo_labels[[x]],
          row_no = r_idx, level_text = "P for trend", outcome_code = y_code,
          or_ci = NA_character_,
          p_str = if (is.null(cf_tr)) NA_character_ else format_p(cf_tr[["p"]]),
          p_num = if (is.null(cf_tr)) NA_real_ else cf_tr[["p"]],
          test_type = "trend", stringsAsFactors = FALSE
        )
        bind_rows(rows)
      }
      
      for (y in bin_outcomes) {
        f_cat <- as.formula(paste(y, "~", paste(c(x, covs), collapse = " + ")))
        f_tr  <- as.formula(paste(y, "~", paste(c("trend_var", covs), collapse = " + ")))
        fit_cat <- tryCatch(glm(f_cat, data = d_sub, family = binomial()), error = function(e) e)
        fit_tr  <- tryCatch(glm(f_tr,  data = d_sub, family = binomial()), error = function(e) e)
        res_list[[paste(m_name, x, y, sep = "|")]] <- add_rows(y, fit_cat, fit_tr)
      }
      
      for (y in multi_outcomes) {
        d_sub[[y]] <- relevel(factor(d_sub[[y]]), ref = "0")
        nonref_cls <- setdiff(levels(d_sub[[y]]), "0")
        f_cat <- as.formula(paste(y, "~", paste(c(x, covs), collapse = " + ")))
        f_tr  <- as.formula(paste(y, "~", paste(c("trend_var", covs), collapse = " + ")))
        fit_cat <- tryCatch(nnet::multinom(f_cat, data = d_sub, trace = FALSE, maxit = 500), error = function(e) e)
        fit_tr  <- tryCatch(nnet::multinom(f_tr,  data = d_sub, trace = FALSE, maxit = 500), error = function(e) e)
        for (cl in nonref_cls) {
          res_list[[paste(m_name, x, y, cl, sep = "|")]] <- add_rows(paste0(y, "_", cl), fit_cat, fit_tr, cl)
        }
      }
    }
  }
  
  res_long <- bind_rows(res_list) %>%
    mutate(family = ifelse(is.na(test_type), NA_character_,
                           paste(test_type, ifelse(grepl("_gr3$", expo_code), "tertile", "quintile"), sep = "_"))) %>%
    group_by(model, family) %>%
    mutate(
      n_tests = ifelse(is.na(family), NA_integer_, sum(!is.na(p_num))),
      p_holm  = adj_na(p_num, "holm"),
      p_fdr   = adj_na(p_num, "BH"),
      p_fdr_str = format_p(p_fdr)
    ) %>%
    ungroup()
  
  wb <- createWorkbook()
  for (m_name in names(models_list)) {
    wide <- res_long %>%
      filter(model == m_name, outcome_code %in% out_cols) %>%
      mutate(outcome_code = factor(outcome_code, levels = out_cols)) %>%
      pivot_wider(id_cols = c(expo_code, expo_group, row_no, level_text),
                  names_from = outcome_code, values_from = c(or_ci, p_str, p_fdr_str),
                  names_glue = "{outcome_code}__{.value}")
    col_order <- unlist(lapply(out_cols, function(o) paste0(o, c("__or_ci", "__p_str", "__p_fdr_str"))))
    final_df <- wide %>%
      mutate(expo_code = factor(expo_code, levels = exposures)) %>%
      arrange(expo_code, row_no) %>%
      mutate(Variable = ifelse(row_no == 1, expo_group, level_text)) %>%
      select(Variable, all_of(col_order))
    final_df[is.na(final_df)] <- ""
    
    hdr1 <- c("Variable", unlist(lapply(out_cols, function(o) c(outcome_labels[[o]], "", ""))))
    hdr2 <- c("", rep(c(ifelse(m_name == "Crude", "Crude OR (95% CI)", "Adjusted OR (95% CI)"),
                        "P", "FDR-adjusted P"), length(out_cols)))
    addWorksheet(wb, m_name)
    writeData(wb, m_name, t(hdr1), startRow = 1, colNames = FALSE)
    writeData(wb, m_name, t(hdr2), startRow = 2, colNames = FALSE)
    writeData(wb, m_name, final_df, startRow = 3, colNames = FALSE)
    for (i in seq_along(out_cols)) mergeCells(wb, m_name, cols = (2 + 3 * (i - 1)):(1 + 3 * i), rows = 1)
  }
  
  fam_sum <- res_long %>%
    filter(!is.na(p_num)) %>%
    group_by(model, family) %>%
    summarise(n_tests = n(), raw_sig = sum(p_num < alpha),
              holm_sig = sum(p_holm < alpha), fdr_sig = sum(p_fdr < alpha), .groups = "drop")
  addWorksheet(wb, "Family_summary")
  writeData(wb, "Family_summary", as.data.frame(fam_sum))
  
  saveWorkbook(wb, file.path(output_dir, paste0(file_prefix, ".xlsx")), overwrite = TRUE)
  write.csv(res_long, file.path(output_dir, paste0(file_prefix, "_long.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

models_3 <- list(
  Crude  = character(0),
  Model1 = c("medu", "mage", "parity", "smoke1", "income", "mbmigr", "sex"),
  Model2 = c("medu", "mage", "parity", "smoke1", "income", "mbmigr", "sex",
             "brefeeding", "gwg_cat", "bweight")
)
exposures_4 <- c("fpg_gr3", "fpg_gr5", "pg1h_gr3", "pg1h_gr5")
expo_labs_4 <- c(fpg_gr3 = "Tertile of FPG", fpg_gr5 = "Quintile of FPG",
                 pg1h_gr3 = "Tertile of 1-hr PG", pg1h_gr5 = "Quintile of 1-hr PG")

dat_reg <- maindata %>%
  mutate(across(c("wtz_traj_model", "htz_traj_model",
                  "bmiz_traj_model", "wcz_traj_model"), ~ {
                    v <- as.numeric(as.character(.x))
                    ifelse(!is.na(v) & min(v, na.rm = TRUE) == 1, v - 1, v)
                  }))

out_labs_traj <- c(
  wtz_traj_model    = "Weight z-score",
  htz_traj_model    = "Height z-score",
  bmiz_traj_model_1 = "BMI z-score: Moderate increase",
  bmiz_traj_model_2 = "BMI z-score: Rapid increase",
  wcz_traj_model    = "WC z-score"
)

run_models_export(
  dat = dat_reg, models_list = models_3, exposures = exposures_4, expo_labels = expo_labs_4,
  bin_outcomes = c("wtz_traj_model", "htz_traj_model", "wcz_traj_model"),
  multi_outcomes = c("bmiz_traj_model"),
  outcome_labels = out_labs_traj,
  out_cols = c("wtz_traj_model", "htz_traj_model",
               "bmiz_traj_model_1", "bmiz_traj_model_2", "wcz_traj_model"),
  file_prefix = "Table3_3models_traj_outcomes"
)

# ------------------ 5. Cox proportional hazards regression (OWO & Obesity) ------------------
cox_outcomes <- c("owo", "obes")
outcome_time_map <- c(
  "owo"  = "owoage",
  "obes" = "obesage"
)

data_cox <- maindata %>%
  mutate(across(all_of(cox_outcomes), ~ case_when(
    as.character(.x) %in% c("Yes", "1") ~ 1,
    as.character(.x) %in% c("No", "0")  ~ 0,
    TRUE                                ~ NA_real_
  ))) %>%
  mutate(across(all_of(unname(outcome_time_map)), as.numeric))

table(data_cox$owo, useNA = "ifany")
table(data_cox$obes, useNA = "ifany")
summary(data_cox[, unname(outcome_time_map)])

trend_show_hr    <- FALSE  # TRUE = P for trend 行同时显示趋势 HR(95%CI)
show_adj_col     <- TRUE   # TRUE = 每个模型增加一列 FDR-adjusted P；FALSE = 仅 HR 和原始 P
family_by_coding <- TRUE   # TRUE = 三分位、五分位各自成族校正
alpha_cox        <- 0.05

outcome_label_cox <- c(
  "owo"  = "Childhood overweight/obesity",
  "obes" = "Childhood obesity"
)

model_label_map <- c(
  "Crude"  = "Crude Model",
  "Model1" = "Model 1",
  "Model2" = "Model 2"
)

get_level_text <- function(expo_code, lev_num) {
  if (expo_code %in% c("fpg_gr3", "pg1h_gr3")) paste0("Tertile ", lev_num)
  else paste0("Quintile ", lev_num)
}

format_hr_ci <- function(est, se, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
          exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
}

build_cox_formula <- function(time_var, status_var, x_term, covs) {
  lhs <- sprintf("Surv(%s, %s)", time_var, status_var)
  as.formula(paste(lhs, "~", paste(c(x_term, covs), collapse = " + ")))
}

get_cox_coef <- function(fit, term) {
  if (inherits(fit, "error")) return(NULL)
  cs <- summary(fit)$coefficients
  if (!term %in% rownames(cs)) return(NULL)
  c(est = cs[term, "coef"], se = cs[term, "se(coef)"], p = cs[term, "Pr(>|z|)"])
}

make_cox_rows <- function(model_name, x, y, ref_lev, levs, fit_cat, fit_trend) {
  n <- 0
  mk <- function(level_text, hr_str, p_str, p_num = NA_real_, test_type = NA_character_) {
    n <<- n + 1
    data.frame(
      model = model_name, expo_code = x, expo_group = unname(expo_labs_4[x]),
      row_no = n, level_text = level_text,
      outcome_code = y, HR_95CI = hr_str, Wald_p = p_str,
      p_num = p_num, test_type = test_type,
      stringsAsFactors = FALSE
    )
  }
  
  rows <- list(
    mk("", NA_character_, NA_character_),
    mk(paste0(get_level_text(x, ref_lev), " (ref)"), "1.00 (ref)", NA_character_)
  )
  for (lev in levs) {
    cf <- get_cox_coef(fit_cat, paste0(x, lev))
    rows[[length(rows) + 1]] <- mk(
      get_level_text(x, lev),
      if (is.null(cf)) NA_character_ else format_hr_ci(cf[["est"]], cf[["se"]]),
      if (is.null(cf)) NA_character_ else format_p(cf[["p"]]),
      if (is.null(cf)) NA_real_ else cf[["p"]],
      "level"
    )
  }
  cf <- get_cox_coef(fit_trend, "trend_var")
  rows[[length(rows) + 1]] <- mk(
    "P for trend",
    if (is.null(cf) || !trend_show_hr) NA_character_ else format_hr_ci(cf[["est"]], cf[["se"]]),
    if (is.null(cf)) NA_character_ else format_p(cf[["p"]]),
    if (is.null(cf)) NA_real_ else cf[["p"]],
    "trend"
  )
  bind_rows(rows)
}

all_cox_results <- list()
ph_checks       <- list()

for (m_name in names(models_3)) {
  covs <- models_3[[m_name]]
  message("正在拟合 Cox 模型：", m_name, " ...")
  
  for (x in exposures_4) {
    dat <- data_cox
    dat[[x]] <- factor(dat[[x]])
    ref_lev     <- levels(dat[[x]])[1]
    nonref_levs <- levels(dat[[x]])[-1]
    dat$trend_var <- as.integer(dat[[x]])
    
    for (y in cox_outcomes) {
      t_var <- unname(outcome_time_map[y])
      
      if (dplyr::n_distinct(na.omit(dat[[y]])) != 2)
        stop(paste(y, "不是二分类状态变量，请检查编码"))
      if (any(!is.na(dat[[t_var]]) & dat[[t_var]] <= 0))
        warning(paste(t_var, "存在 <= 0 的随访时间，请检查是否需要剔除或调整"))
      
      fit_cat <- tryCatch(
        coxph(build_cox_formula(t_var, y, x, covs), data = dat),
        error = function(e) e
      )
      fit_tr <- tryCatch(
        coxph(build_cox_formula(t_var, y, "trend_var", covs), data = dat),
        error = function(e) e
      )
      
      if (inherits(fit_cat, "error")) message("拟合失败：", m_name, " ", x, " ", y, " -> ", conditionMessage(fit_cat))
      if (inherits(fit_tr,  "error")) message("趋势失败：", m_name, " ", x, " ", y, " -> ", conditionMessage(fit_tr))
      
      if (!inherits(fit_cat, "error")) {
        zph <- tryCatch(cox.zph(fit_cat), error = function(e) NULL)
        if (!is.null(zph)) {
          ph_checks[[paste(m_name, x, y, sep = "_")]] <- data.frame(
            model        = m_name,
            expo_code    = x,
            outcome_code = y,
            expo_ph_p    = zph$table[x, "p"],
            global_ph_p  = zph$table["GLOBAL", "p"],
            stringsAsFactors = FALSE
          )
        }
      }
      
      all_cox_results[[paste(m_name, x, y, sep = "_")]] <-
        make_cox_rows(m_name, x, y, ref_lev, nonref_levs, fit_cat, fit_tr)
    }
  }
}

res_cox_long <- bind_rows(all_cox_results) %>%
  mutate(
    coding = ifelse(grepl("_gr3$", expo_code), "tertile", "quintile"),
    family = ifelse(is.na(test_type), NA_character_,
                    if (family_by_coding) paste(test_type, coding, sep = "_") else test_type)
  ) %>%
  group_by(model, family) %>%
  mutate(
    n_tests = ifelse(is.na(family), NA_integer_, sum(!is.na(p_num))),
    p_holm  = adj_na(p_num, "holm"),
    p_fdr   = adj_na(p_num, "BH")
  ) %>%
  ungroup() %>%
  mutate(
    p_fdr_str  = fmt_p_vec(p_fdr),
    p_holm_str = fmt_p_vec(p_holm)
  )

cox_family_summary <- res_cox_long %>%
  filter(!is.na(p_num)) %>%
  group_by(model, family) %>%
  summarise(
    n_tests     = n(),
    bonf_thresh = alpha_cox / n(),
    raw_sig     = sum(p_num  < alpha_cox),
    holm_sig    = sum(p_holm < alpha_cox),
    fdr_sig     = sum(p_fdr  < alpha_cox),
    .groups = "drop"
  )

ph_summary <- bind_rows(ph_checks)

wb_cox <- createWorkbook()
model_names <- names(models_3)
suffixes    <- c("__HR_95CI", "__Wald_p", if (show_adj_col) "__p_fdr_str")
k           <- length(suffixes)

hdr_style <- createStyle(textDecoration = "bold", halign = "center",
                         valign = "center", wrapText = TRUE, border = "bottom")
grp_style <- createStyle(textDecoration = "bold")
ctr_style <- createStyle(halign = "center")

for (y in cox_outcomes) {
  sub_df <- res_cox_long %>% filter(outcome_code == y)
  
  wide <- sub_df %>%
    mutate(model = factor(model, levels = model_names)) %>%
    pivot_wider(
      id_cols     = c(expo_code, expo_group, row_no, level_text),
      names_from  = model,
      values_from = c(HR_95CI, Wald_p, p_fdr_str),
      names_glue  = "{model}__{.value}"
    )
  
  col_order <- unlist(lapply(model_names, function(m) paste0(m, suffixes)))
  
  final <- wide %>%
    mutate(expo_code = factor(expo_code, levels = exposures_4)) %>%
    arrange(expo_code, row_no) %>%
    mutate(Variable = ifelse(row_no == 1, expo_group, level_text)) %>%
    select(Variable, all_of(col_order))
  
  body <- as.data.frame(final, stringsAsFactors = FALSE)
  body[is.na(body)] <- ""
  
  hdr1 <- c(unname(outcome_label_cox[y]),
            unlist(lapply(model_names, function(m) c(unname(model_label_map[m]), rep("", k - 1)))))
  hdr2 <- c("Variable",
            unlist(lapply(model_names, function(m) {
              hr_txt <- ifelse(m == "Crude", "Crude HR (95% CI)", "Adjusted HR (95% CI)")
              c(hr_txt, "P", if (show_adj_col) "FDR-adjusted P")
            })))
  
  sheet_name <- toupper(y)
  addWorksheet(wb_cox, sheet_name)
  writeData(wb_cox, sheet_name, t(hdr1), startRow = 1, colNames = FALSE)
  writeData(wb_cox, sheet_name, t(hdr2), startRow = 2, colNames = FALSE)
  writeData(wb_cox, sheet_name, body,    startRow = 3, colNames = FALSE)
  
  for (i in seq_along(model_names)) {
    mergeCells(wb_cox, sheet_name, cols = (2 + k * (i - 1)):(1 + k * i), rows = 1)
  }
  
  addStyle(wb_cox, sheet_name, hdr_style, rows = 1:2, cols = 1:length(hdr1),
           gridExpand = TRUE, stack = TRUE)
  grp_rows <- which(final$Variable %in% unname(expo_labs_4)) + 2
  addStyle(wb_cox, sheet_name, grp_style, rows = grp_rows, cols = 1, stack = TRUE)
  addStyle(wb_cox, sheet_name, ctr_style, rows = 3:(nrow(body) + 2),
           cols = 2:length(hdr1), gridExpand = TRUE, stack = TRUE)
  setColWidths(wb_cox, sheet_name, cols = 1, widths = 24)
  setColWidths(wb_cox, sheet_name, cols = 2:length(hdr1), widths = 18)
}

addWorksheet(wb_cox, "Family_summary")
writeData(wb_cox, "Family_summary", as.data.frame(cox_family_summary))
addStyle(wb_cox, "Family_summary", createStyle(textDecoration = "bold"),
         rows = 1, cols = 1:ncol(cox_family_summary), gridExpand = TRUE)
setColWidths(wb_cox, "Family_summary", cols = 1:ncol(cox_family_summary), widths = "auto")

if (nrow(ph_summary) > 0) {
  addWorksheet(wb_cox, "PH_test_summary")
  writeData(wb_cox, "PH_test_summary", as.data.frame(ph_summary))
  addStyle(wb_cox, "PH_test_summary", createStyle(textDecoration = "bold"),
           rows = 1, cols = 1:ncol(ph_summary), gridExpand = TRUE)
  setColWidths(wb_cox, "PH_test_summary", cols = 1:ncol(ph_summary), widths = "auto")
}

saveWorkbook(wb_cox, file.path(output_dir, "Table_owo_obes_cox_3models.xlsx"), overwrite = TRUE)
write.csv(res_cox_long, file.path(output_dir, "Table_owo_obes_cox_3models_long.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n===== Cox 模型检验族校正汇总 =====\n")
print(as.data.frame(cox_family_summary))
cat("\n===== 等比例风险假设 (PH) 检验汇总 =====\n")
print(as.data.frame(ph_summary))
message("完成：已输出 Table_owo_obes_cox_3models.xlsx 和 Table_owo_obes_cox_3models_long.csv")

# ------------------ 6. Restricted cubic spline (RCS) dose-response curves ------------------
covars_m2 <- models_3$Model2
expo_xlab <- c(fpg = "Fasting plasma glucose (mmol/L)", h1pg2 = "One hour plasma glucose (mmol/L)")

make_surv <- function(flag, age) {
  e <- ifelse(as.character(flag) == "Yes", 1, 0)
  t <- as.numeric(age)
  t[!is.na(t) & t <= 0] <- NA_real_
  list(t = t, e = e)
}
s_owo  <- make_surv(maindata$owo, maindata$owoage)
s_obes <- make_surv(maindata$obes, maindata$obesage)

dta_rcs <- maindata %>%
  mutate(OWO_t = s_owo$t, OWO_e = s_owo$e, OBES_t = s_obes$t, OBES_e = s_obes$e)

raw_traj_groups <- list(
  wtz_traj_model  = list(lab = "Weight z-score", group = "model"),
  htz_traj_model  = list(lab = "Height z-score", group = "model"),
  bmiz_traj_model = list(lab = "BMI z-score",    group = "model"),
  wcz_traj_model  = list(lab = "WC z-score",     group = "model")
)

outcome_specs <- list(
  OWO  = list(type = "cox", t = "OWO_t",  e = "OWO_e",  lab = "Overweight", group = "ts"),
  OBES = list(type = "cox", t = "OBES_t", e = "OBES_e", lab = "Obesity",    group = "ts")
)

for (y in names(raw_traj_groups)) {
  v <- as.numeric(as.character(dta_rcs[[y]]))
  if (!all(is.na(v)) && min(v, na.rm = TRUE) == 1) v <- v - 1
  dta_rcs[[y]] <- v
  uniq_v <- sort(unique(na.omit(v)))
  sp <- raw_traj_groups[[y]]
  if (length(uniq_v) == 2) {
    col_nm <- paste0(".bin_", y)
    dta_rcs[[col_nm]] <- ifelse(is.na(v), NA_real_, ifelse(v == uniq_v[2], 1, 0))
    outcome_specs[[y]] <- list(type = "logit", y = col_nm, lab = sp$lab, group = sp$group)
  } else {
    for (cl in setdiff(uniq_v, 0)) {
      spec_nm <- paste0(y, "_c", cl, "vs0")
      col_nm  <- paste0(".sub_", spec_nm)
      dta_rcs[[col_nm]] <- dplyr::case_when(is.na(v) ~ NA_real_, v == 0 ~ 0, v == cl ~ 1, TRUE ~ NA_real_)
      outcome_specs[[spec_nm]] <- list(
        type = "logit", y = col_nm,
        lab = paste0(sp$lab, " (class ", cl, " vs 0)"), group = sp$group
      )
    }
  }
}

rcs_res <- list()
for (o in names(outcome_specs)) {
  sp <- outcome_specs[[o]]
  sub_d <- if (sp$type == "cox") dta_rcs[!is.na(dta_rcs[[sp$t]]) & !is.na(dta_rcs[[sp$e]]), ]
  else dta_rcs[!is.na(dta_rcs[[sp$y]]), ]
  dd_sub <- datadist(sub_d[, c(names(expo_xlab), covars_m2)])
  options(datadist = "dd_sub")
  for (x in names(expo_xlab)) {
    rhs <- paste(c(sprintf("rcs(%s, 3)", x), covars_m2), collapse = " + ")
    fit <- if (sp$type == "cox") cph(as.formula(sprintf("Surv(%s, %s) ~ %s", sp$t, sp$e, rhs)), data = sub_d, x = TRUE, y = TRUE)
    else lrm(as.formula(sprintf("%s ~ %s", sp$y, rhs)), data = sub_d, x = TRUE, y = TRUE)
    pred <- as.data.frame(eval(bquote(Predict(fit, .(as.name(x)), ref.zero = TRUE, fun = exp))))
    an <- anova(fit)
    rn <- trimws(rownames(an))
    rcs_res[[paste(o, x, sep = "|")]] <- list(
      outcome = o, expo = x, pred = pred,
      p_overall = as.numeric(an[match(x, rn), "P"]),
      p_nonlinear = as.numeric(an[match("Nonlinear", rn), "P"]),
      n = unname(fit$stats["Obs"]),
      events = if (sp$type == "cox") unname(fit$stats["Events"]) else unname(fit$freq[length(fit$freq)])
    )
  }
}

p_tab <- bind_rows(lapply(rcs_res, function(r) {
  sp <- outcome_specs[[r$outcome]]
  data.frame(key = paste(r$outcome, r$expo, sep = "|"), Outcome_var = r$outcome,
             Group_type = sp$group, Outcome = sp$lab,
             Exposure = ifelse(r$expo == "fpg", "FPG", "1-hr PG"),
             Model = ifelse(sp$type == "cox", "Cox (HR)", "Logistic (OR)"),
             N = as.integer(r$n), Events = as.integer(r$events),
             P_overall = r$p_overall, P_nonlinear = r$p_nonlinear, stringsAsFactors = FALSE)
})) %>%
  mutate(P_overall_Holm = adj_na(P_overall, "holm"),
         P_overall_FDR  = adj_na(P_overall, "BH"),
         P_nonlinear_Holm = adj_na(P_nonlinear, "holm"),
         P_nonlinear_FDR  = adj_na(P_nonlinear, "BH"))

wb_rcs <- createWorkbook()
addWorksheet(wb_rcs, "RCS_P_values")
writeData(wb_rcs, "RCS_P_values", as.data.frame(select(p_tab, -key)))
saveWorkbook(wb_rcs, file.path(output_dir, "RCS_P_values.xlsx"), overwrite = TRUE)

plot_one_rcs <- function(r, x_lab, y_lab, panel_tag) {
  dat <- r$pred[order(r$pred[[r$expo]]), ]
  p_txt <- paste0("P-overall = ", format_p(r$p_overall), "\nP-nonlinear = ", format_p(r$p_nonlinear))
  ggplot(dat, aes(x = .data[[r$expo]], y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#6baed6", alpha = 0.25) +
    geom_line(color = "#0068b5", linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
    annotate("text", x = Inf, y = Inf, label = p_txt, hjust = 1.05, vjust = 1.4, size = 3.2) +
    annotate("text", x = -Inf, y = Inf, label = panel_tag, hjust = -0.3, vjust = 1.3, size = 6, fontface = "bold") +
    scale_x_continuous(labels = function(v) sprintf("%.1f", v), expand = expansion(mult = c(0.02, 0.02))) +
    theme_classic(base_size = 12.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(8, 10, 8, 12)) +
    labs(x = x_lab, y = y_lab)
}

save_rcs_group <- function(grp_name, file_stem, h_per_row, outer_only = TRUE) {
  outs <- names(Filter(function(s) s$group == grp_name, outcome_specs))
  ps <- list(); idx <- 0L
  for (i in seq_along(outs)) {
    for (x in names(expo_xlab)) {
      idx <- idx + 1L
      r <- rcs_res[[paste(outs[i], x, sep = "|")]]
      sp <- outcome_specs[[outs[i]]]
      xl <- if (!outer_only || i == length(outs)) expo_xlab[[x]] else ""
      yl <- if (!outer_only || x == names(expo_xlab)[1]) {
        if (sp$type == "cox") paste0(sp$lab, " risk, HR (95% CI)") else paste0(sp$lab, "\nOR (95% CI)")
      } else ""
      ps[[idx]] <- plot_one_rcs(r, xl, yl, LETTERS[idx])
    }
  }
  fig <- wrap_plots(ps, ncol = length(expo_xlab), byrow = TRUE)
  ggsave(file.path(output_dir, paste0(file_stem, ".tiff")), fig, device = "tiff",
         width = 12, height = h_per_row * length(outs), units = "in", dpi = 300, compression = "lzw")
  ggsave(file.path(output_dir, paste0(file_stem, ".pdf")), fig,
         width = 12, height = h_per_row * length(outs), units = "in")
}

save_rcs_group("ts",    "RCS_curve_OWO_OBES",       5.0, outer_only = FALSE)
save_rcs_group("model", "RCS_curve_traj_model_all", 3.8, outer_only = TRUE)