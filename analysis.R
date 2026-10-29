library(data.table)
library(fixest)
library(minpack.lm)
library(MatchIt)
library(ggplot2)
library(dplyr)
library(tidyr)

# SETTINGS
set.seed(42)
setDTthreads(1)
setFixest_nthreads(1)

DATA_PATH <- "climate_ai_filtered.csv"
OUT_DIR <- "outputs_wedge_complete_cases"

WIN_PROB <- 0.01
PEAKRATE_SMOOTH_K <- 3
PEAKRATE_MIN_YEARS <- 4
RUN_MATCHING <- TRUE
RUN_FULL_CPC <- TRUE
MIN_TOPICS <- 10L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SSC <- fixest::ssc(K.adj = TRUE, K.fixef = "nonnested", G.adj = TRUE, G.df = "min", t.df = "min", K.exact = FALSE)
write_csv <- function(x, name) data.table::fwrite(x, file.path(OUT_DIR, name), na = "NA")

# DATA PREPARATION
required_columns <- c(
  "patent_id", "relevance", "topic", "grant_year", "cpc_sections",
  "AI", "AI_30", "AI_40", "AI_50", "AI_60", "AI_70", "ai_score_max",
  "ai_score_ml", "ai_score_evo", "ai_score_nlp", "ai_score_speech",
  "ai_score_vision", "ai_score_kr", "ai_score_planning", "ai_score_hardware",
  "predict50_ml", "predict50_evo", "predict50_nlp", "predict50_speech",
  "predict50_vision", "predict50_kr", "predict50_planning", "predict50_hardware",
  "ln_xi_real", "ln_cites_3y", "ln_cites_5y",
  "ln_examiner_forward_cites_3y", "ln_examiner_forward_cites_5y",
  "ln_num_claims", "n_cpc_groups", "ln_team_size", "female_inventor_share",
  "ln_n_assignees", "has_company_assignee", "has_government_assignee",
  "n_applicant_backward_citations", "n_examiner_backward_citations",
  "n_applicant_science_refs", "mean_science_confscore", "has_science_ref",
  "primary_assignee_id",
  "self_cite_share", "cross_topic_share", "early_share_3of5",
  "ln_total_cites_edges", "ln_cites_0_3_edges", "ln_cites_3_5_edges",
  "ln_self_cites", "ln_external_cites", "ln_citing_assignee_breadth",
  "ln_citing_topic_breadth", "ln_cross_topic_cites", "ln_within_topic_cites"
)

raw_names <- names(data.table::fread(DATA_PATH, nrows = 0L, showProgress = FALSE))
missing_columns <- setdiff(required_columns, raw_names)
if (length(missing_columns)) stop("Missing input columns: ", paste(missing_columns, collapse = ", "))
base <- as.data.frame(data.table::fread(DATA_PATH, select = required_columns,
  colClasses = list(character = c("patent_id", "topic", "cpc_sections", "primary_assignee_id")),
  showProgress = FALSE))
if (anyNA(base$patent_id) || any(!nzchar(trimws(base$patent_id))) || anyDuplicated(base$patent_id)) stop("patent_id must be nonmissing and unique; resolve duplicates upstream.")
character_vars <- c("patent_id", "topic", "cpc_sections", "primary_assignee_id")
for (v in setdiff(names(base), character_vars)) {
  original <- base[[v]]
  base[[v]] <- suppressWarnings(as.numeric(original))
  if (any(!is.na(original) & is.na(base[[v]]))) stop("Nonnumeric values in ", v)
  base[[v]][!is.finite(base[[v]])] <- NA_real_
}
for (v in c("relevance", "AI", paste0("AI_", c(30, 40, 50, 60, 70)), grep("^predict50_", names(base), value = TRUE),
  "has_company_assignee", "has_government_assignee", "has_science_ref")) {
  if (any(!is.na(base[[v]]) & !base[[v]] %in% c(0, 1))) stop("Expected binary variable: ", v)
}
if (any(!is.na(base$grant_year) & base$grant_year != floor(base$grant_year))) stop("grant_year must be integer-valued.")
base$topic <- trimws(base$topic)
base$topic[!is.na(base$topic) & (base$topic == "-1" | base$topic == "")] <- NA_character_
base$topic <- factor(base$topic)
base$cpc_section <- substr(base$cpc_sections, 1L, 1L)
count_map <- c(ln_cpc_breadth = "n_cpc_groups", ln_applicant_bwd = "n_applicant_backward_citations",
  ln_examiner_bwd = "n_examiner_backward_citations", ln_applicant_sci = "n_applicant_science_refs")
for (v in names(count_map)) {
  x <- base[[count_map[[v]]]]
  if (any(x < 0, na.rm = TRUE)) stop("Negative counts in ", count_map[[v]])
  base[[v]] <- log1p(x)
}
base$mean_science_conf <- base$mean_science_confscore
structural_zero <- !is.na(base$has_science_ref) & base$has_science_ref == 0 & is.na(base$mean_science_conf)
base$mean_science_conf[structural_zero] <- 0
base$AI_bin <- base$AI
base$ai_intensity <- base$ai_score_max
controls <- list(A = c("ln_num_claims", "ln_cpc_breadth"),
  B = c("ln_team_size", "female_inventor_share", "ln_n_assignees", "has_company_assignee", "has_government_assignee"),
  C = c("ln_applicant_bwd", "ln_examiner_bwd"),
  D = c("ln_applicant_sci", "mean_science_conf", "has_science_ref"))
control_vars <- function(block = "ABCD") unname(unlist(controls[strsplit(block, "", fixed = TRUE)[[1]]]))
valid <- function(d, vars) {
  if (length(setdiff(vars, names(d)))) stop("Missing estimation variables: ", paste(setdiff(vars, names(d)), collapse = ", "))
  keep <- rep(TRUE, nrow(d))
  for (v in vars) {
    x <- d[[v]]
    keep <- keep & if (is.numeric(x)) is.finite(x) else !is.na(x) & nzchar(trimws(as.character(x)))
  }
  keep
}
write_csv(data.frame(variable = names(base), missing = vapply(base, function(x) sum(is.na(x)), numeric(1))), "input_missingness.csv")
# MATURITY
winsorise <- function(x, p = 0.01) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

.ma_centered <- function(x, k) {
  n <- length(x)
  if (n < 2 || k <= 1) return(x)
  half <- (k - 1) %/% 2
  vapply(seq_len(n), function(i) {
    mean(x[max(1, i - half):min(n, i + half)], na.rm = TRUE)
  }, numeric(1))
}

compute_maturity_measures <- function(data, smooth_k = PEAKRATE_SMOOTH_K,
                                      min_years_peak = PEAKRATE_MIN_YEARS) {
  by_ty <- data %>%
    count(topic, grant_year, name = "n_patents") %>%
    arrange(topic, grant_year) %>%
    group_by(topic) %>%
    mutate(
      cumulative = cumsum(n_patents),
      total_t = sum(n_patents),
      cumshare = ifelse(total_t > 0, cumulative / total_t, NA_real_),
      t_norm = (grant_year - min(grant_year)) /
        max(1, max(grant_year) - min(grant_year))
    ) %>%
    ungroup()

  fit_logistic_topic <- function(td) {
    t <- td$t_norm
    C <- td$cumulative
    if (length(unique(C[is.finite(C)])) < 4) {
      return(td %>% mutate(maturity = cumulative / max(cumulative)))
    }
    tryCatch({
      fit <- nlsLM(
        C ~ K / (1 + exp(-r * (t - t0))),
        data = td,
        start = list(
          K = max(C) * 1.5,
          r = 5,
          t0 = t[which.min(abs(C - max(C) / 2))]
        ),
        lower = c(K = max(C) * 1.05, r = 0.01, t0 = -1),
        upper = c(K = Inf, r = 50, t0 = 2),
        control = nls.lm.control(maxiter = 200)
      )
      Khat <- as.numeric(coef(fit)["K"])
      td %>% mutate(maturity = pmin(pmax(as.numeric(predict(fit)) / Khat, 0), 1))
    }, error = function(e) {
      td %>% mutate(maturity = cumulative / max(cumulative))
    })
  }

  logistic <- by_ty %>%
    group_by(topic) %>%
    group_modify(~ fit_logistic_topic(.x)) %>%
    ungroup() %>%
    transmute(topic, grant_year, mat_logistic = maturity)

  cumshare <- by_ty %>%
    transmute(topic, grant_year, mat_cumshare = cumshare)

  peakrate <- by_ty %>%
    group_by(topic) %>%
    mutate(
      n_years = dplyr::n(),
      flow_smooth = .ma_centered(n_patents, smooth_k),
      peak_year = grant_year[which.max(flow_smooth)],
      span = max(1, max(grant_year) - min(grant_year)),
      mat_peakrate = if_else(
        n_years < min_years_peak,
        (grant_year - min(grant_year)) / span,
        pmin(pmax(0.5 + (grant_year - peak_year) / (2 * span), 0), 1)
      )
    ) %>%
    ungroup() %>%
    transmute(topic, grant_year, mat_peakrate)

  logistic %>%
    full_join(cumshare, by = c("topic", "grant_year")) %>%
    full_join(peakrate, by = c("topic", "grant_year")) %>%
    arrange(topic, grant_year)
}

climate_rows <- !is.na(base$relevance) & base$relevance == 1
sample_log <- list()
scale_logbook <- list()
prepare <- function(citation, window, climate = TRUE, ids = NULL, max_year = NULL, label = NULL) {
  if (is.null(label)) label <- paste(if (climate) "climate" else "full", citation, sep = "_")
  available <- !is.na(base$relevance) & !is.na(base$ln_xi_real) & !is.na(base[[citation]])
  d <- base[available, , drop = FALSE]
  value_limits <- quantile(d$ln_xi_real, c(WIN_PROB, 1 - WIN_PROB), na.rm = TRUE, names = FALSE)
  cite_limits <- quantile(d[[citation]], c(WIN_PROB, 1 - WIN_PROB), na.rm = TRUE, names = FALSE)
  d$value_w <- winsorise(d$ln_xi_real, WIN_PROB)
  d$cites_w <- winsorise(d[[citation]], WIN_PROB)
  if (climate) d <- d[d$relevance == 1, , drop = FALSE]
  d <- d[valid(d, c("topic", "grant_year")), , drop = FALSE]
  if (!is.null(ids)) d <- d[d$patent_id %in% ids, , drop = FALSE]
  if (!is.null(max_year)) d <- d[!is.na(d$grant_year) & d$grant_year <= max_year, , drop = FALSE]
  population_n <- nrow(d)
  if (!population_n) stop("Empty preprocessing sample: ", label)
  maturity <- compute_maturity_measures(d[, c("topic", "grant_year"), drop = FALSE])
  write_csv(maturity, paste0("maturity_", label, ".csv"))
  idx <- match(paste(d$topic, d$grant_year), paste(maturity$topic, maturity$grant_year))
  for (v in c("mat_logistic", "mat_cumshare", "mat_peakrate")) d[[v]] <- maturity[[v]][idx]
  d$mat <- d$mat_logistic
  value_mean <- mean(d$value_w, na.rm = TRUE)
  cite_mean <- mean(d$cites_w, na.rm = TRUE)
  value_sd <- sd(d$value_w, na.rm = TRUE)
  cite_sd <- sd(d$cites_w, na.rm = TRUE)
  if (!is.finite(value_sd) || !is.finite(cite_sd) || value_sd <= 0 || cite_sd <= 0) stop("Invalid outcome scale: ", label)
  d$y0 <- (d$value_w - value_mean) / value_sd
  d$y1 <- (d$cites_w - cite_mean) / cite_sd
  d$gap <- d$y1 - d$y0
  d$cites_raw <- d[[citation]]
  center <- mean(d$ai_intensity, na.rm = TRUE)
  d$ai_int_c <- d$ai_intensity - center
  needed <- c("patent_id", "topic", "grant_year", "cpc_section", "y0", "y1", "mat", "AI_bin", control_vars())
  keep <- valid(d, needed)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("No complete cases in ", label)
  scale_logbook[[label]] <<- data.frame(reference = label, outcome = c("private_value", citation),
    winsorization_population = sum(available), standardization_population = population_n,
    lower = c(value_limits[1], cite_limits[1]), upper = c(value_limits[2], cite_limits[2]),
    mean = c(value_mean, cite_mean), sd = c(value_sd, cite_sd))
  sample_log[[label]] <<- data.frame(sample = label, window = window, population = population_n,
    complete = sum(keep), excluded_missing = sum(!keep), topics = length(unique(d$topic)), ai_intensity_center = center)
  d
}
d5 <- prepare("ln_cites_5y", 5L)
write_csv(data.frame(patents = nrow(d5), stacked_observations = 2L * nrow(d5),
  topics = length(unique(d5$topic)), structural_science_confidence_zeros = sum(structural_zero)), "baseline_sample_size.csv")
message("Baseline: ", nrow(d5), " patents; ", 2L * nrow(d5), " stacked observations.")
if (length(unique(d5$topic)) < MIN_TOPICS) warning("Fewer than ", MIN_TOPICS, " topics; conventional clustered inference may be unreliable.")

# ESTIMATION
fe_specs <- c(FE1 = "grant_year^out", FE2 = "grant_year^out + cpc_section^out",
  FE3 = "grant_year^out + topic^out", FE5 = "patent_id + topic^out + grant_year^out")
canon <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(z) paste(sort(z), collapse = ":"), character(1))
model_log <- list()
fit_model <- function(fml, d, label, cluster = "topic", expected = character()) {
  f <- as.formula(fml)
  vars <- unique(c("patent_id", "topic", all.vars(f), strsplit(cluster, " + ", fixed = TRUE)[[1]]))
  dat <- d[valid(d, vars), vars, drop = FALSE]
  if (!nrow(dat)) stop("Empty sample: ", label)
  if ("out" %in% names(dat) && (nrow(dat) != 2L * length(unique(dat$patent_id)) || any(table(dat$patent_id) != 2L))) stop("Unbalanced outcome pairs: ", label)
  g <- length(unique(dat$topic))
  if (g < 3L) stop("Fewer than three topics: ", label)
  m <- fixest::feols(f, data = dat, cluster = as.formula(paste("~", cluster)), ssc = SSC,
    lean = TRUE, mem.clean = TRUE, nthreads = 1, fixef.tol = 1e-8, fixef.iter = 20000, fixef.rm = "none")
  if (!is.null(m$convStatus) && !isTRUE(m$convStatus)) stop("Estimation failed to converge: ", label)
  ct <- fixest::coeftable(m)
  terms <- canon(rownames(ct))
  absent <- setdiff(canon(expected), terms)
  if (length(absent)) stop("Required coefficients are not identified in ", label, ": ", paste(absent, collapse = ", "))
  V <- stats::vcov(m)
  dimnames(V) <- list(canon(rownames(V)), canon(colnames(V)))
  df <- as.numeric(fixest::degrees_freedom(m, type = "t"))
  tab <- data.frame(term = terms, estimate = ct[, 1], std_error = ct[, 2], p_value = ct[, 4], row.names = NULL)
  tab$conf_low <- tab$estimate - qt(0.975, df) * tab$std_error
  tab$conf_high <- tab$estimate + qt(0.975, df) * tab$std_error
  ans <- list(coefs = tab, V = V, n = stats::nobs(m), patents = length(unique(dat$patent_id)), topics = g,
    df = df, cluster = cluster, formula = fml, dropped = paste(m$collin.var, collapse = "; "))
  model_log[[label]] <<- data.frame(model = label, formula = fml, observations = ans$n, patents = ans$patents,
    topics = g, inference_df = df, cluster = cluster, dropped_terms = ans$dropped)
  rm(m, dat)
  gc(FALSE)
  ans
}
stack <- function(d, vars) {
  vars <- unique(c("patent_id", "topic", "grant_year", vars))
  d <- d[valid(d, c(vars, "y0", "y1")), c(vars, "y0", "y1"), drop = FALSE]
  a <- d[, vars, drop = FALSE]
  b <- a
  a$out <- 0L
  b$out <- 1L
  a$outcome_std <- d$y0
  b$outcome_std <- d$y1
  rbind(a, b)
}
fit_stack <- function(d, label, fe = "FE5", block = "ABCD", exposure = "AI_bin", cluster = "topic") {
  cv <- control_vars(block)
  terms <- c("mat", exposure, paste0("mat:", exposure), cv)
  rhs <- if (fe == "FE5") paste(paste0(terms, ":out"), collapse = " + ") else paste0("out * (", paste(terms, collapse = " + "), ")")
  vars <- c("mat", exposure, cv, if (fe == "FE2") "cpc_section", strsplit(cluster, " + ", fixed = TRUE)[[1]])
  dat <- stack(d, vars)
  fit_model(paste("outcome_std ~", rhs, "|", fe_specs[[fe]]), dat, label, cluster,
    expected = paste0("mat:", exposure, ":out"))
}
fit_single <- function(d, y, label, block = "ABCD", exposure = "AI_bin", cluster = "topic") {
  rhs <- paste(c("mat", exposure, paste0("mat:", exposure), control_vars(block)), collapse = " + ")
  fit_model(paste(y, "~", rhs, "| topic + grant_year"), d, label, cluster, paste0("mat:", exposure))
}
extract <- function(m, term) {
  i <- match(canon(term), m$coefs$term)
  if (is.na(i)) stop("Missing coefficient: ", term)
  m$coefs[i, , drop = FALSE]
}
tex_escape <- function(x) {
  x <- gsub("_", "\\_", as.character(x), fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x
}
export_models <- function(models, file, labels = NULL, note = "") {
  long <- do.call(rbind, lapply(names(models), function(nm) {
    m <- models[[nm]]
    z <- m$coefs
    z$model <- nm
    z$observations <- m$n
    z$patents <- m$patents
    z$topics <- m$topics
    z$inference_df <- m$df
    z$cluster <- m$cluster
    z
  }))
  write_csv(long, paste0(file, ".csv"))
  saveRDS(models, file.path(OUT_DIR, paste0(file, ".rds")))
  if (is.null(labels)) {
    ts <- unique(long$term)
    labels <- setNames(gsub(":", " x ", ts, fixed = TRUE), ts)
  }
  keys <- canon(names(labels))
  lines <- c("\\begin{table}[!htbp]", "\\centering", "\\begin{adjustbox}{max width=\\textwidth}",
    paste0("\\begin{tabular}{l", paste(rep("c", length(models)), collapse = ""), "}"),
    "\\toprule", paste0(" & ", paste(tex_escape(names(models)), collapse = " & "), " \\\\"), "\\midrule")
  for (i in seq_along(keys)) {
    vals <- ses <- rep("", length(models))
    for (j in seq_along(models)) {
      r <- match(keys[i], models[[j]]$coefs$term)
      if (!is.na(r)) {
        z <- models[[j]]$coefs[r, ]
        star <- if (is.na(z$p_value)) "" else if (z$p_value < .01) "***" else if (z$p_value < .05) "**" else if (z$p_value < .1) "*" else ""
        vals[j] <- paste0("$", sprintf("%.4f", z$estimate), if (nzchar(star)) paste0("^{", star, "}") else "", "$")
        ses[j] <- paste0("$(", sprintf("%.4f", z$std_error), ")$")
      }
    }
    lines <- c(lines, paste0(tex_escape(unname(labels[i])), " & ", paste(vals, collapse = " & "), " \\\\"),
      paste0(" & ", paste(ses, collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\midrule")
  metadata <- c(n = "Observations", patents = "Patents", topics = "Topics", cluster = "Clustering")
  for (field in names(metadata)) {
    values <- vapply(models, function(m) tex_escape(m[[field]]), character(1))
    lines <- c(lines, paste0(metadata[[field]], " & ", paste(values, collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{adjustbox}",
    paste0("\\caption{", tex_escape(note), " Standard errors are in parentheses. Significance: * p<0.10, ** p<0.05, *** p<0.01.}"),
    paste0("\\label{tab:", file, "}"), "\\end{table}")
  writeLines(lines, file.path(OUT_DIR, paste0(file, ".tex")))
  print(long[, c("model", "term", "estimate", "std_error", "p_value", "observations")], row.names = FALSE)
  invisible(long)
}
export_data <- function(d, file, caption) {
  write_csv(d, paste0(file, ".csv"))
  format_cell <- function(x) {
    if (is.numeric(x)) return(ifelse(is.na(x), "", formatC(x, digits = 4, format = "fg", flag = "#")))
    tex_escape(x)
  }
  cells <- as.data.frame(lapply(d, format_cell), stringsAsFactors = FALSE)
  rows <- apply(cells, 1L, function(x) paste0(paste(x, collapse = " & "), " \\\\"))
  writeLines(c("\\begin{table}[!htbp]", "\\centering", "\\begin{adjustbox}{max width=\\textwidth}",
    paste0("\\begin{tabular}{", paste(rep("l", ncol(d)), collapse = ""), "}"), "\\toprule",
    paste0(paste(tex_escape(names(d)), collapse = " & "), " \\\\"), "\\midrule", rows,
    "\\bottomrule", "\\end{tabular}", "\\end{adjustbox}", paste0("\\caption{", tex_escape(caption), "}"),
    paste0("\\label{tab:", file, "}"), "\\end{table}"), file.path(OUT_DIR, paste0(file, ".tex")))
}

stack_labels <- c("mat:out" = "Maturity x Outcome", "AI_bin:out" = "AI x Outcome", "mat:AI_bin:out" = "Maturity x Outcome x AI")
single_labels <- c("mat" = "Maturity", "AI_bin" = "AI", "mat:AI_bin" = "Maturity x AI")

# TABLE A1
summary_vars <- c("ln_xi_real", "ln_cites_5y", "y0", "y1", "gap", "mat", "AI_bin", control_vars())
descriptives <- do.call(rbind, lapply(c("All", "AI", "Non-AI"), function(group) {
  dd <- if (group == "All") d5 else d5[d5$AI_bin == as.integer(group == "AI"), , drop = FALSE]
  do.call(rbind, lapply(summary_vars, function(v) data.frame(group = group, variable = v, n = nrow(dd),
    mean = mean(dd[[v]]), sd = sd(dd[[v]]), min = min(dd[[v]]), median = median(dd[[v]]), max = max(dd[[v]]))))
}))
export_data(descriptives, "tableA1_descriptives", "Descriptive statistics for the climate baseline estimation sample after complete-case selection. Gap denotes standardized citations minus standardized private value.")
write_csv(d5[, c("patent_id", "topic", "grant_year", "AI_bin", "mat", "mat_cumshare", "mat_peakrate", "y0", "y1", "gap")], "main_estimation_sample.csv")

# TABLE A2
fe_data <- d5
fe_models <- setNames(lapply(names(fe_specs), function(fe) fit_stack(fe_data, paste0("A2_", fe), fe = fe)), names(fe_specs))
export_models(fe_models, "tableA2_FE_progression", stack_labels,
  "Fixed-effect progression on a common sample. FE1: year-by-outcome; FE2: year-by-outcome and CPC-section-by-outcome; FE3: year-by-outcome and topic-by-outcome; FE5: patent, year-by-outcome and topic-by-outcome. FE1-FE3 fully interact regressors with outcome; FE5 includes all outcome-specific lower-order terms. All columns include ABCD controls.")
gap_check <- fit_single(fe_data, "gap", "equivalence_gap")
checks <- data.frame(model = c("Fully interacted FE3", "Patent FE5", "Direct difference"),
  estimate = c(extract(fe_models$FE3, "mat:AI_bin:out")$estimate, extract(fe_models$FE5, "mat:AI_bin:out")$estimate, extract(gap_check, "mat:AI_bin")$estimate))
write_csv(checks, "coefficient_equivalence.csv")
if (max(abs(checks$estimate - checks$estimate[3])) > 1e-6) stop("FE3, FE5 and gap estimates disagree; inspect samples and convergence.")
rm(fe_models, fe_data, gap_check)
gc(FALSE)

# TABLE A3
blocks <- c("A", "AB", "ABC", "ABCD")
control_models <- setNames(lapply(blocks, function(b) fit_stack(d5, paste0("A3_", b), block = b)), blocks)
main_model <- control_models$ABCD
stopifnot(main_model$patents == nrow(d5), main_model$n == 2L * nrow(d5))
export_models(control_models, "tableA3_control_progression", stack_labels,
  "Control progression on the same estimation sample after complete-case selection. All columns include patent, topic-by-outcome and year-by-outcome fixed effects. Controls are interacted with outcome. A: quality and scope; B: team and assignee; C: prior art; D: science linkage. Outcome is zero for value and one for citations.")
export_models(list("Baseline" = main_model), "table_main_stacked", stack_labels,
  "Baseline stacked model. Outcomes are separately winsorized log private value and log five-year forward citations, standardized on the baseline sample. Patent, topic-by-outcome and year-by-outcome fixed effects and ABCD controls are included.")
rm(control_models)

# TABLE A4
robust <- list("Validated 5y" = main_model)
run_variant <- function(citation, window, label, climate = TRUE, ids = NULL) {
  dd <- prepare(citation, window, climate, ids = ids, label = label)
  fit_stack(dd, label)
}
if (RUN_FULL_CPC) robust[["Full CPC 5y"]] <- run_variant("ln_cites_5y", 5L, "A4_full", FALSE)
robust[["Validated 3y"]] <- run_variant("ln_cites_3y", 3L, "A4_3y")
robust[["Examiner-added 5y"]] <- run_variant("ln_examiner_forward_cites_5y", 5L, "A4_examiner5")
robust[["Examiner-added 3y"]] <- run_variant("ln_examiner_forward_cites_3y", 3L, "A4_examiner3")
robust[["Through 2020"]] <- fit_stack(prepare("ln_cites_5y", 5L, max_year = 2020, label = "climate_through2020"), "A4_through2020")
robust[["Patent-clustered SE"]] <- fit_stack(d5, "A4_patent_cluster", cluster = "patent_id")
robust[["Topic and year clusters"]] <- fit_stack(d5, "A4_two_way", cluster = "topic + grant_year")
robust[["Continuous AI"]] <- fit_stack(d5, "A4_continuous", exposure = "ai_int_c")
matched_ids <- character()
if (RUN_MATCHING) {
  md <- prepare("ln_cites_5y", 5L, label = "matching_population")
  md <- md[, c("patent_id", "AI_bin", "grant_year", control_vars()), drop = FALSE]
  md$year_bin <- cut(md$grant_year, breaks = seq(min(md$grant_year, na.rm = TRUE),
    max(md$grant_year, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)
  md <- md[!is.na(md$AI_bin) & !is.na(md$year_bin), , drop = FALSE]
  mm <- MatchIt::matchit(AI_bin ~ ln_num_claims + ln_applicant_bwd + ln_examiner_bwd +
    female_inventor_share + ln_n_assignees + has_company_assignee + has_government_assignee,
    data = md, method = "nearest", distance = "glm", exact = ~year_bin, ratio = 1,
    replace = FALSE, estimand = "ATT")
  mb <- summary(mm, standardize = TRUE)
  write_csv(data.frame(variable = rownames(mb$sum.all), mb$sum.all, row.names = NULL), "matching_balance_before.csv")
  write_csv(data.frame(variable = rownames(mb$sum.matched), mb$sum.matched, row.names = NULL), "matching_balance_after.csv")
  matched <- MatchIt::match.data(mm, data = md)
  if (any(abs(matched$weights - 1) > 1e-10)) stop("Matching weights are not uniform; weighted estimation is required.")
  matched_ids <- matched$patent_id
  write_csv(matched[, c("patent_id", "AI_bin", "weights", "subclass")], "matched_patents.csv")
  match5 <- prepare("ln_cites_5y", 5L, ids = matched_ids, label = "matched5")
  match5$match_pair <- as.character(matched$subclass[match(match5$patent_id, matched$patent_id)])
  robust[["Matched 5y"]] <- fit_stack(match5, "A4_matched5", cluster = "topic + match_pair")
  match3 <- prepare("ln_cites_3y", 3L, ids = matched_ids, label = "matched3")
  match3$match_pair <- as.character(matched$subclass[match(match3$patent_id, matched$patent_id)])
  robust[["Matched 3y"]] <- fit_stack(match3, "A4_matched3", cluster = "topic + match_pair")
  rm(md, mm, mb, matched, match5, match3)
}
export_models(robust, "tableA4_robustness", c(stack_labels, "ai_int_c:out" = "AI intensity x Outcome", "mat:ai_int_c:out" = "Maturity x AI intensity x Outcome"),
  "Sensitivity analyses with patent, topic-by-outcome and year-by-outcome fixed effects and ABCD controls. Samples use the available nonmissing precomputed citation outcomes, with no additional calendar-year cutoff. Topic -1 is excluded. Required outcomes, controls and fixed-effect identifiers must be observed; undefined science confidence is zero only when has_science_ref equals zero. Full available-sample winsorization and sample-specific maturity and standardization retain their original order. Matched models cluster by topic and matched pair. Matching is descriptive balance adjustment, not causal identification. Continuous AI has different units from binary AI.")
rm(robust)
gc(FALSE)

# TABLE A5
threshold_vars <- paste0("AI_", c(30, 40, 50, 60, 70))
dt <- d5
threshold_models <- setNames(lapply(threshold_vars, function(v) {
  dd <- dt
  dd$AI_bin <- dd[[v]]
  fit_stack(dd, paste0("A5_", v))
}), threshold_vars)
export_models(threshold_models, "tableA5_thresholds", stack_labels,
  "AI thresholds using the original baseline preprocessing; each model excludes missing values of its own AI definition. Each model includes the relevant AI-by-outcome term, the maturity-by-outcome term and their triple interaction; patent and topic/year-by-outcome fixed effects and ABCD controls.")
rm(dt, threshold_models)

# TABLE A6
maturity_vars <- c("mat_logistic", "mat_cumshare", "mat_peakrate")
maturity_models <- setNames(lapply(maturity_vars, function(v) {
  dd <- d5
  dd$mat <- dd[[v]]
  fit_stack(dd, paste0("A6_", v))
}), c("Logistic (original algorithm)", "Cumulative share", "Peak-rate"))
export_models(maturity_models, "tableA6_maturity", stack_labels,
  "Alternative retrospective maturity measures on a common sample. Maturity uses the original baseline sample, observed topic-year counts and original fitting algorithm. Patent and topic/year-by-outcome fixed effects and ABCD controls are included. Coefficient magnitudes depend on the measure scale. The original algorithm falls back to cumulative share for insufficient observations or fitting errors; fitting warnings must be reviewed.")
rm(maturity_models)

# TABLE A7
single_models <- list("Value standardized" = fit_single(d5, "y0", "A7_value_std"),
  "Citations standardized" = fit_single(d5, "y1", "A7_cites_std"),
  "Value raw log" = fit_single(d5, "ln_xi_real", "A7_value_raw"),
  "Citations raw log" = fit_single(d5, "cites_raw", "A7_cites_raw"))
if (RUN_FULL_CPC) {
  dd <- prepare("ln_cites_5y", 5L, FALSE)
  single_models[["Value Full CPC raw log"]] <- fit_single(dd, "ln_xi_real", "A7_full_value")
  single_models[["Citations Full CPC raw log"]] <- fit_single(dd, "cites_raw", "A7_full_cites")
}
dd <- prepare("ln_cites_3y", 3L)
single_models[["Citations 3y raw log"]] <- fit_single(dd, "cites_raw", "A7_cites3")
if (length(matched_ids)) {
  dd <- prepare("ln_cites_5y", 5L, ids = matched_ids, label = "matched5")
  single_models[["Matched value raw log"]] <- fit_single(dd, "ln_xi_real", "A7_match_value")
  single_models[["Matched citations raw log"]] <- fit_single(dd, "cites_raw", "A7_match_cites")
}
decomposition <- extract(single_models[["Citations standardized"]], "mat:AI_bin")$estimate - extract(single_models[["Value standardized"]], "mat:AI_bin")$estimate
write_csv(data.frame(stacked = extract(main_model, "mat:AI_bin:out")$estimate, separate_difference = decomposition), "separate_outcome_equivalence.csv")
if (abs(decomposition - extract(main_model, "mat:AI_bin:out")$estimate) > 1e-6) stop("Separate standardized outcomes do not reproduce the stacked coefficient.")
export_models(single_models, "tableA7_separate_outcomes", single_labels,
  "Separate outcome regressions with topic and grant-year fixed effects and ABCD controls. Only the standardized baseline columns exactly decompose the stacked estimate. Raw-log columns use the supplied log outcomes without winsorization.")
rm(single_models, dd)

# TABLE A8
diffusion_vars <- c("ln_total_cites_edges", "ln_cites_0_3_edges", "ln_cites_3_5_edges", "early_share_3of5",
  "ln_external_cites", "ln_self_cites", "self_cite_share", "ln_citing_assignee_breadth", "ln_citing_topic_breadth",
  "ln_cross_topic_cites", "ln_within_topic_cites", "cross_topic_share")
for (v in c("early_share_3of5", "self_cite_share", "cross_topic_share")) {
  if (any(d5[[v]] < 0 | d5[[v]] > 1, na.rm = TRUE)) stop("Share outside [0,1]: ", v)
}
diffusion_models <- setNames(lapply(diffusion_vars, function(v) fit_single(d5, v, paste0("A8_", v))), diffusion_vars)
export_models(diffusion_models, "tableA8_citation_mechanisms", single_labels,
  "Associations with citation timing and diffusion, conditional on baseline-sample membership. Topic and grant-year fixed effects and ABCD controls. Outcome-specific missingness is retained; undefined shares are not set to zero. These regressions do not identify causal mechanisms or measure environmental benefits. Verify that the upstream edge measures use the stated citation horizon.")
diffusion_results <- do.call(rbind, lapply(names(diffusion_models), function(v) {
  z <- extract(diffusion_models[[v]], "mat:AI_bin")
  z$outcome <- v
  z$patents <- diffusion_models[[v]]$patents
  z
}))
diffusion_results$p_BH <- p.adjust(diffusion_results$p_value, "BH")
write_csv(diffusion_results, "tableA8_interactions_BH.csv")
rm(diffusion_models)

# TABLE A9
components <- c("ml", "evo", "nlp", "speech", "vision", "kr", "planning", "hardware")
score_vars <- paste0("ai_score_", components)
pred_vars <- paste0("predict50_", components)
component_sample <- d5
S <- as.matrix(component_sample[, score_vars, drop = FALSE])
storage.mode(S) <- "double"
ss <- rowSums(S)
wvars <- paste0("ai_w_", components)
component_sample[wvars] <- as.data.frame(S / ifelse(ss > 0, ss, 1) * component_sample$ai_intensity)
weighted <- fit_stack(component_sample, "A9_weighted", exposure = wvars)
weighted_labels <- setNames(c("Maturity x Outcome", paste(components, "x Outcome"), paste("Maturity x", components, "x Outcome")),
  c("mat:out", paste0(wvars, ":out"), paste0("mat:", wvars, ":out")))
export_models(list("Joint weighted components" = weighted), "tableA9_subtypes_weighted", weighted_labels,
  "Joint model of component scores normalized to shares and multiplied by overall AI intensity. All component-by-outcome and maturity-by-component-by-outcome terms are included, with patent and topic/year-by-outcome fixed effects and ABCD controls. These are score-based exposures, not calibrated probabilities or mutually exclusive technologies.")
wresults <- do.call(rbind, lapply(wvars, function(v) {z <- extract(weighted, paste0("mat:", v, ":out")); z$component <- v; z}))
wresults$p_BH <- p.adjust(wresults$p_value, "BH")
write_csv(wresults, "tableA9_interactions_BH.csv")
rm(S, ss, weighted)

# TABLE A10
indicator_models <- setNames(lapply(pred_vars, function(v) {
  dd <- component_sample
  dd$AI_bin <- dd[[v]]
  fit_stack(dd, paste0("A10_", v))
}), components)
export_models(indicator_models, "tableA10_subtypes_indicators", stack_labels,
  "Separate subtype-indicator models using baseline preprocessing and each model's complete observations. AI denotes the indicated subtype; its comparison group can include other AI subtypes. Each model includes all lower-order terms, patent and topic/year-by-outcome fixed effects and ABCD controls. Estimates are not independent subtype effects.")
indicator_results <- do.call(rbind, lapply(names(indicator_models), function(v) {z <- extract(indicator_models[[v]], "mat:AI_bin:out"); z$component <- v; z}))
indicator_results$p_BH <- p.adjust(indicator_results$p_value, "BH")
write_csv(indicator_results, "tableA10_interactions_BH.csv")
rm(component_sample, indicator_models)
gc(FALSE)

# TABLE A11
loo <- do.call(rbind, lapply(sort(unique(d5$topic)), function(t) {
  m <- fit_stack(d5[d5$topic != t, , drop = FALSE], paste0("A11_exclude_", t))
  z <- extract(m, "mat:AI_bin:out")
  z$excluded_topic <- t
  z$patents <- m$patents
  z
}))
write_csv(loo, "tableA11_LOO_all.csv")
export_data(data.frame(statistic = c("min", "mean", "max"), value = c(min(loo$estimate), mean(loo$estimate), max(loo$estimate))), "tableA11_LOO_summary", "Leave-one-topic-out estimates of the maturity-by-AI difference in the citation-value gradient. Maturity and outcome scaling retain their baseline definitions.")

# TABLE A12
firm_data <- d5
firm_data$primary_assignee <- trimws(sub("[;|,].*$", "", firm_data$primary_assignee_id))
firm_data <- firm_data[valid(firm_data, "primary_assignee"), , drop = FALSE]
portfolio <- as.data.frame(table(firm_data$primary_assignee), stringsAsFactors = FALSE)
names(portfolio) <- c("primary_assignee", "patents_in_analysis_sample")
firm_data$portfolio_n <- portfolio$patents_in_analysis_sample[match(firm_data$primary_assignee, portfolio$primary_assignee)]
cutoffs <- quantile(portfolio$patents_in_analysis_sample, c(.75, .90, .95), names = FALSE)
write_csv(portfolio, "assignee_sample_portfolios.csv")
lincomb <- function(m, terms, weights) {
  keys <- canon(terms)
  ix <- match(keys, m$coefs$term)
  if (anyNA(ix)) stop("Unidentified group contrast.")
  b <- sum(m$coefs$estimate[ix] * weights)
  variance <- as.numeric(t(weights) %*% m$V[keys, keys, drop = FALSE] %*% weights)
  if (!is.finite(variance) || variance <= 0) stop("Nonpositive variance for group contrast.")
  se <- sqrt(variance)
  data.frame(estimate = b, std_error = se, p_value = 2 * pt(-abs(b / se), m$df), conf_low = b - qt(.975, m$df) * se, conf_high = b + qt(.975, m$df) * se)
}
firm_models <- list()
firm_results <- list()
for (i in seq_along(cutoffs)) {
  label <- c("Top 25 percent", "Top 10 percent", "Top 5 percent")[i]
  fd <- firm_data
  fd$large <- as.integer(fd$portfolio_n >= cutoffs[i])
  dat <- stack(fd, c("mat", "AI_bin", "large", control_vars()))
  rhs <- paste(c("mat:out", "AI_bin:out", "large:out", "mat:AI_bin:out", "mat:large:out", "AI_bin:large:out", "mat:AI_bin:large:out", paste0(control_vars(), ":out")), collapse = " + ")
  m <- fit_model(paste("outcome_std ~", rhs, "|", fe_specs[["FE5"]]), dat, paste0("A12_", i),
    expected = c("mat:AI_bin:out", "mat:AI_bin:large:out"))
  firm_models[[label]] <- m
  small <- lincomb(m, "mat:AI_bin:out", 1)
  additional <- lincomb(m, "mat:AI_bin:large:out", 1)
  large <- lincomb(m, c("mat:AI_bin:out", "mat:AI_bin:large:out"), c(1, 1))
  z <- rbind(small, additional, large)
  z$contrast <- c("Smaller portfolio", "Additional larger-portfolio gradient", "Larger portfolio")
  z$threshold <- label
  z$cutoff_patents <- cutoffs[i]
  z$actual_large_assignee_share <- mean(portfolio$patents_in_analysis_sample >= cutoffs[i])
  z$patents <- nrow(fd)
  firm_results[[label]] <- z
}
export_models(firm_models, "tableA12_large_firms", NULL,
  "Assignee portfolio heterogeneity. Portfolio size counts patents in the climate analysis sample, not total firm size or global patent holdings. Missing primary assignees are excluded. All lower-order outcome interactions are included. Quantile thresholds are across assignees; ties can change the realized group shares.")
export_data(do.call(rbind, firm_results), "tableA12_group_contrasts", "Estimated gradients for smaller and larger assignee portfolios, and their difference. Standard errors and intervals for the larger-portfolio gradient use the full covariance of the summed coefficients.")
rm(firm_data, fd, dat, firm_models, firm_results)

fig <- ggplot2::ggplot(
  fig_summary,
  ggplot2::aes(
    x = maturity, y = mean,
    colour = AI, fill = AI, group = AI
  )
) +
  ggplot2::geom_hline(
    yintercept = 0, linetype = "dashed", colour = "grey60"
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = low, ymax = high),
    alpha = .13, colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_colour_manual(
    values = c(AI = "#C49A00", "Non-AI" = "#5B2A86")
  ) +
  ggplot2::scale_fill_manual(
    values = c(AI = "#C49A00", "Non-AI" = "#5B2A86")
  ) +
  ggplot2::labs(
    x = "Technology maturity",
    y = "Private value − technological influence (std.)"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(OUT_DIR, "figure2_wedge_by_maturity.pdf"),
  fig, width = 7, height = 4.5
)

print(fig)