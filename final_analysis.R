# Does AI Change Candidate Trust? The Role of Advertisement Source and Background Color
# Run in RStudio with Source, or: Rscript final_analysis.R [input_csv] [output_directory]
# Uses R 4.x and its standard packages. Default output directory: results.
#
# Professor's instructions: candidate trust, binary treatment coding, eight models
# in order, an eight-column table, and four condition predictions from the interaction.
# Project decisions: categorical controls, common sample, HC1 for figures, HC3
# sensitivity, M8 standardization, equal-weight contrasts, graphics and packaging.
#
# Input is the final CSV rebuilt from the response workbook by scripts/prepare_workbook.py.
# Team-confirmed tabs: 1 AI/blue, 2 AI/red, 3 non-AI/blue, 4 non-AI/red.
# The earlier CSV supplies a private comparison baseline and stable IDs only.

# 1. Validate the analysis file before making any model or plot ----------------
args <- if (interactive()) character() else commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript final_analysis.R [input_csv] [output_directory]")
input <- if (length(args) >= 1L) args[1] else file.path("data", "survey_data.csv")
if (!length(args) && !file.exists(input)) {
  if (file.exists("survey_data.csv")) {
    input <- "survey_data.csv"
  } else if (interactive()) {
    message("Select the final survey_data.csv in the file chooser.")
    input <- file.choose()
  }
}
out <- if (length(args) >= 2L) args[2] else "results"
if (!file.exists(input)) stop("Missing analysis CSV: ", input, ". See the README for local replication.")
raw <- read.csv(input, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                na.strings = c("", "NA"), check.names = FALSE)
demos <- c("age_range", "gender", "political_preference")
outcomes <- c("trust_candidate", "believe_policy", "want_more_ads")
required <- c("response_id", "ad_source", "background", "condition", demos, outcomes, "n_outcomes")
if (anyDuplicated(names(raw)) || !setequal(names(raw), required))
  stop("Unexpected columns. Expected the analysis schema: ", paste(required, collapse = ", "))
if (!nrow(raw)) stop("No survey responses found.")
if (anyNA(raw$response_id) || anyDuplicated(raw$response_id) || any(!nzchar(trimws(raw$response_id))))
  stop("Every response needs a unique, nonblank ID.")
if (!all(raw$ad_source %in% c("ai", "non_ai")) || !all(raw$background %in% c("blue", "red")))
  stop("Invalid source or background. Use ai/non_ai and blue/red.")
if (anyNA(raw$condition) || any(raw$condition != paste(raw$ad_source, raw$background, sep = "_")))
  stop("Condition disagrees with its source and background.")
for (v in outcomes) {
  x <- raw[[v]][!is.na(raw[[v]])]
  if (length(x) && (!is.numeric(x) || !all(x %in% 1:10)))
    stop(v, " contains a noninteger or out-of-range rating (allowed: 1 through 10).")
}
if (anyNA(raw$n_outcomes) || any(raw$n_outcomes != rowSums(!is.na(raw[outcomes]))))
  stop("n_outcomes disagrees with the three outcome fields.")

# These are the verified label spellings. Explicit refusal is a recorded category.
allowed <- list(age_range = c("18-24", "25-34", "45-54"),
                gender = c("Female", "Male", "Non-Binary"),
                political_preference = c("Democratic", "Conservative", "Moderate", "Prefer Not to Say"))
for (v in demos) {
  unknown <- setdiff(stats::na.omit(unique(raw[[v]])), allowed[[v]])
  if (length(unknown)) stop("Unrecognised ", v, " label(s): ", paste(unknown, collapse = ", "))
}

# Eligibility uses trust and all adjustment variables, even in unadjusted models.
# Policy belief and interest are separate outcomes and cannot exclude a trust score.
eligible <- complete.cases(raw[c("trust_candidate", demos)])
d <- raw[eligible, , drop = FALSE]
N <- nrow(d)
if (!N) stop("No complete cases for the eight trust models.")
factor_reference <- function(x, reference) {
  if (!reference %in% x) stop("Reference category absent from the common sample: ", reference)
  if (length(unique(x)) < 2) stop("Categorical control has only one observed level: ", reference)
  relevel(factor(x), ref = reference)
}
for (v in demos) d[[v]] <- factor_reference(d[[v]], allowed[[v]][1L])
options(contrasts = c("contr.treatment", "contr.poly"))

# Professor's binary coding. Red = 0 and AI = 0 identifies blue/non-AI.
d$red <- as.integer(d$background == "red")
d$ai <- as.integer(d$ad_source == "ai")
d$red_ai <- d$red * d$ai
d$cell <- 1L + d$ai + 2L * d$red
cell_codes <- c("non_ai_blue", "ai_blue", "non_ai_red", "ai_red")
counts <- tabulate(d$cell, nbins = 4L)
if (any(counts == 0)) stop("An empty analysis cell prevents identification of the interaction.")

tab <- file.path(out, "tables")
poster <- file.path(out, "poster")
supp <- file.path(out, "supplementary")
for (path in c(tab, poster, supp)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create output directory: ", path)
}
write_table <- function(x, name) write.csv(x, file.path(tab, name), row.names = FALSE, na = "NA")

# 2. Cell sizes and descriptive summaries ------------------------------------
total_counts <- as.integer(table(factor(raw$condition, levels = cell_codes)))
write_table(data.frame(condition = cell_codes, submissions = total_counts, model_n = counts,
                       excluded = total_counts-counts), "condition_counts.csv")
missing_fields <- data.frame(field = required, n_missing = colSums(is.na(raw[required])))
write_table(missing_fields, "missingness.csv")

# Compute descriptives using available answers to each question, independently of
# complete-case selection for the regressions. Separately export exploratory outcomes.
descriptive <- do.call(rbind, lapply(outcomes, function(y) {
  do.call(rbind, lapply(cell_codes, function(g) {
    z <- raw[raw$condition == g, y]
    x <- z[!is.na(z)]
    n <- length(x)
    data.frame(outcome = y, condition = g, submissions = length(z), n = n,
               mean = if (n) mean(x) else NA_real_, sd = sd(x), median = median(x),
               min = if (n) min(x) else NA_real_, max = if (n) max(x) else NA_real_)
  }))
}))
write_table(subset(descriptive, outcome == "trust_candidate"), "trust_descriptives.csv")
write.csv(subset(descriptive, outcome != "trust_candidate"),
          file.path(supp, "exploratory_outcome_descriptives.csv"), row.names = FALSE, na = "NA")

# 3. The professor's eight specifications, in the original sequence -----------
rhs <- c(
  M1 = "red",
  M2 = "red + age_range + gender + political_preference",
  M3 = "ai",
  M4 = "ai + age_range + gender + political_preference",
  M5 = "red + ai",
  M6 = "red + ai + red_ai",
  M7 = "red + ai + age_range + gender + political_preference",
  M8 = "red + ai + red_ai + age_range + gender + political_preference"
)
models <- lapply(rhs, function(x) lm(as.formula(paste("trust_candidate ~", x)),
                                     data = d, na.action = na.fail))
for (id in names(models)) {
  m <- models[[id]]
  if (anyNA(coef(m)) || m$rank != ncol(model.matrix(m)))
    stop(id, " is rank deficient. Review the design matrix before inference.")
  if (df.residual(m) <= 0) stop(id, " has no residual degrees of freedom.")
  if (!identical(rownames(model.frame(m)), rownames(d))) stop(id, " used a different sample.")
}

# OLS estimates are identical across covariance choices. HC1 uses n/(n-rank).
# HC3 divides residuals by (1-h); a unit-leverage row makes the whole estimator
# undefined. Never replace such corrections or influence diagnostics with zero.
is_unit <- function(h) 1 - h < 1e-8
robust_vcov <- function(m, type) {
  X <- model.matrix(m)
  B <- solve(crossprod(X))
  u <- residuals(m)
  if (type == "HC3") {
    h <- hatvalues(m)
    if (any(is_unit(h))) return(NULL)
    u <- u / (1 - h)
  }
  Xu <- X * as.numeric(u)
  V <- B %*% crossprod(Xu) %*% B
  if (type == "HC1") V <- V * nobs(m) / df.residual(m)
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}
covariances <- list(OLS = lapply(models, vcov), HC1 = lapply(models, robust_vcov, type = "HC1"),
                    HC3 = lapply(models, robust_vcov, type = "HC3"))
fit_stats <- do.call(rbind, lapply(names(models), function(id) {
  m <- models[[id]]
  data.frame(model = id, specification = unname(rhs[id]), n = nobs(m), rank = m$rank,
             n_columns = ncol(model.matrix(m)), residual_df = df.residual(m),
             R2 = summary(m)$r.squared, adjusted_R2 = summary(m)$adj.r.squared,
             max_leverage = max(hatvalues(m)), unit_leverage_n = sum(is_unit(hatvalues(m))),
             OLS_available = TRUE, HC1_available = TRUE, HC3_available = !is.null(covariances$HC3[[id]]))
}))
write_table(fit_stats, "model_fit.csv")

# Use the covariance of the complete linear combination including all covariances.
# Inference uses two-sided t tests and pointwise intervals with model residual df.
inference <- function(m, V, L) {
  b <- as.vector(L %*% coef(m))
  variances <- diag(L %*% V %*% t(L))
  if (any(!is.finite(variances)) || any(variances < -1e-9)) stop("Invalid contrast covariance.")
  se <- sqrt(pmax(variances, 0))  # guard only against round-off on a zero variance
  df <- df.residual(m)
  q <- qt(.975, df)
  p <- rep(NA_real_, length(se))
  positive <- se > 0
  p[positive] <- 2 * pt(-abs(b[positive] / se[positive]), df)
  data.frame(estimate = b, se = unname(se), lo = b - q*se, hi = b + q*se, p = p)
}
coefficients_out <- do.call(rbind, lapply(names(covariances), function(method) {
  do.call(rbind, lapply(names(models), function(id) {
    m <- models[[id]]
    V <- covariances[[method]][[id]]
    if (is.null(V)) return(NULL)
    data.frame(model = id, term = names(coef(m)), method = method,
               inference(m, V, diag(m$rank)), n = nobs(m), residual_df = df.residual(m))
  }))
}))
write_table(coefficients_out, "coefficients_all_methods.csv")

# 4. Professor's eight-column table: one model per column, omitted terms blank
pretty_term <- function(x) {
  x[x == "red"] <- "Red vs blue"
  x[x == "ai"] <- "AI vs non-AI"
  x[x == "red_ai"] <- "Red × AI"
  x[x == "(Intercept)"] <- "Intercept"
  x <- sub("^age_range", "Age range · ", x)
  x <- sub("^gender", "Gender · ", x)
  sub("^political_preference", "Political preference · ", x)
}
term_order <- unique(c("(Intercept)", "red", "ai", "red_ai", coefficients_out$term))
make_table <- function(estimator) {
  result <- data.frame(Term = c(pretty_term(term_order), "Observations", "Residual df", "R²", "Adjusted R²"))
  for (id in names(models)) {
    m <- models[[id]]
    b <- coef(m)
    included <- match(names(b), term_order)
    V <- covariances[[estimator]][[id]]
    cells <- rep("", length(term_order))
    if (is.null(V)) {
      # Coefficients are estimated even when ordinary HC3 uncertainty is undefined.
      cells[included] <- sprintf("%.3f (NA)", b)
    } else {
      se <- sqrt(pmax(diag(V), 0))
      cells[included] <- sprintf("%.3f (%.3f)", b, se)
    }
    f <- fit_stats[fit_stats$model == id, ]
    result[[id]] <- c(cells, f$n, f$residual_df, sprintf("%.3f", f$R2), sprintf("%.3f", f$adjusted_R2))
  }
  result
}
escape_html <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE)
html_table <- function(t) {
  header <- paste0("<th>", names(t), "</th>", collapse = "")
  body <- apply(t, 1, function(r) paste0("<tr>", paste0("<td>", escape_html(r), "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(body, collapse = "\n"), "</tbody></table>")
}
html <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Candidate trust regressions</title>',
          '<style>body{font:16px/1.5 Arial,sans-serif;color:#193245;margin:3%}table{border-collapse:collapse;font-size:13px}td,th{padding:7px 9px;border-bottom:1px solid #ddd;text-align:right}td:first-child,th:first-child{text-align:left}th{background:#eaf1f3}p{max-width:78em}</style>',
          '<h1>Candidate trust across eight model specifications</h1>',
          '<p>Cells report estimates in trust points with standard errors in parentheses. Empty cells indicate an omitted term. The machine-readable coefficient table includes confidence intervals and two-sided p-values.</p>',
          '<p>References are blue, non-AI, age 18–24, Female and Democratic. In M6 and M8 the AI coefficient compares sources under blue, and the red coefficient compares colors under non-AI.</p>')
for (method in names(covariances)) {
  wide <- make_table(method)
  write_table(wide, paste0("eight_model_table_", method, ".csv"))
  html <- c(html, paste0("<h2>", method, "</h2>"), html_table(wide))
}
html <- c(html, "<p>A coefficient followed by (NA) has undefined HC3 uncertainty because the model contains unit leverage. All eight models use the same respondents.</p></html>")
writeLines(enc2utf8(html), file.path(tab, "eight_model_tables.html"), useBytes = TRUE)

# 5. Professor's four predictions; project extension to standardized M8 means
# M6 cell order: b0; b0 + bAI; b0 + bRed; b0 + bRed + bAI + bInteraction.
# In M6/M8, bAI compares sources under blue and bRed compares colors under non-AI.
# bInteraction is the red AI gap minus the blue AI gap.
prediction_design <- function(m) {
  do.call(rbind, lapply(1:4, function(g) {
    # Clone the entire common sample for each condition; keep all demographics
    # exactly as observed. Averaging design rows also gives the SE of the mean.
    new <- d
    new$red <- as.integer(g > 2)
    new$ai <- as.integer(g %in% c(2, 4))
    new$red_ai <- new$red * new$ai
    colMeans(model.matrix(delete.response(terms(m)), new, contrasts.arg = m$contrasts, xlev = m$xlevels))
  }))
}
inter_models <- c("M6", "M8")
pred <- do.call(rbind, lapply(inter_models, function(id) {
  m <- models[[id]]
  stats <- inference(m, covariances$HC1[[id]], prediction_design(m))
  data.frame(model = id, condition = cell_codes, observed_n = counts, method = "HC1",
             stats[c("estimate", "se", "lo", "hi")])
}))
write_table(pred, "condition_predictions_HC1.csv")

# Columns below are non-AI/blue, AI/blue, non-AI/red, AI/red, respectively.
# Contrasts with "equal" in the name give each level of the other factor 50%.
weights <- rbind(AI_equal_colors = c(-.5,.5,-.5,.5), AI_blue = c(-1,1,0,0),
                 AI_red = c(0,0,-1,1), interaction = c(1,-1,-1,1),
                 red_equal_sources = c(-.5,-.5,.5,.5))
contrasts_out <- do.call(rbind, lapply(names(covariances), function(method) {
  do.call(rbind, lapply(inter_models, function(id) {
    m <- models[[id]]
    V <- covariances[[method]][[id]]
    if (is.null(V)) return(NULL)
    L <- weights %*% prediction_design(m)
    data.frame(model = id, contrast = rownames(weights), method = method, inference(m, V, L))
  }))
}))
write_table(contrasts_out, "direct_contrasts.csv")

# Cross-check estimands by a second route using treatment coefficients directly.
close_enough <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-9))
cell_means <- as.numeric(tapply(d$trust_candidate, d$cell, mean))
b6 <- coef(models$M6)
formula_means <- unname(b6["(Intercept)"]) + c(0, b6["ai"], b6["red"], sum(b6[c("red", "ai", "red_ai")]))
m6_means <- pred$estimate[pred$model == "M6"]
if (!close_enough(m6_means, formula_means) || !close_enough(m6_means, cell_means))
  stop("M6 predictions disagree with the four coefficient formulas or the observed cell means.")
if (!close_enough(hatvalues(models$M6), 1/counts[d$cell])) stop("M6 leverage check failed.")
for (id in inter_models) {
  b <- coef(models[[id]])
  direct <- c(b["ai"] + .5*b["red_ai"], b["ai"], b["ai"] + b["red_ai"],
              b["red_ai"], b["red"] + .5*b["red_ai"])
  rows <- contrasts_out$model == id & contrasts_out$method == "HC1"
  if (!close_enough(contrasts_out$estimate[rows], direct)) stop(id, " contrast identity failed.")
}

# 6. Diagnostics and sensitivity, without participant-level public labels -----
diagnostics <- lapply(models[inter_models], function(m) {
  h <- hatvalues(m)
  good <- !is_unit(h)
  std <- cook <- rep(NA_real_, N)
  if (summary(m)$sigma > 0) {
    std[good] <- residuals(m)[good] / (summary(m)$sigma * sqrt(1 - h[good]))
    cook[good] <- std[good]^2 * h[good] / (m$rank * (1 - h[good]))
  }
  data.frame(fitted = fitted(m), residual = residuals(m), leverage = h, std = std, cook = cook)
})
finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
diag_summary <- do.call(rbind, lapply(inter_models, function(id) {
  z <- diagnostics[[id]]
  data.frame(model = id, n = N, residual_sd = summary(models[[id]])$sigma,
             max_leverage = max(z$leverage), unit_leverage_n = sum(is_unit(z$leverage)),
             undefined_standardized_residuals = sum(is.na(z$std)), undefined_Cooks_D = sum(is.na(z$cook)),
             largest_defined_abs_standardized_residual = finite_max(abs(z$std)), largest_defined_Cooks_D = finite_max(z$cook))
}))
write_table(diag_summary, "diagnostic_summary.csv")

# Leave-one-out comparison for saturated M6. Replace the affected cell mean by
# its mean after removing one response. If a cell empties, its contrast is undefined.
loo <- vapply(seq_len(N), function(i) {
  g <- d$cell[i]
  if (counts[g] < 2) return(NA_real_)
  means <- cell_means
  means[g] <- (counts[g] * means[g] - d$trust_candidate[i]) / (counts[g] - 1)
  sum(weights["AI_equal_colors", ] * means)
}, numeric(1))
loo_range <- if (all(is.na(loo))) c(NA_real_, NA_real_) else range(loo, na.rm = TRUE)
write_table(data.frame(model = "M6", contrast = "AI_equal_colors", defined_deletions = sum(!is.na(loo)),
                       undefined_deletions = sum(is.na(loo)), min_estimate = loo_range[1], max_estimate = loo_range[2]),
            "leave_one_out_summary.csv")

# Keep the original script's two secondary source comparisons in their own file.
# Use available scores for each secondary question and equal weights for colors.
# Both tests form one exploratory family. Holm adjusts p-values across these two
# planned outputs; the confidence intervals remain pointwise and unadjusted.
secondary <- do.call(rbind, lapply(outcomes[-1], function(y) {
  z <- raw[!is.na(raw[[y]]), ]
  z$red <- as.integer(z$background == "red")
  z$ai <- as.integer(z$ad_source == "ai")
  z$red_ai <- z$red * z$ai
  estimable <- nrow(z) > 4 && length(unique(z$condition)) == 4
  if (estimable) {
    m <- lm(reformulate(c("red", "ai", "red_ai"), response = y), data = z)
    L <- matrix(c(0, 0, 1, .5), nrow = 1)
    stats <- inference(m, robust_vcov(m, "HC1"), L)
  } else {
    stats <- data.frame(estimate = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
  }
  data.frame(outcome = y, contrast = "AI_equal_colors", n = nrow(z), method = "HC1",
             estimable = estimable, stats)
}))
secondary$p_holm_two_tests <- p.adjust(secondary$p, method = "holm", n = 2)
write.csv(secondary, file.path(supp, "exploratory_source_contrasts.csv"), row.names = FALSE, na = "NA")

# 7. Re-read exported tables as the single source of plotted estimates --------
pred_plot <- read.csv(file.path(tab, "condition_predictions_HC1.csv"))
contrast_plot <- subset(read.csv(file.path(tab, "direct_contrasts.csv")), method == "HC1")
coef_plot <- subset(read.csv(file.path(tab, "coefficients_all_methods.csv")), method == "HC1")
for (field in c("estimate", "se", "lo", "hi")) {
  if (!close_enough(pred_plot[[field]], pred[[field]]) ||
      !close_enough(contrast_plot[[field]], contrasts_out[[field]][contrasts_out$method == "HC1"]) ||
      !close_enough(coef_plot[[field]], coefficients_out[[field]][coefficients_out$method == "HC1"]))
    stop("Exported and in-memory values disagree for ", field)
}

# Normalized page coordinates reserve independent space for titles, legends,
# axes and notes on both devices. PDF text and marks stay vector-based.
ink <- "#162F42"
muted <- "#475F70"
gridcol <- "#DDE5E9"
blue <- "#146695"
redcol <- "#AB4841"
model_col <- c(M6 = "#384D61", M8 = "#147769")
txt <- function(x, y, label, size = 12, color = ink, face = "plain", just = "left", rot = 0) {
  if (names(grDevices::dev.cur()) == "pdf") {
    label <- gsub("[–−]", "-", as.character(label))
    label <- gsub("’", "'", label)
  }
  grid::grid.text(label, x = x, y = y, just = just, rot = rot,
                  gp = grid::gpar(fontfamily = "sans", fontface = face, fontsize = size, col = color))
}
seg <- function(x1, y1, x2, y2, color = gridcol, width = 1, type = 1) {
  grid::grid.segments(x1, y1, x2, y2, gp = grid::gpar(col = color, lwd = width, lty = type, lineend = "butt"))
}
mark <- function(x, y, pch, color, size = 3.6) {
  grid::grid.points(x, y, default.units = "npc", pch = pch, size = grid::unit(size, "mm"),
                    gp = grid::gpar(col = color, fill = color))
}
title_line <- function(title, subtitle) {
  txt(.045, .933, title, 25, face = "bold")
  txt(.045, .865, subtitle, 13, muted)
}
signed <- function(x) ifelse(x < 0, sprintf("−%.2f", abs(x)), sprintf("+%.2f", x))
scale_to <- function(x, domain, lo, hi) lo + (x - domain[1])/diff(domain) * (hi-lo)
interval_axis <- function(x) {
  bounds <- range(c(0, x), finite = TRUE)
  ticks <- pretty(bounds + c(-1,1)*max(diff(bounds)*.025, .01), n = 5)
  list(domain = range(c(bounds,ticks)), ticks = ticks)
}
save_plot <- function(folder, stem, draw, width, height) {
  for (ext in c("pdf", "png")) {
    path <- file.path(folder, paste0(stem, ".", ext))
    if (ext == "pdf") {
      if (capabilities("cairo")) grDevices::cairo_pdf(path, width = width, height = height, family = "sans")
      else grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE,
                          encoding = "ISOLatin1.enc")
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = 350, bg = "white")
    }
    tryCatch({ grid::grid.newpage(); draw() }, finally = grDevices::dev.off())
  }
}

# Main figure 1. Both models share the full response scale. The limits also
# expand for any calculated intervals beyond 1–10. Horizontal offsets separate
# the background-specific whiskers at the same source category.
figure_predictions <- function() {
  title_line("Candidate trust across four advertisement conditions",
             "Predicted mean trust on a 1–10 scale. Whiskers show pointwise 95% HC1 confidence intervals.")
  legend_y <- .782
  for (bg in 0:1) {
    x <- c(.045, .295)[bg+1]
    color <- c(blue, redcol)[bg+1]
    seg(x, legend_y, x+.044, legend_y, color, 2.2, bg+1)
    mark(x+.022, legend_y, c(16,17)[bg+1], color)
    txt(x+.055, legend_y, c("Blue background", "Red background")[bg+1], 13)
  }
  bounds <- range(c(1, 10, pred_plot$lo, pred_plot$hi))
  ytick <- sort(unique(c(floor(bounds[1]), 1:10, ceiling(bounds[2]))))
  ys <- function(v) scale_to(v, range(ytick), .196, .63)
  for (j in seq_along(inter_models)) {
    id <- inter_models[j]
    p <- pred_plot[pred_plot$model == id, ]
    p <- p[match(cell_codes, p$condition), ]
    left <- c(.115, .565)[j]
    right <- left + .365
    txt((left+right)/2, .703, c("M6  ·  Without controls", "M8  ·  With controls")[j], 16, face = "bold", just = "centre")
    for (tick in ytick) {
      seg(left, ys(tick), right, ys(tick), if (tick %in% c(1,10)) "#A7B7C0" else gridcol)
      if (j == 1L) txt(left-.015, ys(tick), tick, 11, muted, just = "right")
    }
    centers <- left + c(.085, .280)
    for (bg in 0:1) {
      g <- 1:2 + 2*bg
      x <- centers + c(-.014,.014)[bg+1]
      color <- c(blue, redcol)[bg+1]
      seg(x[1], ys(p$estimate[g[1]]), x[2], ys(p$estimate[g[2]]), color, 2, bg+1)
      seg(x, ys(p$lo[g]), x, ys(p$hi[g]), color, 1.8)
      for (v in c("lo","hi")) seg(x-.006, ys(p[[v]][g]), x+.006, ys(p[[v]][g]), color, 1.8)
      mark(x, ys(p$estimate[g]), c(16,17)[bg+1], color, 4.2)
    }
    txt(centers, .149, c("Non-AI","AI"), 14, just = "centre")
  }
  txt(.044, .412, "Candidate trust (points)", 13, rot = 90, just = "centre")
  seg(.045, .107, .955, .107)
  txt(.045, .058, "M8 averages over the same respondents in all four conditions. Lines connect source categories.", 12, muted)
}

# Main figure 2. The printed estimates and plotted positions use identical rows.
# The symmetric pair of marks around each label separates M6 from M8.
figure_contrasts <- function() {
  title_line("Direct comparisons of candidate trust",
             "Differences in trust points. Whiskers show pointwise 95% HC1 confidence intervals.")
  legend_y <- .787
  for (j in 1:2) {
    id <- inter_models[j]
    x <- c(.046,.273)[j]
    seg(x, legend_y, x+.032, legend_y, model_col[id], 2)
    mark(x+.016, legend_y, c(16,18)[j], model_col[id], 4)
    txt(x+.044, legend_y, c("M6  Without controls", "M8  With controls")[j], 12)
  }
  left <- .398; right <- .775
  a <- interval_axis(c(contrast_plot$lo, contrast_plot$hi))
  xs <- function(v) scale_to(v, a$domain, left, right)
  y <- seq(.670, .260, length.out = 5)
  for (tick in a$ticks) {
    seg(xs(tick), .214, xs(tick), .719, if (tick == 0) "#6D7B87" else gridcol,
        if (tick == 0) 1.5 else .8, if (tick == 0) 2 else 1)
    txt(xs(tick), .176, as.character(tick), 11, muted, just = "centre")
  }
  hdr_x <- c(.839, .933)
  txt(hdr_x, .742, inter_models, 12, color = muted, face = "bold", just = "centre")
  labels <- c("AI minus non-AI", "AI minus non-AI", "AI minus non-AI",
              "Source-by-color interaction", "Red minus blue")
  context <- c("Equally averaged across colors", "Blue background", "Red background",
               "Red source gap minus blue source gap", "Equally averaged across sources")
  for (i in seq_along(y)) {
    txt(.046, y[i]+.014, labels[i], 14, face = "bold")
    txt(.046, y[i]-.025, context[i], 11, muted)
  }
  for (j in 1:2) {
    id <- inter_models[j]
    z <- contrast_plot[contrast_plot$model == id, ]
    z <- z[match(rownames(weights), z$contrast), ]
    yy <- y + c(.013,-.013)[j]
    seg(xs(z$lo), yy, xs(z$hi), yy, model_col[id], 2.0)
    for (v in c("lo","hi")) seg(xs(z[[v]]), yy-.006, xs(z[[v]]), yy+.006, model_col[id], 1.45)
    mark(xs(z$estimate), yy, c(16,18)[j], model_col[id], 4)
    txt(hdr_x[j], y, signed(z$estimate), 13, model_col[id], "bold", just = "centre")
  }
  txt((left+right)/2, .125, "Difference in trust points", 12, just = "centre")
  seg(.045, .091, .955, .091)
  txt(.045, .047, "Equal averages give each color or source 50% weight. Zero indicates equal predicted trust.", 12, muted)
}
save_plot(poster, "01_condition_trust", figure_predictions, 10.8, 6.1)
save_plot(poster, "02_direct_comparisons", figure_contrasts, 11.2, 6.5)

# The eight-model plot has one panel per coefficient. Model rows with a term
# excluded remain empty. Its axes share a scale to avoid visual exaggeration.
figure_coefficients <- function() {
  title_line("Treatment coefficients across eight models",
             "OLS estimates in trust points with pointwise 95% HC1 confidence intervals.")
  term_set <- c("red", "ai", "red_ai")
  plot_data <- coef_plot[coef_plot$term %in% term_set, ]
  a <- interval_axis(c(plot_data$lo, plot_data$hi))
  yy <- seq(.695, .255, length.out = 8)
  txt(.047, yy, names(models), 13, face = "bold")
  for (j in 1:3) {
    left <- c(.120, .418, .716)[j]
    right <- left+.238
    xs <- function(v) scale_to(v, a$domain, left, right)
    txt((left+right)/2, .777, c("Red minus blue","AI minus non-AI","Red × AI")[j], 15, face = "bold", just = "centre")
    seg(left, yy, right, yy)
    seg(xs(0), .222, xs(0), .728, muted, 1.2, 2)
    z <- plot_data[plot_data$term == term_set[j], ]
    y <- yy[match(z$model, names(models))]
    seg(xs(z$lo), y, xs(z$hi), y, model_col["M6"], 2)
    mark(xs(z$estimate), y, 16, model_col["M6"])
    txt(xs(a$ticks), .182, a$ticks, 11, muted, just = "centre")
  }
  seg(.045, .142, .955, .142)
  txt(.045, .094, "M1–M4 pool over the omitted factor. M5 and M7 impose a common difference across levels of the other factor.", 11.5, muted)
  txt(.045, .047, "M6 and M8 compare AI under blue and red under non-AI. Blank entries indicate that the term was omitted.", 11.5, muted)
}
save_plot(supp, "03_eight_model_coefficients", figure_coefficients, 11.2, 6.4)

# Fitted-versus-residual and leverage distributions. The public figure has
# no respondent labels, demographic coloring, row-level export or ID crosswalk.
figure_diagnostics <- function() {
  title_line("Residuals and leverage in the trust models",
             "Raw residuals show model fit. Leverage measures how unusual a respondent’s predictors are.")
  # One shared scale per column. Pad to keep symbols away from the frame.
  all_fit <- unlist(lapply(diagnostics, "[[", "fitted"))
  all_res <- unlist(lapply(diagnostics, "[[", "residual"))
  fx <- pretty(range(all_fit), n = 5)
  ry <- pretty(range(c(0,all_res)), n = 4)
  cuts <- c(0,.1,.2,.4,.6,.8,1+1e-8)
  leverage_bins <- lapply(diagnostics, function(z) tabulate(findInterval(z$leverage, cuts), nbins = 6))
  max_count <- max(unlist(leverage_bins))
  tops <- c(.72, .407); bottoms <- c(.505,.192)
  txt(.317, .79, "Raw residual vs fitted trust", 15, face = "bold", just = "centre")
  txt(.767, .79, "Responses by leverage band", 15, face = "bold", just = "centre")
  for (j in 1:2) {
    id <- inter_models[j]
    z <- diagnostics[[id]]
    y0 <- bottoms[j]; y1 <- tops[j]
    xs <- function(x) scale_to(x, range(fx), .148, .481)
    ys <- function(x) scale_to(x, range(ry), y0, y1)
    txt(.046, (y0+y1)/2, id, 16, face = "bold")
    for (v in ry) {
      seg(.148, ys(v), .481, ys(v), if (v == 0) muted else gridcol, 1, if (v == 0) 2 else 1)
      txt(.136, ys(v), v, 10, muted, just = "right")
    }
    txt(xs(fx), y0-.028, fx, 10, muted, just = "centre")
    mark(xs(z$fitted), ys(z$residual), 1, model_col[id], 2.9)
    bins <- leverage_bins[[id]]
    cx <- seq(.601,.922,length.out = 6)
    ybar <- y0 + (y1-y0)*.84*bins/max_count
    seg(.565,y0,.954,y0)
    for (k in 1:6) {
      grid::grid.rect(cx[k], y0, width=.037, height=ybar[k]-y0, just=c("centre","bottom"),
                      gp=grid::gpar(fill=model_col[id],col=NA))
      txt(cx[k], ybar[k]+.019, bins[k], 11, just = "centre")
    }
    txt(cx, y0-.030, c("0–.1",".1–.2",".2–.4",".4–.6",".6–.8",".8–1"), 10, muted, just = "centre")
  }
  seg(.045,.114,.955,.114)
  txt(.045,.080, "M6 has four fitted condition means. Points at identical coordinates overlap. Leverage bands include their lower endpoint.", 11, muted)
  txt(.045,.039, sprintf("Unit-leverage observations: M6 = %d, M8 = %d. Their standardized residuals and Cook’s distances are undefined.",
                         diag_summary$unit_leverage_n[1], diag_summary$unit_leverage_n[2]), 11, muted)
}
save_plot(supp, "04_model_diagnostics", figure_diagnostics, 11.2, 7.1)

# 8. Validation and reproducibility record -----------------------------------
# Merge workbook reconciliation only when its recorded checksum identifies the
# exact analysis CSV. A different input never inherits the release's audit claims.
audit_file <- file.path(dirname(input), "workbook_audit.txt")
digest <- unname(tools::md5sum(input))
audit <- if (file.exists(audit_file)) readLines(audit_file, encoding = "UTF-8") else character(0)
if (!paste("Analysis CSV MD5:", digest) %in% audit)
  audit <- c("WORKBOOK RECONCILIATION", "No matching workbook reconciliation supplied for this CSV.")
result_line <- function(id, k, method = "HC1") {
  z <- contrasts_out[contrasts_out$model == id & contrasts_out$contrast == k & contrasts_out$method == method, ]
  if (!nrow(z)) return(paste(id, k, method, "unavailable"))
  sprintf("%s %s (%s): %+.3f [%.3f, %.3f], p = %.4f.", id, k, method, z$estimate, z$lo, z$hi, z$p)
}
report <- c(audit, "", "ANALYSIS VALIDATION",
            sprintf("All submissions: %d. Common trust sample: %d. Excluded: %d.", nrow(raw), N, nrow(raw)-N),
            paste("Missing values by field:", paste(paste0(missing_fields$field,"=",missing_fields$n_missing),collapse=", ")),
            "All recorded ratings are numeric integers in [1,10]; condition strings and recorded outcome counts agree.",
            paste("Reference levels:", paste(sapply(demos, function(k) levels(d[[k]])[1]), collapse = ", ")),
            sapply(demos, function(k) paste(k, paste(names(table(d[[k]])), as.integer(table(d[[k]])), sep="=", collapse=", "))),
            sprintf("All eight fits are full rank and use identical observations. Range of residual df: %d to %d.",
                    min(fit_stats$residual_df), max(fit_stats$residual_df)),
            paste("Unit leverage per model:", paste(fit_stats$model, fit_stats$unit_leverage_n, sep="=", collapse=", ")),
            paste("Ordinary HC3 available:", paste(fit_stats$model[fit_stats$HC3_available], collapse=", ")),
            "Unit-leverage responses are retained; their standardized residuals and Cook's distances are undefined.",
            "HC1 avoids division by 1-h. Sparse categories still limit adjusted inference.",
            "M6 fitted means match the four professor-specified formulas and observed cell means. Leverages equal inverse cell sizes.",
            "M6 and M8 contrasts agree with direct coefficient formulas. M8 uses one empirical demographic distribution.",
            "Plotted estimates and interval endpoints were reread from the exported tables and verified.",
            "Intervals use residual-df t critical values; p-values are two-sided and pointwise comparisons are unadjusted.",
            "", "SENSITIVITY AND INTERPRETATION",
            unlist(lapply(inter_models, function(id) sapply(rownames(weights), function(k) result_line(id, k)))),
            "Equal-color AI contrasts under conventional OLS covariance and under HC3, where defined:",
            unlist(lapply(inter_models, function(id) sapply(c("OLS","HC3"), function(v) result_line(id,"AI_equal_colors",v)))),
            sprintf("M6 leave-one-out equal-color AI difference: %.3f to %.3f across %d defined deletions.", loo_range[1], loo_range[2], sum(!is.na(loo))),
            "The eight specifications reuse one dataset. Intervals spanning zero do not establish equivalence.",
            "Random assignment was proposed; implementation, independent participation and timing of controls were not verified.")
writeLines(enc2utf8(report), file.path(out, "validation_report.txt"), useBytes = TRUE)

# Retain platform, R and package versions. Replace machine-specific library paths
# occasionally printed for BLAS/LAPACK with basenames in the shareable record.
session <- capture.output(sessionInfo())
session <- sub("^(BLAS|LAPACK):[[:space:]]*([^ ]*/)", "\\1: ", session)
writeLines(session, file.path(out, "session_info.txt"))
cat(sprintf("Completed %d responses, %d common-sample observations and eight trust regressions.\n", nrow(raw), N))
print(data.frame(condition = cell_codes, n = counts), row.names = FALSE)
print(contrast_plot[c("model", "contrast", "estimate", "lo", "hi", "p")], row.names = FALSE, digits = 4)
cat("Eight-column tables:", file.path(tab, "eight_model_tables.html"), "\n")