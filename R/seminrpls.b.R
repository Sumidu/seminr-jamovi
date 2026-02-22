#' @importFrom R6 R6Class
#' @importFrom jmvcore Analysis
SeminrPLSClass <- if (requireNamespace("jmvcore", quietly = TRUE))
    R6::R6Class(
        "SeminrPLSClass",
        inherit = SeminrPLSBase,
        private = list(

            .init = function() {
                # Show a helpful HTML message about model setup
                html <- self$results$modelSpec
                if (!private$.hasValidSpec()) {
                    html$setContent(private$.buildInstructionsHtml())
                } else {
                    html$setContent(private$.buildModelSummaryHtml())
                }
            },

            .run = function() {
                if (!private$.hasValidSpec())
                    return()

                if (!requireNamespace("seminr", quietly = TRUE))
                    jmvcore::reject(
                        "The 'seminr' package is required for this analysis. ",
                        "Install it with: install.packages('seminr')")

                data     <- self$data
                options  <- self$options

                # Build measurement model
                mm <- private$.buildMeasurementModel(data)

                # Build structural model
                sm <- private$.buildStructuralModel()

                # Determine inner weighting scheme
                inner_weights <- if (options$innerWeights == "path_weighting")
                    seminr::path_weighting
                else
                    seminr::factor_weighting

                # Determine missing data handling
                missing_fn <- if (options$missing == "mean_replacement")
                    seminr::mean_replacement
                else
                    seminr::listwise_deletion

                # Estimate PLS model
                pls_model <- tryCatch(
                    seminr::estimate_pls(
                        data              = data,
                        measurement_model = mm,
                        structural_model  = sm,
                        inner_weights     = inner_weights,
                        missing           = missing_fn),
                    error = function(e) {
                        jmvcore::reject(paste0("PLS estimation failed: ", e$message))
                    })

                model_summary <- seminr::summary(pls_model)

                # Populate measurement model table
                private$.fillMeasurementTable(pls_model, model_summary, NULL)

                # Populate structural model table (without bootstrap)
                private$.fillStructuralTable(pls_model, model_summary, NULL)

                # Populate R-squared table
                private$.fillRSquaredTable(model_summary)

                # Populate model fit table
                private$.fillModelFitTable(model_summary)

                # Bootstrap if requested
                if (options$bootstrap) {
                    set.seed(options$seed)
                    boot_model <- tryCatch(
                        seminr::bootstrap_model(
                            seminr_model = pls_model,
                            nboot        = options$nboot,
                            cores        = 1L),
                        error = function(e) {
                            jmvcore::reject(paste0("Bootstrapping failed: ", e$message))
                        })

                    boot_summary <- seminr::summary(boot_model)

                    # Re-fill tables with bootstrap results
                    private$.fillMeasurementTable(pls_model, model_summary, boot_summary)
                    private$.fillStructuralTable(pls_model, model_summary, boot_summary)
                }
            },

            # -- Helper: check if the user has specified a valid model --------

            .hasValidSpec = function() {
                options <- self$options
                constructs <- options$constructs
                paths      <- options$paths

                if (length(constructs) == 0 || length(paths) == 0)
                    return(FALSE)

                # At least one construct must have a name and indicators
                has_construct <- any(vapply(constructs, function(c) {
                    nchar(trimws(c$name)) > 0 && length(c$indicators) > 0
                }, logical(1)))

                # At least one path must have from and to filled
                has_path <- any(vapply(paths, function(p) {
                    nchar(trimws(p$from)) > 0 && nchar(trimws(p$to)) > 0
                }, logical(1)))

                has_construct && has_path
            },

            # -- Helper: build seminr measurement model -----------------------

            .buildMeasurementModel = function(data) {
                options    <- self$options
                constructs <- options$constructs

                construct_list <- lapply(constructs, function(c) {
                    cname <- trimws(c$name)
                    ctype <- c$type
                    inds  <- c$indicators

                    if (nchar(cname) == 0 || length(inds) == 0)
                        return(NULL)

                    # Validate that all indicator columns exist in data
                    missing_cols <- setdiff(inds, colnames(data))
                    if (length(missing_cols) > 0)
                        jmvcore::reject(paste0(
                            "Variables not found in dataset: ",
                            paste(missing_cols, collapse = ", ")))

                    if (ctype == "reflective") {
                        seminr::reflective(cname, inds)
                    } else if (ctype == "composite_a") {
                        seminr::composite(cname, inds, weights = seminr::mode_A)
                    } else {
                        seminr::composite(cname, inds, weights = seminr::mode_B)
                    }
                })

                # Remove NULLs (incomplete constructs)
                construct_list <- Filter(Negate(is.null), construct_list)

                if (length(construct_list) == 0)
                    jmvcore::reject(
                        "No valid constructs defined. ",
                        "Please specify at least one construct with a name and indicators.")

                do.call(seminr::constructs, construct_list)
            },

            # -- Helper: build seminr structural model ------------------------

            .buildStructuralModel = function() {
                options <- self$options
                paths   <- options$paths

                path_list <- lapply(paths, function(p) {
                    pfrom <- trimws(p$from)
                    pto   <- trimws(p$to)

                    if (nchar(pfrom) == 0 || nchar(pto) == 0)
                        return(NULL)

                    seminr::paths(from = pfrom, to = pto)
                })

                path_list <- Filter(Negate(is.null), path_list)

                if (length(path_list) == 0)
                    jmvcore::reject(
                        "No valid structural paths defined. ",
                        "Please specify at least one path with From and To construct names.")

                do.call(seminr::relationships, path_list)
            },

            # -- Helper: fill measurement model table -------------------------

            .fillMeasurementTable = function(pls_model, model_summary, boot_summary) {
                table    <- self$results$measurement
                options  <- self$options
                constructs <- options$constructs

                # Clear existing rows
                table$deleteRows()

                # Get loadings/weights matrix from model summary
                # model_summary$loadings contains loadings for reflective
                # model_summary$weights contains weights for composite
                loadings <- model_summary$loadings
                weights  <- model_summary$weights

                row_idx <- 0L

                for (c in constructs) {
                    cname <- trimws(c$name)
                    ctype <- c$type
                    inds  <- c$indicators

                    if (nchar(cname) == 0 || length(inds) == 0)
                        next

                    type_label <- switch(ctype,
                        "reflective"  = "Reflective",
                        "composite_a" = "Composite (A)",
                        "composite_b" = "Composite (B)")

                    for (ind in inds) {
                        row_idx <- row_idx + 1L

                        # Get the loading or weight value
                        value <- NA_real_
                        if (ctype == "reflective" && !is.null(loadings)) {
                            if (ind %in% rownames(loadings) && cname %in% colnames(loadings))
                                value <- loadings[ind, cname]
                        } else if (!is.null(weights)) {
                            if (ind %in% rownames(weights) && cname %in% colnames(weights))
                                value <- weights[ind, cname]
                        }

                        row <- list(
                            construct = cname,
                            indicator = ind,
                            type      = type_label,
                            loading   = value)

                        # Add bootstrap results if available
                        if (!is.null(boot_summary)) {
                            # seminr boot_summary structure for loadings/weights
                            boot_row <- private$.getBootstrapMeasurement(
                                boot_summary, cname, ind, ctype)
                            row$se     <- boot_row$se
                            row$tstat  <- boot_row$tstat
                            row$pvalue <- boot_row$pvalue
                        }

                        table$addRow(rowKey = row_idx, values = row)
                    }
                }
            },

            # -- Helper: fill structural model table --------------------------

            .fillStructuralTable = function(pls_model, model_summary, boot_summary) {
                table   <- self$results$structural
                options <- self$options
                paths   <- options$paths

                table$deleteRows()

                # Path coefficients are in model_summary$paths
                path_coefs <- model_summary$paths

                row_idx <- 0L

                for (p in paths) {
                    pfrom <- trimws(p$from)
                    pto   <- trimws(p$to)

                    if (nchar(pfrom) == 0 || nchar(pto) == 0)
                        next

                    row_idx <- row_idx + 1L

                    value <- NA_real_
                    if (!is.null(path_coefs)) {
                        if (pfrom %in% rownames(path_coefs) &&
                            pto   %in% colnames(path_coefs))
                            value <- path_coefs[pfrom, pto]
                    }

                    row <- list(
                        from     = pfrom,
                        to       = pto,
                        estimate = value)

                    # Add bootstrap results if available
                    if (!is.null(boot_summary)) {
                        boot_row <- private$.getBootstrapPath(
                            boot_summary, pfrom, pto)
                        row$se       <- boot_row$se
                        row$tstat    <- boot_row$tstat
                        row$pvalue   <- boot_row$pvalue
                        row$ci_lower <- boot_row$ci_lower
                        row$ci_upper <- boot_row$ci_upper
                    }

                    table$addRow(rowKey = row_idx, values = row)
                }
            },

            # -- Helper: fill R-squared table ---------------------------------

            .fillRSquaredTable = function(model_summary) {
                table <- self$results$rSquared
                table$deleteRows()

                rsq <- model_summary$paths  # seminr stores R² info here too
                # Actually, R² is in model_summary$it_criteria or summary$r_squared
                # Try to extract from summary
                r2_vals <- tryCatch({
                    # seminr 2.x: summary has $r_squared
                    if (!is.null(model_summary$r_squared))
                        model_summary$r_squared
                    else if (!is.null(model_summary$paths))
                        # Fallback: try to get R² from the it_criteria
                        NULL
                    else
                        NULL
                }, error = function(e) NULL)

                if (is.null(r2_vals))
                    return()

                # r2_vals is typically a named matrix or vector
                if (is.matrix(r2_vals)) {
                    row_idx <- 0L
                    for (construct in colnames(r2_vals)) {
                        r2     <- r2_vals["R^2", construct]
                        r2_adj <- if ("R^2 Adj." %in% rownames(r2_vals))
                            r2_vals["R^2 Adj.", construct]
                        else
                            NA_real_

                        if (!is.na(r2) && r2 > 0) {
                            row_idx <- row_idx + 1L
                            table$addRow(rowKey = row_idx, values = list(
                                construct = construct,
                                rsq       = r2,
                                rsq_adj   = r2_adj))
                        }
                    }
                }
            },

            # -- Helper: fill model fit table ---------------------------------

            .fillModelFitTable = function(model_summary) {
                table <- self$results$modelFit

                srmr      <- NA_real_
                rms_theta <- NA_real_
                nfl       <- NA_real_

                # seminr 2.x: model_summary$validity$fl_criteria or similar
                # model fit is in model_summary$validity$model_fit
                tryCatch({
                    fit <- model_summary$validity$model_fit
                    if (!is.null(fit)) {
                        srmr      <- fit$SRMR
                        rms_theta <- fit$RMS_theta
                    }
                }, error = function(e) {})

                # Some versions use $it_criteria for fit
                tryCatch({
                    if (is.na(srmr) && !is.null(model_summary$it_criteria)) {
                        ic <- model_summary$it_criteria
                        srmr <- ic["SRMR", 1]
                    }
                }, error = function(e) {})

                table$setRow(rowKey = 1, values = list(
                    srmr      = srmr,
                    rms_theta = rms_theta,
                    nfl       = nfl))
            },

            # -- Helper: extract bootstrap measurement results ----------------

            .getBootstrapMeasurement = function(boot_summary, construct, indicator, type) {
                result <- list(se = NA_real_, tstat = NA_real_, pvalue = NA_real_)

                tryCatch({
                    # seminr bootstrap summary has $bootstrapped_loadings or
                    # $bootstrapped_weights
                    if (type == "reflective") {
                        bt <- boot_summary$bootstrapped_loadings
                    } else {
                        bt <- boot_summary$bootstrapped_weights
                    }

                    if (!is.null(bt)) {
                        # bt is typically a 3D array or list
                        # Row format: c(indicator, construct, ...)
                        key <- paste0(indicator, " -> ", construct)
                        if (!is.null(bt[key, "SE", drop = FALSE])) {
                            result$se    <- bt[key, "SE"]
                            result$tstat <- bt[key, "T Stat."]
                            result$pvalue <- private$.tstat_to_pvalue(bt[key, "T Stat."])
                        }
                    }
                }, error = function(e) {})

                result
            },

            # -- Helper: extract bootstrap path results -----------------------

            .getBootstrapPath = function(boot_summary, from, to) {
                result <- list(
                    se       = NA_real_,
                    tstat    = NA_real_,
                    pvalue   = NA_real_,
                    ci_lower = NA_real_,
                    ci_upper = NA_real_)

                tryCatch({
                    bt <- boot_summary$bootstrapped_paths
                    if (!is.null(bt)) {
                        key <- paste0(from, " -> ", to)
                        if (key %in% rownames(bt)) {
                            result$se     <- bt[key, "SE"]
                            result$tstat  <- bt[key, "T Stat."]
                            result$pvalue <- private$.tstat_to_pvalue(bt[key, "T Stat."])
                            if ("2.5% CI" %in% colnames(bt))
                                result$ci_lower <- bt[key, "2.5% CI"]
                            if ("97.5% CI" %in% colnames(bt))
                                result$ci_upper <- bt[key, "97.5% CI"]
                        }
                    }
                }, error = function(e) {})

                result
            },

            # -- Helper: convert t-statistic to two-tailed p-value -----------

            .tstat_to_pvalue = function(tstat) {
                if (is.na(tstat)) return(NA_real_)
                2 * pt(-abs(tstat), df = Inf)
            },

            # -- Helper: build HTML instructions ------------------------------

            .buildInstructionsHtml = function() {
                paste0(
                    "<p><b>PLS-SEM Setup Instructions</b></p>",
                    "<p>To run a PLS-SEM analysis, you need to:</p>",
                    "<ol>",
                    "<li><b>Define Constructs:</b> For each latent construct, ",
                    "specify a name, measurement type (reflective or composite), ",
                    "and select the indicator variables that measure it.</li>",
                    "<li><b>Define Structural Paths:</b> Specify the hypothesized ",
                    "paths between constructs using the exact construct names ",
                    "defined above (e.g., 'Image' &rarr; 'Loyalty').</li>",
                    "<li><b>Set Options:</b> Choose the inner weighting scheme ",
                    "and whether to bootstrap for inference.</li>",
                    "</ol>",
                    "<p><i>Note: The seminr package must be installed. ",
                    "Use <code>install.packages('seminr')</code> to install it.</i></p>")
            },

            # -- Helper: build HTML model summary ----------------------------

            .buildModelSummaryHtml = function() {
                options    <- self$options
                constructs <- options$constructs
                paths      <- options$paths

                # Build constructs summary
                c_rows <- ""
                for (c in constructs) {
                    cname <- trimws(c$name)
                    if (nchar(cname) == 0) next
                    type_label <- switch(c$type,
                        "reflective"  = "Reflective",
                        "composite_a" = "Composite (Mode A)",
                        "composite_b" = "Composite (Mode B)",
                        c$type)
                    inds <- paste(c$indicators, collapse = ", ")
                    c_rows <- paste0(c_rows,
                        "<tr><td><b>", cname, "</b></td>",
                        "<td>", type_label, "</td>",
                        "<td>", inds, "</td></tr>")
                }

                # Build paths summary
                p_rows <- ""
                for (p in paths) {
                    pfrom <- trimws(p$from)
                    pto   <- trimws(p$to)
                    if (nchar(pfrom) == 0 || nchar(pto) == 0) next
                    p_rows <- paste0(p_rows,
                        "<tr><td>", pfrom, "</td>",
                        "<td>&rarr;</td>",
                        "<td>", pto, "</td></tr>")
                }

                paste0(
                    "<b>Measurement Model</b>",
                    "<table border='1' cellpadding='4' style='border-collapse:collapse'>",
                    "<tr><th>Construct</th><th>Type</th><th>Indicators</th></tr>",
                    c_rows,
                    "</table><br>",
                    "<b>Structural Model</b>",
                    "<table border='1' cellpadding='4' style='border-collapse:collapse'>",
                    "<tr><th>From</th><th></th><th>To</th></tr>",
                    p_rows,
                    "</table>")
            }

        ) # end private
    ) # end R6Class
