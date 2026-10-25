library(scMMR)

local({
  failures <- character()
  checks <- 0L
  check <- function(label, code) {
    checks <<- checks + 1L
    error <- tryCatch({ force(code); NULL }, error = conditionMessage)
    if (!is.null(error)) {
      failures <<- c(failures, paste(label, error, sep = ": "))
    }
  }

  d <- data.frame(y = 1:8, x = c(1, 2, 3, 5, 4, 6, 8, 7))
  check("Kendall coefficient and p-value", {
    result <- RunCorrelation(d, "y", method = "kendall",
                             min.cells = 3, verbose = FALSE)
    reference <- stats::cor.test(d$y, d$x, method = "kendall", exact = FALSE)
    stopifnot(isTRUE(all.equal(result$score, unname(reference$estimate))),
              isTRUE(all.equal(result$pvalue, reference$p.value)))
  })

  d <- data.frame(y = 1:8, a = c(NA, 1:7),
                  b = c(8, 5, NA, 3, 2, 4, 1, 0))
  for (method in c("pearson", "spearman", "kendall")) {
    check(paste(method, "uses complete pairs"), {
      result <- RunCorrelation(d, "y", method = method,
                               min.cells = 3, verbose = FALSE)
      for (feature in c("a", "b")) {
        complete <- is.finite(d$y) & is.finite(d[[feature]])
        reference <- stats::cor.test(d$y[complete], d[[feature]][complete],
                                     method = method, exact = FALSE)
        row <- result[result$gene == feature, ]
        stopifnot(nrow(row) == 1L, row$n_cells == sum(complete),
                  isTRUE(all.equal(row$score, unname(reference$estimate))),
                  isTRUE(all.equal(row$pvalue, reference$p.value)))
      }
    })
  }

  check("Kendall p-values account for ties", {
    d <- data.frame(y = 1:8, x = c(1, 1, 2, 2, 4, 4, 7, 8))
    result <- RunCorrelation(d, "y", method = "kendall",
                             min.cells = 3, verbose = FALSE)
    reference <- stats::cor.test(d$y, d$x, method = "kendall", exact = FALSE)
    stopifnot(isTRUE(all.equal(result$score, unname(reference$estimate))),
              isTRUE(all.equal(result$pvalue, reference$p.value)))
  })

  check("non-finite values are excluded pairwise", {
    d <- data.frame(y = c(1:7, Inf), x = c(NA, 1:5, Inf, 7))
    result <- RunCorrelation(d, "y", method = "spearman",
                             min.cells = 3, verbose = FALSE)
    stopifnot(result$n_cells == 5L, abs(result$score - 1) < 1e-12)
  })

  check("numeric group labels are names rather than positions", {
    d <- data.frame(y = rep(1:8, 2), x = rep(c(1, 2, 3, 5, 4, 6, 8, 7), 2),
                    g = rep(c(10, 20), each = 8))
    result <- RunCorrelation(d, "y", "x", group.by = "g",
                             min.cells = 3, verbose = FALSE)
    stopifnot(setequal(result$group, c("10", "20")),
              all(result$n_cells == 8L))
    d$g <- factor(d$g, levels = c(20, 10, 30))
    factor_result <- RunCorrelation(d, "y", "x", group.by = "g",
                                    min.cells = 3, verbose = FALSE)
    stopifnot(setequal(factor_result$group, c("10", "20")),
              all(factor_result$n_cells == 8L))
  })

  if (requireNamespace("data.table", quietly = TRUE)) {
    check("data.table input follows the data.frame path", {
      d <- data.frame(y = 1:8, x = c(1, 2, 3, 5, 4, 6, 8, 7))
      expected <- RunCorrelation(d, "y", min.cells = 3, verbose = FALSE)
      result <- RunCorrelation(data.table::as.data.table(d), "y",
                               min.cells = 3, verbose = FALSE)
      stopifnot(identical(result, expected))
    })
  }

  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    check("min.pct filters genes across cells", {
      expression <- rbind(target = 1:30, rare = c(1, rep(0, 29)),
                          common = rep(c(1, 3, 2), 10))
      colnames(expression) <- paste0("cell", 1:30)
      object <- SeuratObject::CreateSeuratObject(
        counts = Matrix::Matrix(expression, sparse = TRUE)
      )
      result <- RunCorrelation(object, "target", c("rare", "common"),
                               layer = "counts", min.pct = 0.5,
                               min.cells = 3, verbose = FALSE)
      stopifnot(identical(result$gene, "common"))
    })
  }

  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)
  cat("RunCorrelation:", checks, "regression cases passed (8-30 rows).\n")
})
