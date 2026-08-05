# balancing does not define every S3 method its objects dispatch to. An `ipw`
# result prints and tabulates through causalgenerics, and a `bw` vector reaches
# the shared `causal_wts` methods for math, summaries, and subsetting. A
# borrowed method is invisible to a search of the sources, so these specs pin
# where each one is registered, which packages balancing imports from, and,
# from a fresh R process, that a whole `ipw()` workflow runs with propensity
# absent from the loaded namespaces.

# ---- Reading the S3 registration tables -----------------------------------

# R records a registered S3 method in the `.__S3MethodsTable__.` of the
# environment where its generic is defined, so reading that table names the
# package a method actually comes from. `getS3method()` is not a substitute: it
# returns `NULL` for a generic that is not visible on the search path, which
# under `R CMD check` would make an absence assertion pass without testing
# anything.
registered_method <- function(name, where) {
  table <- get(".__S3MethodsTable__.", envir = where)
  if (!exists(name, envir = table, inherits = FALSE)) {
    return(NULL)
  }
  get(name, envir = table, inherits = FALSE)
}

# The package a registered method was defined in, or `NA` when the table holds
# no method by that name.
method_source <- function(name, where) {
  method <- registered_method(name, where)
  if (is.null(method)) {
    return(NA_character_)
  }
  environmentName(environment(method))
}

# `UseMethod()` searches the environment its generic was called from as well as
# the registration table, so a method left behind in balancing's namespace would
# still shadow an inherited one after its registration was removed. Every
# absence claim has to be checked in both places.
defined_in_balancing <- function(name) {
  exists(name, envir = asNamespace("balancing"), inherits = FALSE)
}

# The symbols balancing imports from `package`, read from the namespace's own
# import record rather than from DESCRIPTION, because that record is the edge
# that loads the other namespace. A package imported from in more than one
# roxygen block appears once per block, so the entries are pooled.
imported_from <- function(package) {
  imports <- getNamespaceImports(asNamespace("balancing"))
  sort(unique(unlist(lapply(imports[names(imports) == package], names))))
}

# The packages listed in one DESCRIPTION dependency field, with any version
# constraint stripped.
declared_packages <- function(field) {
  value <- utils::packageDescription("balancing", fields = field)
  if (length(value) != 1L || is.na(value)) {
    return(character())
  }
  entries <- trimws(strsplit(value, ",")[[1]])
  trimws(sub("\\(.*\\)$", "", entries[nzchar(entries)]))
}

# ---- Running balancing in a fresh R process -------------------------------

# Which namespaces loading balancing pulls in can only be observed from outside
# this session, because the suite itself uses propensity for the psw interop
# specs: `loadedNamespaces()` here reports propensity no matter what balancing
# declares. Checking the search path instead would prove nothing either, since a
# package in `Imports` is loaded without being attached.

# The copy of balancing under test is an installed one under `R CMD check` and a
# source directory under `devtools::test()`. The subprocess has to load that
# same copy; loading whichever version happens to sit in the library would
# report on a different package.
balancing_loader <- function() {
  path <- getNamespaceInfo("balancing", "path")
  if (file.exists(file.path(path, "Meta", "package.rds"))) {
    return(sprintf("library(balancing, lib.loc = %s)", deparse1(dirname(path))))
  }
  skip_if_not_installed("pkgload")
  sprintf(
    paste0(
      "pkgload::load_all(%s, compile = FALSE, export_all = FALSE, ",
      "helpers = FALSE, attach_testthat = FALSE, quiet = TRUE)"
    ),
    deparse1(path)
  )
}

# Run R code in a fresh process that starts from this session's library paths
# and reads no startup files, and return its exit status alongside the combined
# standard output and error.
run_fresh_r <- function(code) {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", deparse1(.libPaths())), code), script)
  # `R CMD check` points R_TESTS at a startup file that belongs to the session
  # it launched; leaving it set makes every nested R process fail at startup.
  withr::local_envvar(R_TESTS = "")
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)),
    stdout = TRUE,
    stderr = TRUE
  ))
  list(status = as.integer(attr(output, "status") %||% 0L), output = output)
}

# The subprocess reports each fact as a `name: value` line so an assertion can
# name the one thing it reads instead of matching against a whole transcript.
# A field the subprocess never reported says nothing on its own, and comparing
# the absence against the expected value would report only that one side was
# missing, so the transcript is what gets shown instead.
subprocess_field <- function(run, name) {
  pattern <- paste0("^", name, ": ")
  line <- grep(pattern, run$output, value = TRUE)
  if (length(line) != 1L) {
    fail(paste0(
      "The subprocess reported ",
      length(line),
      " lines for the field \"",
      name,
      "\", expected one:\n",
      paste(run$output, collapse = "\n")
    ))
    return(NA_character_)
  }
  sub(pattern, "", line)
}

# ---- The ipw() workflow without propensity --------------------------------

test_that("the ipw() workflow runs with propensity absent from the namespaces", {
  run <- run_fresh_r(c(
    balancing_loader(),
    # Two readings of the same fact bracket the workflow. This one is taken
    # before any balancing code runs, so it reports only what loading the
    # namespace pulled in on its own.
    'cat(sprintf("propensity_at_load: %s\\n", "propensity" %in% loadedNamespaces()))',
    'options(balancing.quiet = TRUE)',
    'set.seed(101)',
    'n <- 200',
    'x1 <- rnorm(n)',
    'x2 <- rnorm(n)',
    'z <- rbinom(n, 1L, plogis(0.7 * x1 - 0.5 * x2))',
    'y <- rbinom(n, 1L, plogis(-0.3 + 0.5 * z + 0.4 * x1))',
    'data <- data.frame(exposure = z, x1 = x1, x2 = x2, y = y)',
    paste0(
      'fit <- balance(data, exposure, c(x1, x2), ',
      'method = bw_entropy(), estimand = "ate")'
    ),
    'data$.wts <- as.numeric(weights(fit))',
    paste0(
      'model <- suppressWarnings(glm(y ~ exposure, data = data, ',
      'family = binomial(), weights = .wts))'
    ),
    'result <- ipw(fit, model)',
    'printed <- capture.output(print(result))',
    'estimates <- as.data.frame(result)',
    # Fitting, weighting, and estimating never coerce, combine, or take the
    # prototype of a `bw` vector, so the workflow above reaches none of the
    # vctrs methods the class registers. Those methods are where a weight
    # vector's estimand is read, and that read is the last place balancing can
    # reach back into propensity, so a namespace assertion that never ran them
    # would report a clean namespace with the dependency still in place. The
    # path is a mainstream one rather than a corner of the class: composing
    # sampling weights onto fitted weights restores through it, and so does
    # every arithmetic operation on the result.
    #
    # The second reading of `loadedNamespaces()` comes after this block, and
    # the pair is what tells the two surviving forms of the dependency apart.
    # Both readings true means the namespace still declares propensity and
    # loading balancing pulls it in. A false followed by a true means the
    # declaration is gone but the estimand reads below still resolve into it,
    # which is the form that leaves every other assertion here satisfied.
    'w <- weights(fit)',
    'invisible(w * 2)',
    'invisible(w + w)',
    'invisible(cumsum(w))',
    'invisible(vctrs::vec_ptype_abbr(w))',
    'invisible(vctrs::vec_ptype_full(w))',
    'invisible(vctrs::vec_ptype2(w, w))',
    'invisible(vctrs::vec_cast(1, w))',
    'invisible(vctrs::vec_cast(1L, w))',
    'cat(sprintf("print_header: %s\\n", printed[[1]]))',
    'cat(sprintf("effects: %s\\n", paste(estimates$effect, collapse = ",")))',
    'cat(sprintf("propensity_after_use: %s\\n", "propensity" %in% loadedNamespaces()))',
    'cat(sprintf("causalgenerics: %s\\n", "causalgenerics" %in% loadedNamespaces()))'
  ))

  # A failing expectation does not end the block, so reading the fields after a
  # subprocess that never got far enough to report them would repeat the same
  # transcript once per field on top of the status failure.
  if (!identical(run$status, 0L)) {
    fail(paste0(
      "The subprocess exited with status ",
      run$status,
      ":\n",
      paste(run$output, collapse = "\n")
    ))
  } else {
    expect_identical(subprocess_field(run, "propensity_at_load"), "FALSE")
    expect_identical(subprocess_field(run, "propensity_after_use"), "FALSE")
    expect_identical(subprocess_field(run, "causalgenerics"), "TRUE")
    expect_identical(
      subprocess_field(run, "print_header"),
      "Inverse Probability Weight Estimator"
    )
    expect_identical(subprocess_field(run, "effects"), "rd,log(rr),log(or)")
  }
})

# ---- The declared dependency surface --------------------------------------

test_that("the balancing namespace imports nothing from propensity", {
  imports <- names(getNamespaceImports(asNamespace("balancing")))
  expect_false("propensity" %in% imports)
  expect_true("causalgenerics" %in% imports)
})

test_that("the causal-weight generics are imported from causalgenerics", {
  expect_true(
    all(
      c("estimand", "ipw", "is_causal_wt") %in% imported_from("causalgenerics")
    )
  )
})

test_that("propensity is a suggested package rather than an imported one", {
  expect_false("propensity" %in% declared_packages("Imports"))
  expect_true("propensity" %in% declared_packages("Suggests"))
  expect_true("causalgenerics" %in% declared_packages("Imports"))
})

# ---- Where the borrowed methods come from ---------------------------------

test_that("the exported causal-weight generics are the causalgenerics ones", {
  expect_identical(estimand, causalgenerics::estimand)
  expect_identical(is_causal_wt, causalgenerics::is_causal_wt)
  expect_identical(ipw, causalgenerics::ipw)
})

# The two generics that move an `ipw` result between its readings belong to
# causalgenerics, which owns the result class and the `effects` field they set.
# They join the re-export block for the reason `ipw()` is in it: a user who has
# attached balancing alone writes `as_conditional(result)` on a result
# `ipw()` built without a second attachment.
test_that("the presentation mode generics are exported and imported", {
  expect_contains(
    getNamespaceExports("balancing"),
    c("as_marginal", "as_conditional")
  )
  expect_true(
    all(c("as_marginal", "as_conditional") %in% imported_from("causalgenerics"))
  )
})

test_that("the exported mode generics are the causalgenerics ones", {
  expect_identical(as_marginal, causalgenerics::as_marginal)
  expect_identical(as_conditional, causalgenerics::as_conditional)
})

test_that("the ipw result methods are the ones causalgenerics registers", {
  expect_identical(method_source("print.ipw", baseenv()), "causalgenerics")
  expect_identical(
    method_source("as.data.frame.ipw", baseenv()),
    "causalgenerics"
  )
  expect_false(defined_in_balancing("print.ipw"))
  expect_false(defined_in_balancing("as.data.frame.ipw"))
})

# `UseMethod()` takes the first match down the class vector, so a method left
# registered on `bw` shadows the `causal_wts` one completely and the inherited
# implementation is never reached. Each of these has to be absent from both the
# registration table and the namespace for the shared method to run.
test_that("the shared causal_wts methods are not registered on bw", {
  specs <- list(
    vec_math.bw = asNamespace("vctrs"),
    Summary.bw = baseenv(),
    min.bw = baseenv(),
    max.bw = baseenv(),
    `[.bw` = baseenv(),
    median.bw = asNamespace("stats"),
    quantile.bw = asNamespace("stats")
  )
  sources <- vapply(
    names(specs),
    function(name) method_source(name, specs[[name]]),
    character(1)
  )
  expect_identical(
    sources,
    stats::setNames(rep(NA_character_, length(specs)), names(specs))
  )
})

test_that("no shadowing copy of a shared method survives in the namespace", {
  names <- c(
    "vec_math.bw",
    "Summary.bw",
    "min.bw",
    "max.bw",
    "[.bw",
    "median.bw",
    "quantile.bw"
  )
  defined <- vapply(names, defined_in_balancing, logical(1))
  expect_identical(
    defined,
    stats::setNames(rep(FALSE, length(names)), names)
  )
})

# The other half of the split. `vec_ptype2`, `vec_cast`, `vec_ptype_abbr`, and
# `vec_ptype_full` are resolved by vctrs on the first element of the class
# vector alone, so they cannot be inherited from `causal_wts` at all;
# `vec_restore` and `vec_arith` could be, but only the concrete class knows that
# `groups` is index-typed and that combining a `bw` with a `psw` must stay an
# error. All of them stay balancing's own.
test_that("the concrete vctrs methods stay registered on bw", {
  specs <- c(
    "vec_restore.bw",
    "vec_ptype_abbr.bw",
    "vec_ptype_full.bw",
    "vec_arith.bw",
    "vec_arith.numeric.bw",
    "vec_ptype2.bw.bw",
    "vec_ptype2.bw.double",
    "vec_ptype2.double.bw",
    "vec_ptype2.bw.integer",
    "vec_ptype2.integer.bw",
    "vec_ptype2.bw.character",
    "vec_ptype2.character.bw",
    "vec_ptype2.bw.psw",
    "vec_ptype2.psw.bw",
    "vec_cast.bw.bw",
    "vec_cast.bw.double",
    "vec_cast.double.bw",
    "vec_cast.bw.integer",
    "vec_cast.integer.bw",
    "vec_cast.character.bw"
  )
  sources <- vapply(
    specs,
    method_source,
    character(1),
    where = asNamespace("vctrs")
  )
  expect_identical(
    sources,
    stats::setNames(rep("balancing", length(specs)), specs)
  )
})

# `vec_arith.bw` is itself a generic balancing defines, so its methods are
# recorded in balancing's own table rather than in vctrs'.
test_that("the vec_arith.bw methods stay registered in balancing's table", {
  specs <- c(
    "vec_arith.bw.bw",
    "vec_arith.bw.default",
    "vec_arith.bw.integer",
    "vec_arith.bw.numeric",
    "vec_arith.bw.MISSING"
  )
  sources <- vapply(
    specs,
    method_source,
    character(1),
    where = asNamespace("balancing")
  )
  expect_identical(
    sources,
    stats::setNames(rep("balancing", length(specs)), specs)
  )
})

# ---- What the inherited methods must still do -----------------------------

# The shared implementations have to produce exactly what balancing's own did,
# on a weight vector carrying the metadata a fitted vector carries: reductions
# drop to a plain numeric, cumulative math keeps the class, and subsetting keeps
# the estimand.
test_that("the inherited causal_wts methods reproduce the previous behavior", {
  w <- new_bw(
    c(3, 1, 2, 4),
    estimand = "ate",
    groups = list(`0` = 1:2, `1` = 3:4)
  )

  expect_false(is_bw(sum(w)))
  expect_equal(sum(w), 10)
  expect_equal(min(w), 1)
  expect_equal(max(w), 4)
  expect_equal(range(w), c(1, 4))
  expect_false(is_bw(min(w)))

  expect_equal(median(w), 2.5)
  expect_false(is_bw(median(w)))
  expect_equal(
    quantile(w, probs = c(0.25, 0.75)),
    quantile(c(3, 1, 2, 4), probs = c(0.25, 0.75))
  )
  expect_false(is_bw(quantile(w)))

  running <- cumsum(w)
  expect_true(is_bw(running))
  expect_equal(vctrs::vec_data(running), c(3, 4, 6, 10))
  expect_false(is_bw(sqrt(w)))

  first_two <- w[1:2]
  expect_true(is_bw(first_two))
  expect_identical(estimand(first_two), "ate")
  expect_equal(vctrs::vec_data(first_two), c(3, 1))
  expect_true(is_bw(w[]))
  expect_false(is_bw(w[matrix(c(1L, 3L), ncol = 1)]))
})
