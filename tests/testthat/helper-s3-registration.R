# Readers for R's S3 registration tables. Two files pin where a borrowed method
# comes from, test-dependencies.R for the printers and the shared weight-vector
# methods and test-ipw-accessors.R for the result accessors, and both ask the
# same question of the same tables.

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

# Every method in one table whose name matches `pattern`, named by the package
# that defined it. A method under a generic a spec does not name is still
# recorded in the table belonging to that generic's package, so a spec ruling
# out any registration at all reads the table whole rather than only at the
# names it expects.
registered_sources <- function(where, pattern) {
  table <- get(".__S3MethodsTable__.", envir = where)
  names <- grep(pattern, ls(table, all.names = TRUE), value = TRUE)
  vapply(names, method_source, character(1), where = where)
}

# `UseMethod()` searches the environment its generic was called from as well as
# the registration table, so a method left behind in balancing's namespace would
# still shadow an inherited one after its registration was removed. Every
# absence claim has to be checked in both places.
defined_in_balancing <- function(name) {
  exists(name, envir = asNamespace("balancing"), inherits = FALSE)
}
