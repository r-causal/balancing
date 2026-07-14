# .onLoad() registers the S7 methods so that generics dispatch on the package's
# classes after the namespace loads. Registration is idempotent, so calling it
# again is safe and lets the test exercise the load hook directly.

test_that(".onLoad() registers S7 methods and returns invisibly", {
  expect_invisible(.onLoad("libname", "balancing"))
})

test_that("S7 generics dispatch on balancing results after loading", {
  data <- sim_binary(n = 150)
  fit <- balance(
    data,
    exposure,
    c(x1, x2),
    method = bal_entropy(),
    estimand = "ate"
  )
  # print() dispatches through the registered S7 method rather than erroring.
  expect_output(print(fit), "Entropy balancing")
})
