feat <- readRDS(test_path("fixtures/unprocessed/rel_ab_table.rds"))
feat <- feat[apply(feat, 1, sd) != 0,]
log.n0 <- 1e-6
sd.min.q <- 0.1

test_that("log.clr sums to 0", {
  expect_true(all.equal(unname(colSums(log.clr(feat, 1e-3))), rep(0, ncol(feat))))
})

test_that("log.clr values", {
  actual <- log.clr(feat, log.n0)
  expected <- readRDS(test_path("fixtures/normalize_features/log_clr.rds"))
  # expected is in natural log, but I switched to log10 for consistency
  # with other normalisations
  expected <- expected / log(10)
  expect_equal(actual, expected)
})

test_that("rank.unit values", {
  actual <- rank.unit(as.matrix(feat))
  expected <- readRDS(test_path("fixtures/normalize_features/rank_unit.rds"))
  expect_equal(actual, expected)
})

test_that("log.std values", {
  actual <- log.std(as.matrix(feat), log.n0, par=list(), sd.min.q=sd.min.q)$feat.norm
  expected <- readRDS(test_path("fixtures/normalize_features/log_std.rds"))
  expect_equal(actual, expected)
})

test_that("rank.std values", {
  actual <- rank.std(as.matrix(feat), par=list(), sd.min.q=sd.min.q)$feat.norm
  expected <- readRDS(test_path("fixtures/normalize_features/rank_std.rds"))
  expect_equal(actual, expected)
})