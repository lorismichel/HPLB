test_that("inputs of function dWitMatrix are checked correctly", {

  # data inputs
  ordering <- array(1,dim = c(2,2,2))
  labels <- c(1,2)

  # checking seed
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, seed = "seed"), "seed")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, seed = c(0,1)), "seed")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, seed = NA), "seed")

  # checking labels and ordering.array
  expect_error(dWitMatrix(labels = c(1), ordering.array = ordering, seed = 0), "labels")
  expect_error(dWitMatrix(labels = c(1,NA), ordering.array = ordering, seed = 0), "labels")
  expect_error(dWitMatrix(labels = c(1,2,3), ordering.array = ordering, seed = 0), "labels")
  expect_error(dWitMatrix(labels = matrix(c(1,2,3,4),ncol=2), ordering.array = ordering, seed = 0), "labels")

  expect_error(dWitMatrix(labels = labels, ordering.array = ordering[,,1], seed = 0), "ordering.array")
  expect_error(dWitMatrix(labels = labels, ordering.array = NA, seed = 0), "ordering.array")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering[1:2,,], seed = 0), "ordering.array")

  # checking alpha
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, alpha = -0.5), "alpha")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, alpha = 1.5), "alpha")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, alpha = c(1,2)), "alpha")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, alpha = "a"), "alpha")

  # checking computation.type
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, computation.type = "a"), "computation.type")
  expect_error(dWitMatrix(labels = labels, ordering.array = ordering, computation.type = NA), "computation.type")
})
