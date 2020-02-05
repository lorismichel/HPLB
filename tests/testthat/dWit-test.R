test_that("inputs of function dWit are checked correctly", {

  # checking seed
  expect_error(dWit(t = c(1,2), rho = c(1,2), seed = "seed"), "seed")
  expect_error(dWit(t = c(1,2), rho = c(1,2), seed = c(0,1)), "seed")
  expect_error(dWit(t = c(1,2), rho = c(1,2), seed = NA), "seed")

  # checking t
  expect_error(dWit(t = c(1,NA), rho = c(1,2)), "t")
  expect_error(dWit(t = c("a", "b"), rho = c(1,2)), "t")
  expect_error(dWit(t = data.frame(t=c(1,2)), rho = c(1,2)), "t")

  # checking s
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5,NA)), "s")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5,3)), "s")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5,-0.5)), "s")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c("a", "b")), "s")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = data.frame(s=c(1.5,1.6)), "s"))

  # checking z
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1,2), z = c(1,NA)), "z")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1,2), z = c(1,3)), "z")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1,2), z = c(1,-0.5)), "z")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1,2), z = c("a", "b")), "z")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1,2), z = data.frame(s=c(1.5,1.6)), "z"))

  # checking rho
  expect_error(dWit(t = c(1,2), rho = c(1.5,NA)), "rho")
  expect_error(dWit(t = c(1,2), rho = c(1.5,3)), "rho")
  expect_error(dWit(t = c(1,2), rho = c(1.5,-0.5)), "rho")
  expect_error(dWit(t = c(1,2), rho = c("a", "b")), "rho")
  expect_error(dWit(t = c(1,2), rho = data.frame(rho=c(1.5,1.6)), "rho"))

  # checking dimensions between rho, t and s
  expect_error(dWit(t = c(1,2,3), rho = matrix(c(1,2,3,4),ncol=2), s = c(1,2)), "dimension")
  expect_error(dWit(t = c(1,2), rho = matrix(c(1,2,3,4),ncol=2), s = c(1)), "dimension")
  expect_error(dWit(t = c(1,2), rho = c(1,2,3), s = c(1)), "dimension")

  # checking estimator.type
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), estimator.type = 1), "estimator")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), estimator.type = c("a")), "estimator")

  # checking alpha
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), alpha = -0.5), "alpha")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), alpha = 1.5), "alpha")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), alpha = c(1,2)), "alpha")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), alpha = "a"), "alpha")

  # checking tv.seq
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), tv.seq = c(NA,0.1)), "tv")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), tv.seq = c(-0.1,0.1)), "tv")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), tv.seq = data.frame(tv=c(0.1,0.2))), "tv")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), tv.seq = c(0.1, 0.1)), "tv")

  # checking direction
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), direction = c(NA)), "direction")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), direction = c(1)), "direction")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5), direction = c("a")), "direction")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5,1.7), direction = c("left", "a")),
               "direction")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5,1.7), direction = c("left")),
               "direction")
  expect_error(dWit(t = c(1,2), rho = c(1,2), s = c(1.5,1.7), direction = c("left", "a")),
               "direction")

  # checking user-defined bounding functions
  # TODO: implement these checks
})
