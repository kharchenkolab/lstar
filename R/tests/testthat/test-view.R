# view() is a soft delegate to the pagoda3 viewer (one-way optional dependency).

test_that("view() is exported and a function", {
  expect_true(is.function(view))
})

test_that("view() errors with an install hint when pagoda3 is absent", {
  local_mocked_bindings(.pagoda3_installed = function() FALSE)
  expect_error(view("x.lstar.zarr"), "pagoda3")
  expect_error(view("x.lstar.zarr"), "install.packages")
})

test_that("view() forwards to pagoda3 when present", {
  # When pagoda3 is installed, view() must delegate to pagoda3::view. An unviewable object (integer)
  # is rejected by pagoda3's coercion BEFORE any server starts, so the error proves the call reached
  # pagoda3 — and is NOT our own "install pagoda3" hint.
  skip_if_not_installed("pagoda3")
  err <- tryCatch(view(42L), error = function(e) conditionMessage(e))
  expect_false(grepl("install.packages", err))
})
