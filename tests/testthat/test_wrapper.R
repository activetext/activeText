# Tests active_label() wrapper.


set.seed(9937) ## From random.org, 2024/01/01

## ---------
## Run tests
## ---------
test_that("regression test", {

  mockr::local_mock(classification_gui = function(options, documents, param_to_save, selections) {
    # return (c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1));
    return (c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2));
  })

  expect_equal(active_label(docs = activeText::bbc_data_all,
                    labels = c("Not Political", "Political"), seed=9937), readRDS("test_comparison_all_pol.RDS"))

})
