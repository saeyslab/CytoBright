test_that("No ISAC recommended short names are altered", {
  url <- "https://raw.githubusercontent.com/ISAC-DSTF/ProbeTagDictionary/master/ProbeTagDictionary.json"
  ProbeTagDictionary <- rjson::fromJSON(file = url)

  shortNames <- lapply(ProbeTagDictionary[[1]], function(x){
    sapply(x[[1]], function(y) y$`Preferred Short Name`)
  })
  names(shortNames) <- sapply(ProbeTagDictionary[[1]], names)

  for(type in names(shortNames)){
    cleaned <- clean_fluorochromes(shortNames[[type]])
    changed <- cleaned != shortNames[[type]]
    if(sum(changed) > 0)
      print(data.frame(original = shortNames[[type]],
                       cleaned = cleaned)[changed,])

    expect_equal(sum(changed), 0)
  }
})

test_that("BV variations", {
  # Dashes
  expect_equal(clean_fluorochromes("BV-421"),
               clean_fluorochromes("BV421"))
  # Capitalization
  expect_equal(clean_fluorochromes("Bv-421"),
               clean_fluorochromes("BV421"))
  # Abbreviated
  expect_equal(clean_fluorochromes("Brilliant violet 421"),
               clean_fluorochromes("BV421"))

})

test_that("Pacific", {
  expect_equal(clean_fluorochromes("PacBlue"),
               clean_fluorochromes("Pacific Blue"))

  expect_equal(clean_fluorochromes("PacBlue"),
               clean_fluorochromes("Pacific blue"))

  expect_equal(clean_fluorochromes("PacOrange"),
               clean_fluorochromes("Pacific Orange"))
})
