##
##
## run the tests.

tests <- c("testGenominatorFxs.R", "testGoodnessOfFit.R", "testHasOverlap.R",
           "testImport.R")

sapply(tests, function(nm) {
    print("-------------------------------------------------------------------------------")
    print(paste("-------------------------------  Running Test:", nm, "-----------"))
    print("-------------------------------------------------------------------------------")
    source(nm)
})
