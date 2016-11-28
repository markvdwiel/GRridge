

testthat::test_that("Input for a partiton is a continuous variable", {
  
  genset <- sapply(1:100,function(x) paste("Gene",x))
  signature <- sapply(seq(1,100,by=2),function(x) paste("Gene",x))
  SignaturePartition <- GRridge::CreatePartition(signature,varnamesdata=genset)
  
  testthat::expect_length(SignaturePartition,length(SignaturePartition))
  
  })


