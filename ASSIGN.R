##Combat
library(sva)

#Read in input data
nonzerolog <-read.table("merged.rpkmlog", header=TRUE, row.names=1)

#Read in pheno.txt file that contains phenotype information, example available in example_files/pheno.txt
pheno <- read.table("pheno.txt", header=TRUE, row.names=1)

#Create a model matrix for ComBat, just an intercept term because we are just adjusting for batch
modcombat = model.matrix(~1, data=pheno)

#Create an array of integers to indicate batches for the samples (modify the batch numbers to match your samples)
batch <- c(rep(1,53),rep(2,192))

#Perform ComBat
combat_rpkmlog <- ComBat(dat=nonzerolog, batch=batch, mod=modcombat, numCovs=NULL)

#ASSIGN
library(ASSIGN)
#Create an output directory
dir.create("output")
tempdir <- "output"

#Create a subset of the matrix containing only the training data (modify to match your samples)
sub_train <- combat_rpkmlog[, 1:53]

#Create a list of labels for the training data in the way ASSIGN expects
trainingLabel1 <- list(control = list(sig1=1:9,sig2=1:9), sig1=10:15,sig2=16:21)

#Create a subset of the matrix containing only the test data (modify to match your samples)
sub_test_sample <- combat_rpkmlog[,54:245]

#Create labels for the test data (modify to match your samples)
testLabel1 <- c(rep('tumor1',53),rep('tumor2',192))

#ASSIGN - Preprocess
processed.data <- assign.preprocess(trainingData=sub_train,
                                            testData=sub_test_sample, 
                                            trainingLabel=trainingLabel1, 
                                            geneList=NULL, n_sigGene=rep(200,5))
#ASSIGN - MCMC, you may want to try adaptive_S as well
mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub,
                                  Bg = processed.data$B_vector,
                                  X=processed.data$S_matrix,
                                  Delta_prior_p = processed.data$Pi_matrix,
                                  iter = 50000, adaptive_B=TRUE,
                                  adaptive_S=FALSE, mixture_beta=TRUE)
#ASSIGN - Convergence
trace.plot <- assign.convergence(test=mcmc.chain, burn_in=10000, iter=50000,
                                         parameter="B", whichGene=1,
                                         whichSample=NA, whichPath=NA)
#ASSIGN - Summary
mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=10000,
                                        iter=50000, adaptive_B=TRUE,
                                        adaptive_S=FALSE,mixture_beta=TRUE)
#ASSIGN - Output
assign.output(processed.data=processed.data,
              mcmc.pos.mean.testData=mcmc.pos.mean,
              trainingData=sub_train, testData=sub_test_sample,
              trainingLabel=trainingLabel1,
              testLabel=testLabel1, geneList=NULL,
              adaptive_B=TRUE, adaptive_S=FALSE,
              mixture_beta=TRUE, outputDir=tempdir)