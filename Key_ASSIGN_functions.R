library(devtools)

gatherFile<-function(baseDir){
  ###gathers all the prediction from baseDir/*/*/pathway_activity_testset* format
  setwd(baseDir)
  filenames<-system("ls */*/pathway_activity_testset*", intern=TRUE)
  filenames
  data=NULL
  for(i in 1:length(filenames)){
    ###reading in the filess one at a time
    f<-read.csv(filenames[i], header=1,row.names=1)
    colnames(f)<-paste(filenames[i],colnames(f),sep='/')
    if(i==1){
      data<-f
    }
    else{
      data<-cbind(data,f)
    }
  }
  return(data)
}

resistant_or_sensitive=function(value,cutoff){
  if(value>cutoff){
    status='S'
  }
  else{
    status='R'
  }
  return (status)
}

short_to_long_TCGA_id=function(longnames=NULL,shortnames=NULL){
  counter=0
  for (j in 1:length(shortnames)){
    if(!is.na(pmatch(shortnames[j],longnames))){
      shortnames[j]<-longnames[pmatch(shortnames[j],longnames, duplicates.ok=F)]  
      counter=counter+1
    }
  }
  print(paste(counter,"names have been changed",sep= " "))
  return(shortnames)
}

merge_drop<-function(x,y,by=0,...){
  new_m<-merge(x,y,by=by,...)
  rownames(new_m)<-new_m$Row.names
  return(new_m[,2:length(colnames(new_m))])
}

pcaplot<-function(mat,sub,center=T,scale=T){
  if(sum(sub)!=length(mat)){
    print("verify the subscripts...exiting now")
  }
  else{
    pca_mat <- prcomp(t(mat), center=center,scale=scale)
    plot(pca_mat)
    plot(pca_mat$x[,1],pca_mat$x[,2])
    abline(h=0,v=0)
    for(i in length(sub):1){
      print(i)
      if(i!=1){
        points(pca_mat$x[sum(sub[1:i-1]):sum(sub[1:i])],pca_mat$x[sum(sub[1:i-1]):sum(sub[1:i]),2],col=i,pch=i)
      }
      else{
        points(pca_mat$x[1:sub[i]],pca_mat$x[1:sub[i],2],col=i,pch=i)
      }
    }
  }
}

assign_easy_multi<-function(trainingData=train, testData=test, trainingLabel1=NULL,g=100,out_dir_base="~/Desktop/tmp",cov=0, geneList=NULL,single=0, sigma_sZero = 0.01, sigma_sNonZero = 1, iter=100000, burn_in=50000){
  if(cov==0 & single==0){
    adapB_folder<-paste(out_dir_base,paste( "adapB_multi",sigma_sZero,sigma_sNonZero,sep='_'),sep='/')
    dir.create(file.path(out_dir_base,paste( "adapB_multi",sigma_sZero,sigma_sNonZero,sep='_')))
    adap_adap_folder<-paste(out_dir_base,paste( "adap_adap_multi",sigma_sZero,sigma_sNonZero,sep='_'),sep='/')
    dir.create(file.path(out_dir_base,paste( "adap_adap_multi",sigma_sZero,sigma_sNonZero,sep='_')))
  }
  else if (cov==0 & single==1){
    adapB_folder<-paste(out_dir_base,paste( "adapB_single",sigma_sZero,sigma_sNonZero,sep='_'),sep='/')
    dir.create(file.path(out_dir_base,paste( "adapB_single",sigma_sZero,sigma_sNonZero,sep='_')))
    adap_adap_folder<-paste(out_dir_base,paste( "adap_adap_single",sigma_sZero,sigma_sNonZero,sep='_'),sep='/')
    dir.create(file.path(out_dir_base,paste( "adap_adap_single",sigma_sZero,sigma_sNonZero,sep='_')))
  }

  if(is.null(geneList)){
    set.seed(1234)
    assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1,
                   geneList=NULL, n_sigGene=g, adaptive_B=T, adaptive_S=F, mixture_beta=F,
                   outputDir=adapB_folder, sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero, iter=iter, burn_in=burn_in)

    set.seed(1234)
    assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1,
                   geneList=NULL, n_sigGene=g, adaptive_B=T, adaptive_S=T, mixture_beta=F,
                   outputDir=adap_adap_folder, sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero, iter=iter, burn_in=burn_in)
  }
  else{
    set.seed(1234)
    assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1,
                   geneList=geneList, n_sigGene=g, adaptive_B=T, adaptive_S=F, mixture_beta=F,
                   outputDir=adapB_folder, sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero, iter=iter, burn_in=burn_in)

    set.seed(1234)
    assign.wrapper(trainingData=trainingData, testData=testData, trainingLabel=trainingLabel1,
                   geneList=geneList, n_sigGene=g, adaptive_B=T, adaptive_S=T, mixture_beta=F,
                   outputDir=adap_adap_folder, sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero, iter=iter, burn_in=burn_in)
  }
}

testSig <- function(sigProtein, numGenes=NA, geneList =NULL, trainingData, testData, trainingLabels, sigma_sZero = 0.01, sigma_sNonZero = 1){
  names(sigProtein)=names(geneList)=strsplit(sigProtein,"_")[[1]][1]
  trainingLabel<-list(control=list(sigProtein=1:trainingLabels[1]),sigProtein=(trainingLabels[1]+1):(trainingLabels[1]+trainingLabels[2]))
  names(trainingLabel$control)=names(trainingLabel)[2]=names(sigProtein)
  if(is.na(numGenes)){
    sub_dir<-paste(basedir,paste(sigProtein,"gene_list", sep="_"),sep='/')
  }
  else{
    sub_dir<-paste(basedir,paste(sigProtein,numGenes,"gene_list", sep="_"),sep='/')
  }
  dir.create(sub_dir)
  assign_easy_multi(trainingData = trainingData,test=testData,trainingLabel1 = trainingLabel,
                    g=numGenes,geneList = geneList,out_dir_base = sub_dir,single = 1,
                    sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero)
}

getGeneList = function(rDataPath){
  load(rDataPath)
  #for a gene list
  output.data$processed.data$diffGeneList
  #signature matrix with coefficients
  output.data$processed.data$S_matrix
}

writeFile = function(variable, filename){
  write.table(variable, filename ,sep='\t', col.names = NA,quote=F)
}

testSig_multi <- function(sigProteins, numGenes=NA,geneList =NULL, trainingData, testData, trainingLabel){
  sub_dir<-paste(basedir,paste(sigProteins,"gene_list", sep="_"),sep='/')
  dir.create(sub_dir)
  assign_easy_multi(trainingData = trainingData,test=testData,trainingLabel1 = trainingLabel,
                    g=numGenes,geneList = geneList,out_dir_base = sub_dir)
}
