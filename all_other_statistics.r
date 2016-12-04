# preprocessing
d=read.table("Hufford_et_al._2012_10kb_statistics.txt",sep="\t",header=T)
popgen_metrics=c("fst_teo_LR",      "fst_teo_mz",      "fst_LR_mz",      "fst_mex_mz",      "fst_mex_teo", "fst_mex_LR","ThetaWteo",       "ThetaPiteo",     "ThetaHteo","TajDteo",         "Rminteo","ThetaWMZ",        "ThetaPiMZ",       "ThetaHMZ",        "TajDMZ",          "RminMZ")
recombi_measures=c("rho_teo","rho_MZ","rho_LR")

#pairwise computation of correlation
## 1.fit genwin for every population-genetics metrics 
## 2.evaluate correlation between window size and every kind of rho

for(j in 1:length(popgen_metrics)){
    
    colIdx=which(colnames(d) == popgen_metrics[j])
    this_d=d[which(!is.na(d[,colIdx])),] # clean up
    this_d_sorted= this_d[order(this_d$winstart),]  # sort for further spline analysis
    
    #split by chromosome
    chrNum=10 
    d_chr_spec=list()
    for(i in 1:chrNum){
        d_chr_spec[[i]]=this_d_sorted[which(this_d_sorted[,"chr"]==i),]
    }
    
    #fit genwin
    wins=list()
    for(i in 1:chrNum){
      this_spline=splineAnalyze(Y=d_chr_spec[[i]][,colIdx],map=rowMeans(cbind(d_chr_spec[[i]]$winstart,d_chr_spec[[i]]$winend)),smoothness=100,plotRaw=F,plotWindows=F,method=4)
      wins[[i]]=this_spline$windowData[,c("WindowStart","WindowStop")]    
    }
    
    ##correlation
     ##1. compute the midpoint of genwin window and map back to original Huffold window
     ##2. skip if no corresponding windows can be map or no corresponding rho value exisit
     ##3. take the corresponding rho values
     ##4. compute correlations and significance test 
    result=NULL
    for( i in 1:length(recombi_measures)){
      fileNameTemp=paste(popgen_metrics[j],recombi_measures[i],sep="_")    
      rhoIdx=which(colnames(d) == recombi_measures[i])
      for(k in 1:chrNum){
         theseGenWinSize=wins[[k]][,"WindowStop"]-wins[[k]][,"WindowStart"]
         theseDomWinMidPoint=rowMeans(wins[[k]]) 
         theseDomRecomBiRate=NULL
         dom_rm_idx=NULL
         for(m in 1:length(theseDomWinMidPoint)){
           if(length(intersect(which(d_chr_spec[[k]]$winstart <= theseDomWinMidPoint[m]),which(d_chr_spec[[k]]$winend >= theseDomWinMidPoint[m] )))==0){
             dom_rm_idx=c(dom_rm_idx,m)
             next
           }
              thisDomRecomBiRate=d_chr_spec[[k]][intersect(which(d_chr_spec[[k]]$winstart <= theseDomWinMidPoint[m]),which(d_chr_spec[[k]]$winend >= theseDomWinMidPoint[m] )),rhoIdx]
              theseDomRecomBiRate=c(theseDomRecomBiRate,thisDomRecomBiRate)
           
         }
     
     
        theseGenWinSize=theseGenWinSize[-dom_rm_idx]
        print("---")
        print(length(dom_rm_idx))
        print("---")
        print(length(theseGenWinSize))
        print("---")
        print(length(theseDomRecomBiRate))
        print("---")

        this_pearson_dom=cor(theseGenWinSize,theseDomRecomBiRate)
        
        this_pearson_dom_test_left=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "less")
        this_pearson_dom_test_right=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "greater")
        this_pearson_dom_test_two=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "two.sided")
        #this_pearson_dom_pv=this_pearson_dom_test$p.value
        thisRow=c(k,this_pearson_dom,this_pearson_dom_test_left,this_pearson_dom_test_right,this_pearson_dom_test_two)
        result=rbind(result,thisRow)
        colnames(result)=c("chr","correlation","pv_less","pv_greater","pv","pv")
        csv_filename=paste(fileNameTemp,"csv",sep=".")
        write.csv(result,csv_filename,row.names=F)

    }
  }    
    
}