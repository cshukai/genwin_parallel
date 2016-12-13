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
    
    for( i in 1:length(recombi_measures)){
    
      fileNameTemp=paste(popgen_metrics[j],recombi_measures[i],sep="_")    
      rhoIdx=which(colnames(d) == recombi_measures[i])
      result=NULL
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
     
        if(!is.null(dom_rm_idx)){
            theseGenWinSize=theseGenWinSize[-dom_rm_idx]
        }
        
        
        rho_na_idx=which(is.na(theseDomRecomBiRate))
        if(!is.null(rho_na_idx)){
            theseGenWinSize=theseGenWinSize[-rho_na_idx]
            theseDomRecomBiRate=theseDomRecomBiRate[-rho_na_idx]
        }
        
        print("---")
        print(length(dom_rm_idx))
        print("---")
        print(length(theseGenWinSize))
        print("---")
        print(length(theseDomRecomBiRate))
        print("---")

        this_pearson_dom=cor(theseGenWinSize,theseDomRecomBiRate) #theseDomRecomBiRate contains NA values now so cor can't be computed
        this_spearman_dom=cor(theseGenWinSize,theseDomRecomBiRate,method="spearman") #theseDomRecomBiRate contains NA values now so cor can't be computed
        this_pearson_dom_test_left=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "less")$p.value
        this_pearson_dom_test_right=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "greater")$p.value
        this_pearson_dom_test_two=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "two.sided")$p.value
        this_spearman_dom_test_left=cor.test(theseGenWinSize,theseDomRecomBiRate,method="spearman",alternative = "less")$p.value
        this_spearman_dom_test_right=cor.test(theseGenWinSize,theseDomRecomBiRate,method="spearman",alternative = "greater")$p.value
        this_spearman_dom_test_two=cor.test(theseGenWinSize,theseDomRecomBiRate,method="spearman",alternative = "two.sided")$p.value
        
        
        
        #this_pearson_dom_pv=this_pearson_dom_test$p.value
        thisRow=c(k,this_pearson_dom,this_pearson_dom_test_left,this_pearson_dom_test_right,this_pearson_dom_test_two,this_spearman_dom,this_spearman_dom_test_left,this_spearman_dom_test_right,this_spearman_dom_test_two)
        result=rbind(result,thisRow)
        colnames(result)=c("chr","pearson","pv_less","pv_greater","pv_two_side","spearman","pv_less","pv_greater","pv_two_side")
        csv_filename=paste(fileNameTemp,"csv",sep=".")
        write.csv(result,csv_filename,row.names=F)

    }
   }    
    
}

csvFiles=Sys.glob("*.csv")

for(i in 1:length(csvFiles)){
    thisCsv=read.csv(csvFiles[i])
    nice.pearson=NULL
    nice.spearman=NULL
    for(j in 1:nrow(thisCsv)){
        this_spearman=thisCsv[j,"spearman"]
        this_pearson=thisCsv[j,"pearson"]
        print(this_spearman)
        print(this_pearson)
        if(this_pearson>0){
            if(thisCsv[j,"pv_greater"]<=0.05){
                 nice.pearson=rbind(nice.pearson,c(j,this_pearson,thisCsv[j,"pv_greater"]))
            }
            
        }
        
        if(this_pearson<0){
            if(thisCsv[j,"pv_less"]<=0.05){
                 nice.pearson=rbind(nice.pearson,c(j,this_pearson,thisCsv[j,"pv_less"]))
            }
            
        }
        
        if(this_spearman>0){
            if(thisCsv[j,"pv_greater.1"]<=0.05){
                 nice.spearman=rbind(nice.spearman,c(j,this_spearman,thisCsv[j,"pv_greater.1"]))
            }
            
        }
        
        if(this_spearman<0){
            if(thisCsv[j,"pv_less.1"]<=0.05){
                 nice.spearman=rbind(nice.spearman,c(j,this_spearman,thisCsv[j,"pv_less.1"]))
            }
            
        }
    
        
        
    }
    
        thisSpearName=paste(csvFiles[i],"spearman.csv",sep=".")
        thisPeasonName=paste(csvFiles[i],"pearson.csv",sep=".")
    
      if(!is.null(nice.spearman)){
           write.csv(nice.spearman,thisSpearName,row.names=F)
        }
        
        if(!is.null(nice.pearson)){
           write.csv(nice.pearson,thisPeasonName,row.names=F)
        }
}


resultFiles=Sys.glob("*csv*.csv") 
bigTable=NULL
for(i in 1:length(resultFiles)){
    metric=unlist(strsplit(x=resultFiles[i],split="_rho"))[1]
    
    cor_type=gsub(x=unlist(strsplit(x=resultFiles[i],split="csv"))[2],pattern="\\.",replacement="")
    
    tmp= unlist(strsplit(x=resultFiles[i],split="rho"))[2]
    rhoSpecies=sub(x=unlist(strsplit(x=tmp,split=".csv"))[1],pattern="_",replacement="")
    
    thisResult=read.csv(resultFiles[i])
   
    for(j in 1:nrow(thisResult)){
        bigTable=rbind(bigTable,cbind(thisResult[j,],metric,cor_type,rhoSpecies))
    }
}

colnames(bigTable)=c("chr","correlation","pv","metrics","correlation_typ","rho_species")
write.csv(bigTable,"bigTable.csv",row.names=F)