
#####################input & data preprocessing##################
#read data
d=read.table("Hufford_et_al._2012_10kb_statistics.txt",sep="\t",header=T)
d_dom=d[intersect(which(!is.na(d[,"XPCLR_DOM"])),which(!is.na(d[,"rho_MZ"]))),]
d_imp=d[intersect(which(!is.na(d[,"XPCLR_IMP"])),which(!is.na(d[,"rho_MZ"]))),]
#sort by position
d_dom_sorted=d_dom[order(d_dom$winstart),]
d_imp_sorted=d_imp[order(d_imp$winstart),]
#split by chromosome
chrNum=10
chr_spec_dom=list()
chr_spec_imp=list()
for(i in 1:chrNum){
    chr_spec_dom[[i]]=d_dom_sorted[which(d_dom_sorted[,"chr"]==i),]
    chr_spec_imp[[i]]=d_imp_sorted[which(d_imp_sorted[,"chr"]==i),]
    
}


######################fit genwin ######################
dom_win=list() # for futher retrieval of window size
imp_win=list() 
for(i in 1:chrNum){
    # taking the midpoint as position
    domSpline=splineAnalyze(Y=chr_spec_dom[[i]]$XPCLR_DOM,map=rowMeans(cbind(chr_spec_dom[[i]]$winstart,chr_spec_dom[[i]]$winend)),smoothness=100,plotRaw=F,plotWindows=F,method=4)
    dom_win[[i]]=domSpline$windowData[,c("WindowStart","WindowStop")]

    impSpline=splineAnalyze(Y=chr_spec_imp[[i]]$XPCLR_IMP,map=rowMeans(cbind(chr_spec_imp[[i]]$winstart,chr_spec_imp[[i]]$winend)),smoothness=100,plotRaw=F,plotWindows=F,method=4)
    imp_win[[i]]=impSpline$windowData[,c("WindowStart","WindowStop")]
    
}

##########################compute correlation between recombination rate and window size
#take the midpoint of spined-window to represent each XP-CLR region and map back to original window in paper

pdf("winsize_vs_rho.pdf")
result=NULL
for(i in 1:chrNum){
   #dom correlation
    theseGenWinSize=dom_win[[i]][,"WindowStop"]-dom_win[[i]][,"WindowStart"]
    theseDomWinMidPoint=rowMeans(dom_win[[i]])
    theseDomRecomBiRate=NULL
    dom_rm_idx=NULL
    for(j in 1:length(theseDomWinMidPoint)){
        if(length(intersect(which(chr_spec_dom[[i]]$winstart <= theseDomWinMidPoint[j]),which(chr_spec_dom[[i]]$winend >= theseDomWinMidPoint[j] )))==0){
           dom_rm_idx=c(dom_rm_idx,j)
           next
        }
    
        thisDomRecomBiRate=chr_spec_dom[[i]][intersect(which(chr_spec_dom[[i]]$winstart <= theseDomWinMidPoint[j]),which(chr_spec_dom[[i]]$winend >= theseDomWinMidPoint[j] )),"rho_MZ"]
        theseDomRecomBiRate=c(theseDomRecomBiRate,thisDomRecomBiRate)
        
    }   
    print(length(dom_rm_idx)/length(theseGenWinSize))
    theseGenWinSize=theseGenWinSize[-dom_rm_idx]
    this_pearson_dom=cor(theseGenWinSize,theseDomRecomBiRate)
    this_pearson_dom_test=cor.test(theseGenWinSize,theseDomRecomBiRate,method="pearson",alternative = "greater")
    this_pearson_dom_pv=this_pearson_dom_test$p.value
    
    this_spearman_dom=cor(theseGenWinSize,theseDomRecomBiRate,method="spearman")
    this_spearman_dom_test=cor.test(theseGenWinSize,theseDomRecomBiRate,method="spearman",alternative = "greater")
    this_spearnman_dom_pv=this_spearman_dom_test$p.value


    this_chr=paste("chr",i,sep=" ")
    plot(x=theseGenWinSize,y=theseDomRecomBiRate,xlab="GenWin Window Size", ylab="Rho",main=this_chr,sub="domestic")
    abline(lm(theseDomRecomBiRate ~ theseGenWinSize),col="red")
    
    #imp correlation
    theseGenWinSize=imp_win[[i]][,"WindowStop"]-imp_win[[i]][,"WindowStart"]
    theseImpWinMidPoint=rowMeans(imp_win[[i]])
    theseImpRecomBiRate=NULL
    imp_rm_idx=NULL
    for(j in 1:length(theseImpWinMidPoint)){
        if(length(intersect(which(chr_spec_imp[[i]]$winstart <= theseImpWinMidPoint[j]),which(chr_spec_imp[[i]]$winend >= theseImpWinMidPoint[j] )))==0){
            imp_rm_idx=c(imp_rm_idx,j)
            next
        }
       
        thisImpRecomBiRate=chr_spec_imp[[i]][intersect(which(chr_spec_imp[[i]]$winstart <= theseImpWinMidPoint[j]),which(chr_spec_imp[[i]]$winend >= theseImpWinMidPoint[j] )),"rho_MZ"]
        theseImpRecomBiRate=c(theseImpRecomBiRate,thisImpRecomBiRate)
        
    }
    print(length(imp_rm_idx)/length(theseGenWinSize))

    theseGenWinSize=theseGenWinSize[-imp_rm_idx]

    this_pearson_Imp=cor(theseGenWinSize,theseImpRecomBiRate)
    this_spearman_Imp=cor(theseGenWinSize,theseImpRecomBiRate,method="spearman")
    
    this_pearson_imp_test=cor.test(theseGenWinSize,theseImpRecomBiRate,method="pearson",alternative = "greater")
    this_pearson_imp_pv=this_pearson_imp_test$p.value
    this_spearman_imp_test=cor.test(theseGenWinSize,theseImpRecomBiRate,method="spearman",alternative = "greater")
    this_spearnman_imp_pv=this_spearman_imp_test$p.value
    
    
    plot(x=theseGenWinSize,y=theseImpRecomBiRate,xlab="GenWin Window Size", ylab="Rho",main=this_chr,sub="improved")
    abline(lm(theseImpRecomBiRate ~ theseGenWinSize),col="red")
    thisRow=c(i,this_pearson_dom,this_pearson_dom_pv,this_spearman_dom, this_spearnman_dom_pv ,this_pearson_Imp,this_pearson_imp_pv,this_spearman_Imp,this_spearnman_imp_pv)
    result=rbind(result,thisRow)
} 
dev.off() 
colnames(result)=c("chr","DOM_Pearson","pv","DOM_Spearman","pv","Imp_Pearson","pv","IMP_Spearman","pv")
save.image("xpclr.RData")
write.csv(result,"xpclr_recombination.csv",row.names=F)


#combine chromosome
allWindowSize=NULL
allDomCombinationRate=NULL
for(i in 1:chrNum){
   #dom correlation
    theseGenWinSize=dom_win[[i]][,"WindowStop"]-dom_win[[i]][,"WindowStart"]
    theseDomWinMidPoint=rowMeans(dom_win[[i]])
    theseDomRecomBiRate=NULL
    dom_rm_idx=NULL
    for(j in 1:length(theseDomWinMidPoint)){
        if(length(intersect(which(chr_spec_dom[[i]]$winstart <= theseDomWinMidPoint[j]),which(chr_spec_dom[[i]]$winend >= theseDomWinMidPoint[j] )))==0){
           dom_rm_idx=c(dom_rm_idx,j)
           next
        }
    
        thisDomRecomBiRate=chr_spec_dom[[i]][intersect(which(chr_spec_dom[[i]]$winstart <= theseDomWinMidPoint[j]),which(chr_spec_dom[[i]]$winend >= theseDomWinMidPoint[j] )),"rho_MZ"]
        theseDomRecomBiRate=c(theseDomRecomBiRate,thisDomRecomBiRate)
        
    }   
    theseGenWinSize=theseGenWinSize[-dom_rm_idx]
    allWindowSize=c(allWindowSize,theseGenWinSize)
    allDomCombinationRate=c(allDomCombinationRate,theseDomRecomBiRate)
}
all_dom_pearson=cor(allDomCombinationRate,allWindowSize)
all_dom_spearman=cor(allDomCombinationRate,allWindowSize,method="spearman")
all_dom_pearson_test=cor.test(allDomCombinationRate,allWindowSize,method="pearson",alternative="greater")
all_dom_spearman_test=cor.test(allDomCombinationRate,allWindowSize,method="spearman",alternative="greater")

allWindowSize=NULL
allImpCombinationRate=NULL
for(i in 1:chrNum){
    #imp correlation
    theseGenWinSize=imp_win[[i]][,"WindowStop"]-imp_win[[i]][,"WindowStart"]
    theseImpWinMidPoint=rowMeans(imp_win[[i]])
    theseImpRecomBiRate=NULL
    imp_rm_idx=NULL
    for(j in 1:length(theseImpWinMidPoint)){
        if(length(intersect(which(chr_spec_imp[[i]]$winstart <= theseImpWinMidPoint[j]),which(chr_spec_imp[[i]]$winend >= theseImpWinMidPoint[j] )))==0){
            imp_rm_idx=c(imp_rm_idx,j)
            next
        }
       
        thisImpRecomBiRate=chr_spec_imp[[i]][intersect(which(chr_spec_imp[[i]]$winstart <= theseImpWinMidPoint[j]),which(chr_spec_imp[[i]]$winend >= theseImpWinMidPoint[j] )),"rho_MZ"]
        theseImpRecomBiRate=c(theseImpRecomBiRate,thisImpRecomBiRate)
        
    }

    theseGenWinSize=theseGenWinSize[-imp_rm_idx]
    allWindowSize=c(allWindowSize,theseGenWinSize)
    allImpCombinationRate=c(allImpCombinationRate,theseImpRecomBiRate)
   
} 


all_imp_pearson=cor(allImpCombinationRate,allWindowSize)
all_imp_spearman=cor(allImpCombinationRate,allWindowSize,method="spearman")
all_imp_pearson_test=cor.test(allImpCombinationRate,allWindowSize,method="pearson",alternative="greater")
all_imp_spearman_test=cor.test(allImpCombinationRate,allWindowSize,method="spearman",alternative="greater")
