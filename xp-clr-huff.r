
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
    this_spearman_dom=cor(theseGenWinSize,theseDomRecomBiRate,method="spearman")
    
    
    
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
    
    thisRow=c(i,this_pearson_dom,this_spearman_dom,this_pearson_Imp,this_spearman_Imp)
    result=rbind(result,thisRow)
} 
 
 colnames(result)=c("chr","DOM_Pearson","DOM_Spearman","Imp_Pearson","IMP_Spearman")
save.image("xpclr.RData")
write.csv(result,"xpclr_recombination.csv",row.names=F)