#################useful fields##########################
#chr: AGPv1 chromosome
#XPCLR_DOM: XP-CLR statistic for selection during domestication (Chen etal. 2010)
#XPCLR_IMP: XP-CLR statistic for selection during improvement (Chen et al. 2010)
#winstart: AGPv1 position for window start site
#winend: AGPv1 position for window end site
#genicbp: number of genic base pairs in window

#####################input & data preprocessing##################
#read data
d=read.table("Hufford_et_al._2012_10kb_statistics.txt",sep="\t",header=T)
d_dom=d[which(!is.na(d[,"XPCLR_DOM"])),]
d_imp=d[which(!is.na(d[,"XPCLR_IMP"])),]
#sort for later spline fit
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
domSpline=splineAnalyze(Y=d_dom_sorted$XPCLR_DOM,map=mean(d_dom_sorted$winstart+d_dom_sorted$winend),smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)

save.image("xpclr.RData")