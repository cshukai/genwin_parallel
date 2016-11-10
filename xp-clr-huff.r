#################useful fields##########################
#chr: AGPv1 chromosome
#XPCLR_DOM: XP-CLR statistic for selection during domestication (Chen etal. 2010)
#XPCLR_IMP: XP-CLR statistic for selection during improvement (Chen et al. 2010)
#winstart: AGPv1 position for window start site
#winend: AGPv1 position for window end site
#genicbp: number of genic base pairs in window

#####################input & data preprocessing##################
d=read.table("Hufford_et_al._2012_10kb_statistics.txt",sep="\t",header=T)
d_dom=d[which(!is.na(d[,"XPCLR_DOM"])),]
d_imp=d[which(!is.na(d[,"XPCLR_IMP"])),]

######################fit genwin ######################


save.image("xpclr.RData")