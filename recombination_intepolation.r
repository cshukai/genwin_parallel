rec=read.table("GeneticMapFromNam.txt",header=T)
fst=read.table("FstTable.txt",header=T)


#clean out missing data and correct misplaced start/end  
colnames(rec)=c("p.rs.","chr","start","end","cM")
rec2=NULL
for(i in 1:nrow(rec)){
  if(!is.na(rec[i,"start"]) && !is.na(rec[i,"end"]) ){
     if(rec[i,"start"]<rec[i,"end"]){
       thisRow=c(rec[i,"p.rs."],rec[i,"chr"],rec[i,"start"],rec[i,"end"],rec[i,"cM"])
       rec2=rbind(rec2,thisRow)
     }
     else{
         thisRow=c(rec[i,"p.rs."],rec[i,"chr"],rec[i,"end"],rec[i,"start"],rec[i,"cM"])
         rec2=rbind(rec2,thisRow)
     }
  }
}
#transform accumulative centimorgan so that all centimorgan values are positive
colnames(rec2)=c("p.rs.","chr","start","end","cM")

chrNum=10
for(i in 1:length(chrNum)){
  these_cm=rec[which(rec[,"chr"]==i),"cM"]
  if(min(these_cm)<1){
     mag=1-min(these_cm)
     for(j in 1:nrow(rec2)){
        if(rec2[j,"chr"]==i){
           rec2[j,"cM"]=rec2[j,"cM"]+mag
        }
     }
  }
}