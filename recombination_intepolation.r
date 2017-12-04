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
