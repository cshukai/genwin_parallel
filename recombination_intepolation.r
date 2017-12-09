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
for(i in 1:chrNum){
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


#map window to markers with FST values
library(doSNOW)
threads= 5
cl=makeCluster(threads)
registerDoSNOW(cl)


windows=NULL
foreach(i = 1:nrow(rec2))  %dopar% {
  this_chr=rec2[i,"chr"]
  this_cm=rec2[i,"cM"]
  this_target=fst[which(fst[,"Chromosome"]==this_chr),]
  
  win_start=rec2[i,"start"]
  win_end=rec2[i,"end"]
  fst_mapped=NULL
  for(j in 1:nrow(this_target)){
     if(this_target[j,"Position"]>=win_start && this_target[j,"Position"]<win_end){
       windows=rbind(i,this_chr,win_start,win_end,this_cm,this_target[j,"Position"],this_target[,"Fst"])
     }
  }
  

    
  
  
}
stopCluster(cl)

# fit genwin to every window
