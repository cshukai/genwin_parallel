#####################data preporcessing#####################
#combining idv files
files=Sys.glob(file.path( "*.out"))
pos_col_idx=3
d=NULL
for(i in 1:length(files)){
    thisfile=read.table(files[i],sep="\t",skip=1)
    thisfile2=thisfile[order(thisfile[,pos_col_idx]),]#sort table on snp location 
    d=rbind(d,thisfile2)
}





#########################HPC############################
splineAnalyze.idvChrom<-function(d,hpc_mode){
    
}

splineAnalyze.multiChrom<-function(d,hpc_mode){
    
}
#sort by chromosome
#sort by snp position

