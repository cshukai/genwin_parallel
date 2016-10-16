##########################test for single chromosome#######################
tic()
chr6Spline <- splineAnalyze(Y=chr6$Fst,map=chr6$Position,smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
toc()



threads= 12
cl=makeCluster(threads)
registerDoSNOW(cl)

tic()
chr6Spline2 <- splineAnalyze2(Y=chr6$Fst,map=chr6$Position,smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
toc()
stopCluster(cl)


# without optimization
    #170.4 sec elapsed  8g laptop
    #49.945 sec elapsed 32g server datascience-2
    # 33.042 sec elapsed 46g  server lindberg

# for loop optimization : 10.403 sec elapsed


####################test for multiple chromosome##################
#combining idv files
files=Sys.glob(file.path( "*.out"))
pos_col_idx=3
d=list()
for(i in 1:length(files)){
    thisfile=read.table(files[i],sep="\t",skip=1)
    thisfile2=thisfile[order(thisfile[,pos_col_idx]),]#sort table on snp location 
    d[[i]]=thisfile2
}



allChrSplines=list()

# without optimization 
tic()
for(i in 1:length(d)){
    thisfile=d[[i]]
    allChrSplines[[i]]= splineAnalyze(Y=thisfile[,4],map=thisfile[,3],smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
}
toc()
#66.305 sec elapsed on lindberg
# with optimization

threads= 12
cl=makeCluster(threads)
registerDoSNOW(cl)

tic()
allChrSplines2=list()
foreach(k=1:length(length(d)) %dopar% 
chr6Spline2 <- splineAnalyze2(Y=chr6$Fst,map=chr6$Position,smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
toc()
stopCluster(cl)
