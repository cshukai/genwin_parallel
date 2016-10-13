tic()
chr6Spline <- splineAnalyze(Y=chr6$Fst,map=chr6$Position,smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
toc()
#170.4 sec elapsed