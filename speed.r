##########################test 1#######################
tic()
chr6Spline <- splineAnalyze(Y=chr6$Fst,map=chr6$Position,smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
toc()
#170.4 sec elapsed  8g laptop
#49.945 sec elapsed 32g server datascience-2
# 33.042 sec elapsed 46g  server lindberg


######################test2################################
tic()
chr6Spline <- splineAnalyze2(Y=chr6$Fst,map=chr6$Position,smoothness=100,plotRaw=TRUE,plotWindows=TRUE,method=4)
toc()
#stopCluster(cl)
##