 library(GenWin)
 library(pspline)
 data(chr6)
 
 pspline_gcv = smooth.Pspline(this.d.this.chr[,"winstart"]:this.d.this.chr[,"winend"], this.d.this.chr[, 4], norder = 2, method = 3)     
