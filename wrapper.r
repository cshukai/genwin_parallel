# combining idv files
files=Sys.glob(file.path( "*.out"))

d=NULL
for(i in 1:length(files)){
    thisfile=read.table(files[i],sep="\t",skip=1)
    d=rbind(d,thisfile)
}
#sort by chromosome
#sort by snp position

