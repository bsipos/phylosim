library("R.oo");

dest.path<-"./pkg/man";
sf       <-list.files()
Rdoc$compile(filename="PhyloSimSource.R", destPath=dest.path,verbose=TRUE, source=TRUE);
warnings();


