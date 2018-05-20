sysuse auto, clear
describe, full
saveold myauto.dta, replace
rsource, terminator(END_OF_R)
.libPaths()
.libPaths("c:/Users/kkranker/documents/R/myRlib" )
.libPaths()
 
x<-2
x
q();
END_OF_R
