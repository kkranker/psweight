dir <- c("C:\\Users\\kkranker\\Documents\\Stata\\Ado\\Devel\\gmatch\\")
mydata <- read.csv(paste0(dir, "testfile.csv"), stringsAsFactors = F) 
library(CBPS)

fit_ATE       <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 0, method='exact')
fit_ATET      <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="exact")
fit_ATET_over <- CBPS(treat ~ x1 + x1 + X__000000 + X__000001 +x4 + x5 + x6 +x7 +x90 +x91+ x92 +x93 +x94+ x95, data = mydata, ATT = 1, method="over")
summary(fit_ATE)
summary(fit_ATET)
summary(fit_ATET_over)




# . sysuse auto, clear
# . describe, full
# . saveold myauto.dta, replace
# . rsource, terminator(END_OF_R)
# . library(foreign);
# . rauto<-read.dta("myauto.dta", convert.f=TRUE);
# . rauto;
# . attributes(rauto);
# . q();
# . END_OF_R
