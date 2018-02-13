version 15.1
clear mata
cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\"
do gmatchclass.mata
//findfile gmatchclass.mata
//do `"`r(fn)'"'
lmbuild lgmatch, replace dir(PERSONAL)
mata: mata describe using lgmatch
