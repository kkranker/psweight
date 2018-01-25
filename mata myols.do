mata:
void myols(string scalar varlist,
           string scalar touse,
           string scalar wgtvar,
           real   scalar addconstant,
         | string scalar residvar)
{
  real colvector b
  real matrix M, CP
  M=.

  st_view(M, ., varlist, touse)

  if (wgtvar=="") {
    CP = quadcross(M, addconstant, M, addconstant)
  }
  else {
    real colvector w
    w=.
    st_view(w, ., wgtvar, touse)
    CP = quadcross(M, addconstant, w, M, addconstant)
  }
  b  = lusolve(CP[|2,2 \ .,.|],CP[|2,1 \ .,1|])  

  // return the b' and N in r(b) and r(N), respectively
  st_rclear()
  st_matrix("r(b)", b')
  st_numscalar("r(N)", rows(M))

  // optionally, export residuals to new variable
  if (args()==5) {
    real vector e; real scalar idx
    if (addconstant & length(b)<2) e = M[|1,1 \ .,1|] :- b[length(b)]
    else if (addconstant)          e = (M[|1,1 \ .,1|] - M[|1,2 \ .,.|]*b[|1 \ (length(b)-1)|]) :- b[length(b)]
    else                           e = (M[|1,1 \ .,1|] - M[|1,2 \ .,.|]*b)
    idx = st_addvar("double", residvar)
    st_global("r(resid)", residvar)
    st_store(., idx, touse, e)
  }
}
end
