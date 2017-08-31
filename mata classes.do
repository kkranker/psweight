clear all
sysuse auto

mata:

class coord {
  real scalar x, y
  real scalar length(), angle()
}

real scalar coord::length()
{
  return(sqrt(x^2 + y^2))
}

real scalar coord::angle()
{
  return(atan2(y, x)*360/(2*pi()))
}

a=coord()
a.x = a.y = 1
a.angle(); a.length()

class gmatch {
  real matrix X0, X1, y0, y1, w0, w1
  real rowvector diff()
  void new()
}

void gmatch::new()
{
  w0 = w1 = 1
}

real rowvector gmatch::diff()
{
  return(mean(X0) :- mean(X1))
}

mydta = gmatch()
mydta.X0 = (1\2)
mydta.X0
mydta.X1 = (4\3)
mydta.diff()

st_view(mydta.X0, ., "price length","foreign")

// does view still work?
mydta.X0[1..4,.]
stata("replace price=2 in 1")
mydta.X0[1..4,.]


mata describe
mydta.X1= st_data(., "price length","foreign")
mydta.X1[1..4,.]
mata describe

end
