"msep" <-
function(n,s,p.tr)
{

expected=0
for(Y in 0:n)
  {
  expected=expected+( (1-(1-Y/n)^(1/s)) * choose(n,Y) * ((1-(1-p.tr)^s)^Y) * ((1-p.tr)^(s*(n-Y)) ) )
  }
expected


varsum=0
for(Y in 0:n)
  {
  varsum=varsum+ ( ( (1-Y/n)^(2/s) ) * choose(n,n-Y) * ( ((1-p.tr)^s)^(n-Y) ) * ((1-(1-p.tr)^s)^Y) ) 
  }
varp=varsum - (1-expected)^2

expp=expected
bias=expected-p.tr
mse=varp + bias^2


list(varp=varp,
  mse=mse,
  bias=bias,
exp=expected)
}

