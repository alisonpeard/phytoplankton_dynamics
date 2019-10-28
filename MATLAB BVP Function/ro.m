function y = ro(q,Rd)
global romax qmax qmin m
y = romax*(qmax-q)/(qmax-qmin)*(Rd/(m+Rd));
