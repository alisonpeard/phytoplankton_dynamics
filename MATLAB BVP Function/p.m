function y = p(I,q)
global mumax qmin h
y = mumax*(1 - qmin/q)*(I/(h+I));