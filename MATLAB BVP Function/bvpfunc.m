function dydz = bvpfunc(x,y) % equation to solve
global v lbg k kbg d
%dydz = zeros(7,1);
A1=y(1);A2=y(2);Rb1=y(3);Rb2=y(4);Rd1=y(5);Rd2=y(6);I=y(7);
q = A1/Rb1;
dydz = [A2
       (v*A2 - p(I,q)*A1 + lbg*A1) / d
       Rb2
       (v*Rb2 - ro(q,Rd1)*A1 + lbg*Rb1) / d
       Rd2
       (ro(q,Rd1)*A1 - lbg*Rb1)/ d
       (k*A1 + kbg)*I];
end