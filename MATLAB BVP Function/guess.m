function g = guess(z)
global I0 kbg
zmax = z(end);
mu = zmax/2;
sigma2 = zmax/10;

g = [exp(-(z-mu).*(z-mu)/sigma2) %A1
     -2*(z-mu).*exp(-(z-mu).*(z-mu)/sigma2) %A2
     exp(-(z-mu).*(z-mu)/sigma2)./sigma2 %Rb1
     -2*(z-mu).*exp(-(z-mu).*(z-mu))./sigma2 %Rb2
     0.5*z %Rd1
     0.5*zeros(1,length(z)) %Rd2
     I0*exp(-kbg*z)];
end