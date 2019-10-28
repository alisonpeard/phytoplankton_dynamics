% Using MATLABS bvp function
% https://www.mathworks.com/help/matlab/ref/bvp4c.html
% sol = bvp4c(odefun,bcfun,solinit)
global v I0 lbg k kbg I0 mumax romax m qmax qmin h d
v=0.25; lbg=0.1; k=0.0003; kbg=0.4; I0=300; mumax=1.2;
romax=0.2;m=15; qmax=0.04;qmin=0.004; h=120; d=10;
%check I0 not redefined anywhere
zmax=30;
zmesh = linspace(0,zmax);

solinit = bvpinit(zmesh, @guess); % returns struct data type

sol = bvp4c(@bvpfunc, @bcfcn, solinit);
plot(sol.x, sol.y, '-o')

%need to change bc func