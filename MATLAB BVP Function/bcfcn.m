function res = bcfcn(ya,yb)
global I0 v d
A1_0=ya(1);A2_0=ya(2);Rb1_0=ya(3);Rb2_0=ya(4);Rd2_0=ya(6);I_0=ya(7);
A2_m=yb(2);Rb1_m=yb(3);Rb2_m=yb(4);Rd2_m=yb(6);
% vector of all things that should equal zero?
% first half for z=0, second for z=zmax
% possibly don't need a BC for each value??
res = [I_0 - I0
    v*A1_0 - d*A2_0
    v*Rb1_0 - d*Rb2_0
    Rd2_0
    % z=zmax conditions - unsure of these
    A2_m
    Rb2_m
    d*Rd2_m - v*Rb1_m];
end