function yp = LotkaVolterra(t,y)
b = 0.5;
p = 0.04;
r = 0.02;
d = 1;

yp = diag([b - p*y(2), r*y(1) - d])*y;
return