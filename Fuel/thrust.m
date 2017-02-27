function F = thrust(VE,mf,PE,PA,AE,a,n,mox)
mtot = mf+mox;
P = PE - PA;
F = mtot*VE + P*AE;
end