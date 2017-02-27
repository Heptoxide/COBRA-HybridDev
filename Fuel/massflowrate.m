function mf = massflowrate(VE,PE,PA,AE,N,L,mox,rho,Ri,t,n,a)
blah1 = pi*N;
blah2 = mox/blah1;
blah3 = blah2^n;
blah4 = 2*pi*N*rho*L*a*blah3;
blah5 = a*(2*n+1)*blah3*t+Ri^(2*n+1);
blah6 = blah5^(1-2*n)/(1+2*n);
mf = blah4*blah6;
end