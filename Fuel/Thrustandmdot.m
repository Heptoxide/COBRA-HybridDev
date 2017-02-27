disp('Welcome to the thrust and mass flow rate calculator.');
disp('If there are any bugs, please email yale6459@colorado.edu and report it.');
disp('No freedom units allowed whilst using this program');
VE = input('please input the exit velocity:');
while VE <= 0 
    disp('exit velocity cannot be less than or equal to zero');
    disp('Please use common sense');
    VE = input('please input the exit velocity:');
end
PE = input('please input the exit pressure:');
while PE<0 
    disp('pressure cannot be negative');
    disp('use your brain');
    PE = input('please input the exit pressure:');
end
PA = input('Please input ambient pressure:');
while PA<0 
    disp('pressure cannot be negative');
    disp('use your brain');
    PA = input('Please input ambient pressure:');
end
AE = input('please input the exit area:');
while AE<=0
    disp('area cannot be less than or equal to zero');
    disp('stop wasting my time');
    AE = input('please input the exit area:');
end
N = input('please input the number of ports in the fuel grain:');
N = ceil(N);
while N<1 
        disp('you cannot have a fuel grain with no or negative ports');
        N = input('please input the number of ports in the fuel grain:');
        N = ceil(N);
end
L = input('please input the length of the fuel grain:');
while L<=0 
    disp('you cannot have a fuel grain with zero or negative length');
    disp('stop wasting my time');
    L = input('please input the length of the fuel grain:');
end
mox = input('please input the mass flow rate of oxidizer flows:');
while mox <= 0
    disp('oxider flow rate cannot be less than or equal to zero');
    mox = input('please input the mass flow rate of oxidizer flow:');
end
rho = input('please input the fuel grain density:');
while rho <= 0 
    disp('fuel density cannot be less than or equal to zero');
    rho = input('please input the fuel grain density:');
end
Ri = input('please input the initial combustion prot radius:');
while Ri <=0 
    disp('initial combustion port radius cannot be less than or equal to zero')
    Ri = input('please input the initial combustion prot radius:');
end
t = input('please input the burn time:');
while t <=0 
    disp('burn time cannot be less than or equal to zero');
    t = input('please input the burn time:');
end
a = input('please input regression rate coefficient "a" according to the correct unit');
while a<=0 
    disp('a cannot be less than or equal to zero');
    a = input('please input regression rate coefficient "a" according to the correct unit');
end
n = input ('please input regression rate exponent "n" according to the correct unit');
while n<=0
    disp('n cannot be less than or equal to 0');
    n = input ('please input regression rate exponent "n" according to the correct unit');
end
t2 = linspace(0,t,10000);
mf = massflowrate(VE,PE,PA,AE,N,L,mox,rho,Ri,t,n,a);
F = thrust(VE,mf,PE,PA,AE,a,n,mox);
disp('mass flow rate of fuel is:');
disp(mf);
disp('thrust at time = t is:');
disp(F);
mf2 = zeros(1,numel(t2));
F2 = zeros(1,numel(t2));
for i = numel(t2)
    mf2(i) = massflowrate(VE,PE,PA,AE,N,L,mox,rho,Ri,t2(i),n,a);
    F2(i) = thrust(VE,mf2(i),PE,PA,AE,a,n,mox);
end
plot(t2,mf2);
hold on;
plot(t2,F2);