

A = [-1*D_T;D_O1;D_O2;D_O3;D_O4];
col = size(A,2);
n = size(A,1);
nt = size(D_T,1);
zt = zeros(n,1);
zt(1:nt) =-1;
A = [A zt];
zo1 = zeros(n,1);
no1 = size(D_O1,1);
zo1(nt+1:no1) = -1;

zo2 = zeros(n,1);
no2 = size(D_O2,1);
zo2(no1+1:no2) = -1;

zo3 = zeros(n,1);
no3 = size(D_O3,1);
zo1(no2+1:no3) = -1;

zo4 = zeros(n,1);
no4 = size(D_O4,1);
zo1(no3+1:no4) = -1;

A = [A zo1 zo2 zo3 zo4];

A = 1e-3 * A;

k = size(A,2);
c1 = zeros(1,k);
c1(col+1) = 1;
c2 = zeros(1,k);
c2(col+2) = 1;
c3 = zeros(1,k);
c3(col+3) = 1;
c4 = zeros(1,k);
c4(col+4) = 1;
c5 = zeros(1,k);
c5(col+5) = 1;
C = [c1;c2;c3;c4;c5];

bt = ones(nt,1) * 1.05*27*0.001*dV*0.5*-1;
bo = ones(no1+no2+no3+no4,1)*0.1*20.1*0.001*dV;
b =[bt;bo];
lb = zeros(1,k);

fitnessfcn = @(x)[x*C'];
options = optimoptions('gamultiobj','UseVectorized',true);
x = gamultiobj(fitnessfcn,k,A,b,[],[],lb,[],options);

