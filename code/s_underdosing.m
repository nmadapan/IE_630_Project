bo1 = 0.1*20.1*0.001*(4/3)* pi * (1*0.5*2)
bo2 = 0.1*20.1*0.001*(4/3) * pi * (0.5*0.5*0.5)
bo3 = 0.1*20.1*0.001*(4/3) * pi * (0.5*(1/3)*0.25)
bn1 = 0.2*17.2*0.001*(10*2*2)
bn3 = 0.2*17.2*0.001*(5*2*2)
bn2 = 0.2*17.2*0.001*(6*2*2)
bn = 0.2*17.2*0.001*(pi*20*10*2);
Dt_N = Dt_N1+Dt_N2 + Dt_N3;
bt = 1.05*27*0.001*(4/3) * pi * (1)*(1/2)

lb = zeros(size(C_N)); % Lower bound

b = [bo1;bo2;bo3;bn;-1*bt];
% C= [C_O;C_N; -1 * Dt_T];
A = 1e-2*[Dt_O1; Dt_O2; Dt_O3;Dt_N;-1 * Dt_T];

I = eye(5);
I = -1*I;
A = [A I];
C = zeros(5,75);
I = eye(5);
C = [C I];

fitnessfcn = @(x)[x*C'];
options = optimoptions('gamultiobj','UseVectorized',true);
x = gamultiobj(fitnessfcn,80,A,b,[],[],lb,[],options);