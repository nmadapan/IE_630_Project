
%max Cx
% Ax <= b
%b column vector
% C = [6,4,5;0,0,1];
% A= [1,1,2;1,2,1;2,1,1];
% b=[12;12;12];

% C = [1,1,0;1,0,-2;-1,0,1];
% A=[1,1,0;0,1,0;1,-1,1];
% b =[1;2;4];
C= -1*[C_O;C_N; -1 * Dt_T];

T = T(2:end,:);

rw = size(T,1);
% col = size(A,2);
% I = eye(rw);
% T = [A I];C;
C = -1*C;
k = size(C,1);
z = zeros(k,rw);
C = [C z];
c_s = sum(C);
C = [c_s;C];
 z = zeros(k+1,1);
C = [C z];
% T = [T b];
T = [C; T];
% B = [rw+1:(rw+col)];
n = size(T,2)-1;



op=1;
J = [0];
V=[];
U=[];
itr = 0;
while op >0
    
    if size(V,1)==4
        break;
    end

    
    c = eff_chk(T,B,k,n,1);
    if c ==0
        while c==0
            [index,u] = lp_pid(T,B,k);
            [T,B] = pivot(T,B,index,k,1);
            c = eff_chk(T,B,k,n,1);
        end
    end
    if c==1
        tmp = zeros(1,n);
        tmp(B) = 1;
        if length(J)==1
            J = tmp;
            V = sol_up(T,B,k,n,1,V);
        else
            disp('Chk')
            disp(size(tmp));
            disp(size(J));
            if ~ismember(tmp,J,'rows')
                J = [J;tmp];
                V = sol_up(T,B,k,n,1,V);
            end
        end
    end
    [index,U_i] = eff_pid(T,B,k,n,1,J)
    if ~isempty(U_i)
        if isempty(U)
            
            U = [U;U_i];
        else
        uc = ~ismember(U_i,U,'rows');
        ad = U_i(uc,:);
        
        U = [U;ad];
        end
    end
    if ~isempty(index)
        [T,B] = pivot(T,B,index,k,1);
        
    else
        if ~isempty(U)
            disp(U);
            break;
        else
            op=0;
        end
    end
end
disp(V)
        

function c = eff_chk(T,B,k,n,s_k)
%n = size(T,2) -1; add this into final code - n # of total variables
a = [1:n];
col = setdiff(a,B);

if s_k
    A = T(2:k+1,col);
    b = T(1,:);
    c = any(all(A<0,2));
    if ~c
        c = ~any(b>0);
    end
else
    A = T(1:k,col);
    c = any(all(A<0,2));
end
end

function [index,u] = lp_pid(T,B,k)
index = [0,0];
c = T(1,1:end-1);
a = find(c == max(c));
j = a(1);


W=T(k+2:end,j);
if W<=0                      %checking for unboundedness. expression evaluates the whole vector.
    u =1;index=[1,1];return; % u, to check if unbounded.
end

R = T(k+2:end,end)./T(k+2:end,j);
R(R<=0) = inf;

k = find(R==min(R(:)));

i = k(1);
u =0;
index = [i,j];
end

function res = x_greater_y(x, y, x_greater)
%%%%%%%%%%%%%%%%%%%%%%%
% Verifies if t * x </> y
%
% Inputs:
% x: 1D array of real numbers
% y: 1D array of real numbers
% x_greater: boolean variable. If true, it evaluates x>y, else x<y.
%   Default value is set to true. 
%
% Return:
%   true or false. 
%%%%%%%%%%%%%%%%%%%%%%%
    if(nargin == 2)
        x_greater = true;
    end
    F = @(m, a, b) (a*m'-b);
    N = numel(x);

    % Find all the points to check
    t = sort([y(:) ./ x(:); 0]); % 0 is added to the list
    nt = 0.5 * ( t(1:end-1) + t(2:end) );
    values = unique([t; nt; max(t)+1]); % 1 + max value is also added. 
    values = values(values > 0);

    res = F(values(:), x(:), y(:));
    if(~x_greater)
        res(res > 0) = inf;
        res(res < 0) = 1;
        res = any(and(sum(res) > 0, sum(res) <= N));
    else
        res(res > 0) = 1;
        res(res < 0) = inf;
        res = any(and(sum(res) > 0, sum(res) <= N));    
    end
end

function [index,U_i] = eff_pid(T,B,k,n,s_k,J)
index = [];
U_i = [];
a = [1:n];
N = setdiff(a,B);
B_c = B;
disp('B:');
disp(B);
disp(T);
i =1;
if s_k    
    C = T(2:k+1,N);
    i = k+2;
else
    C = T(1:k,N);
    i = k+1;
end

b = sign(C);
b1 = b>0;
b1 = sum(b1);
b2 = b<0;
b2=sum(b2);
c = b1 & b2;
C_m = C(:,c);

N_m = N(c);
R = T(i:end,end);
A_m = T(i:end,N);
A_m = A_m(:,c);

if length(N_m)==1
    
    h = R./A_m;
    h(h<=0)=inf;
    t = min(h);
    if t>0
        id = find(h ==t,1)
        
        B_c(id) = N_m(1)
        tmp = zeros(1,n);
        tmp(B_c) = 1
        if ~ismember(tmp,J,'rows')
            index =[id,N_m(1)]
        end
         
    end
else
R = ones(length(B),size(C_m,2)).*R;

h = R ./A_m;

h(h<=0)=inf;
[t,ids] = min(h)
disp(C_m);
N_e = N_m;
ids_e = ids;

for i = 1:length(N_m)
    for chk = i+1:length(N_m)
        if x_greater_y(C_m(:,i)', C_m(:,chk)', true) 
            N_e(chk) =0;
            ids_e(chk)=0;

        elseif x_greater_y(C_m(:,i)', C_m(:,chk)', false) 
            N_e(i) =0;
            ids_e(i)=0;
        end
    end
end
c = (N_e ~=0);
N_e = N_e(c);
ids_e = ids_e(c);

for i=1:length(N_e)
    n_b = B;
    n_b(ids_e(i))=N_e(i);
    tmp = zeros(1,n);
    tmp(n_b) = 1;
    if ~ismember(tmp,J,'rows')
        U_i = [U_i; tmp];
        index(1) = ids_e(i);
        index(2) = N_e(i);
    end
end

if ~isempty(U_i)
disp('ui')
disp(U_i);
U_i(end,:) = [];
end
end
end


function V = sol_up(T,B,k,n,s_k,V)
if s_k
    sol = T(k+2:end,end);
else
    sol = T(k+1:end,end);
end
tmp = zeros(1,n);
tmp(B) = sol;
V = [V;tmp];
end

function [T,B] = pivot(T,B,index,k,s_f)
i = index(1);
j = index(2);
B(i)=j;
i = i+k;
m = size(T,1);
if s_f
    i=i+1;
end
a = T(i,j);
s = eye(m);
c = T(:,j);
s(:,i)= -c(:)/a;
s(i,i) = 1/a;
T = s*T;
end
        