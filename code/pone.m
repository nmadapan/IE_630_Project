% A = [D_T;D_O1;D_O2;D_O3;D_O4];
% zero_ids = find(sum(A, 1) == 0);
% nzero_ids = setdiff(1:size(A,2), zero_ids);
% A = A(:,nzero_ids);
% col = size(A,2);
% n = size(A,1);
% nt = size(D_T,1);
% zt = zeros(n,1);
% zt(1:nt) =1;
% A = [A zt];
% zo1 = zeros(n,1);
% no1 = size(D_O1,1);
% zo1(nt+1:nt+no1) = -1;
% 
% zo2 = zeros(n,1);
% no2 = size(D_O2,1);
% zo2(nt+no1+1:nt+no1+no2) = -1;
% 
% zo3 = zeros(n,1);
% no3 = size(D_O3,1);
% zo1(nt+no1+no2+1:nt+no1+no2+no3) = -1;
% 
% zo4 = zeros(n,1);
% no4 = size(D_O4,1);
% zo1(nt+no1+no2+no3+1:nt+no1+no2+no3+no4) = -1;
% 
% A = [A zo1 zo2 zo3 zo4];
% 
% % N = [D_N1;D_N2;D_N3];
% % ns = 
% 
% A = 1e-3 * A;
% % A = A*100;
% 
% k = size(A,2);
% c1 = zeros(1,k);
% c1(col+1) = 1;
% c2 = zeros(1,k);
% c2(col+2) = 1;
% c3 = zeros(1,k);
% c3(col+3) = 1;
% c4 = zeros(1,k);
% c4(col+4) = 1;
% c5 = zeros(1,k);
% c5(col+5) = 1;
% C = [c1;c2;c3;c4;c5];
% 
% bt = ones(nt,1) * 1.05*27*0.001*dV*0.5;
% bo = ones(no1+no2+no3+no4,1)*0.1*20.1*0.001*dV;
% b =[bt;bo];

% bo1 = 0.1*20.1*0.001*(4/3)* pi * (1*0.5*2)
% bo2 = 0.1*20.1*0.001*(4/3) * pi * (0.5*0.5*0.5)
% bo3 = 0.1*20.1*0.001*(4/3) * pi * (0.5*(1/3)*0.25)
% bn1 = 0.2*17.2*0.001*(10*2*2)
% bn3 = 0.2*17.2*0.001*(5*2*2)
% bn2 = 0.2*17.2*0.001*(6*2*2)
% bn = 0.2*17.2*0.001*(pi*20*10*2);
% Dt_N = Dt_N1+Dt_N2 + Dt_N3;
% bt = 1.05*27*0.001*(4/3) * pi * (1)*(1/2)
% 
% 
% 
% b = [bo1;bo2;bo3;bn;bt];
% 
% A = 1e-2*[Dt_O1; Dt_O2; Dt_O3;Dt_N;Dt_T];
% C= [C_O;C_N; -1 * Dt_T];
% nc = size(A,1);
% I = eye(nc);
% I(1:nt,:) =-1*I(1:nt,:);
bo1 = 0.1*20.1*0.001*(4/3)* pi * (1*0.5*2)
bo2 = 0.1*20.1*0.001*(4/3) * pi * (0.5*0.5*0.5)
bo3 = 0.1*20.1*0.001*(4/3) * pi * (0.5*(1/3)*0.25)
bo4 = 0.1*20.1*0.001*(4/3) * pi * (0.5*0.5*0.5)
bn1 = 0.2*17.2*0.001*(10*2*2)
bn3 = 0.2*17.2*0.001*(5*2*2)
bn2 = 0.2*17.2*0.001*(6*2*2)
bn = 0.2*17.2*0.001*(pi*20*10*2);
Dt_N = Dt_N1+Dt_N2 + Dt_N3 +Dt_N4;
bt = 1.05*27*0.001*(4/3) * pi * (1)*(1/2)

lb = zeros(size(C_N)); % Lower bound

b = [bo1;bo2;bo3;bo4;bn; bt];
% C= [C_O;C_N; -1 * Dt_T];
A =  1e-2*[Dt_O1; Dt_O2; Dt_O3;Dt_O4;Dt_N;Dt_T];

% I = eye(6);
% I = -1*I;
% I(end,end) =1;
% A = [A I];
% C = zeros(6,100);
% I = eye(6);
% C = [C I];

I = eye(6);
I(end,end) =-1;
A = [A I];
c=[]

% nc = size(A,1);
% I = eye(nc);
% A = [A I];
% c = zeros(391,1);
% c(220:end)=1;
% lb = zeros(391,1);

% [x,fval] = linprog(c,A,b,[],[],lb)




%lin_solve is my LP solver. It returns the optimum cost and the optimum x value.

[T,B] = lin_solve(A,b,c)


%In lin_solve, each part of the model is handeled by a function.


function [T,B] = lin_solve(A,b,c)
m = size(A,1);
n=size(A,2);

[T,B] = initialtab(A,b); %initialtab() adds the artificial variables and initializes the initial tableau.
                         % T is the tableau and B is an array whic has the indices of the basic variables.

[T,B,f] =ph_one(T,B);   %phase one of the two phase approach. It also checks for feasibility and returns f=1 
                        %when feasible

if f==1
[T,B] = red_remove(T,B,m,n); % checks and removes the redundant constraints if any. It reurns the tableau without
                             % the artificial variables


end
end








function [T,B] = initialtab(A,b) %adds the artificial variables and initializes the initial tableau.
m = size(A,1);
n=size(A,2);
i = eye(m); %artificial variables
A = [A,i];
c_b = ones(m,1);
c = [zeros(n,1);c_b];

top = (c_b)'*A - c';
f = [top;A];
cost = (c_b)'*b;
r =[cost;b];
T =[f,r];
B = n+1 : n+m;
end



function [T,B,f]= ph_one(T,B)
o =0;
while o==0
    o=chkopt(T);
    if (o==0)
        index = findpivot(T,B);
        
        [T,B] = pivot(T,index,B);
        
    else 
        if(T(1,end)<1e-6)  %similar to T(1,end)==0, but only till 10^-6.
            disp("LP is Feasible"); %check for feasibility, cost=0 at the end of phase 1.
            f=1;
        else
            disp("LP is not Feasible");
            f=0;
        end
        
        
    end
end
end
    


function [T,B] = red_remove(T,B,m,n)
%if write functions for removing red. constraints
art = n+1:m+n;
k = ismember(B,art);
if k==0
    disp("No redundant Constraints");
else
    disp("Redundant Constraints found");
    
    l = find(k==1);
    for dum=1:length(l)
        p = ismember(B,art);
        i = find(p==1,1);
        
        if T(i+1,1:n)==0
            T(i+1,:)=[];
            B(i)=[];
        else
            Y=T(i+1,:);
            j = find(Y~=0,1);
            index = [i+1,j];
            [T,B] = pivot(T,index,B);
        end
    end
    disp("Redundant Constraints Removed");
    disp(T);disp(B);
end

T(:,n+1:end-1)=[];
end


function [T,B] = phtwoinitialize(c,T,B) %phase 2 initialization - reduced costs
T(1,1:end-1) = -1*c';
T(1,end) = 0;
for j = B
    i =find(j==B(:));
    x = T(1,j)/T(i+1,j);
    T(1,:) = T(1,:)-x*T(i+1,:);
end
end

    

function [T,B,u] = ph_two(T,B) % phase 2
o = 0;
while o==0
    o=chkopt(T);               %chkopt() checks if the tableau is optimal and returns 1 when optimal.
    if o==0
        [index,u] = findpivot_ptwo(T,B); %findpivot_ptwo() returns the index of the pivot element
        if u==1                          %It also checks if the problem is unbounded, if unbounded, 
            disp("Unbounded Problem");   %it returns u=1.
            return;
        end
        [T,B] = pivot(T,index,B);
    else
        disp("Optimal Solution found.") %when tableau is optimal
        
    end
end
end



function [o_c,x]= disp_sol(T,B) % displays the solution
o_c = T(1,end);
n = size(T,2)-1;
x = zeros(n,1);
for k = B
    i = find(k==B(:));
    x(k) = T(i+1,end);
end
disp("Optimal Cost:");  %optimal cost
disp(o_c);
disp("solution:");      %optimal x value
disp(x);
end




function [T,B] = pivot(T,index,B) %does the pivot operation on the tableau.
i = index(1);
j = index(2);
m = size(T,1);
B(i-1) = j;  %because there is an extra reduced cost row in the tableau
a = T(i,j);
T(i,:) = T(i,:)/a;
for l = 1:m
    if l~=i
        x = T(l,j);
        T(l,:) = T(l,:) - x*T(i,:);
    end
end

end




%Finding pivots index using bland's rule. B is the xi's in the basis.
function [index] = findpivot(T,B) %finds the pivot element based on Bland's rule. 
A = T(1,1:end-1);

j = find(A>0,1);
disp(A);

R = T(2:end,end)./T(2:end,j);
R(R<=0) = inf;
k = find(R==min(R(:)));
if length(k)>1
    i = find(B==min(B(k)));

else
    i=k;
end

index = [i+1,j];
end
    


    
function [index,u] = findpivot_ptwo(T,B) %finds the pivot element based on Bland's rule.
	                                     %also checks for unboundedness.
A = T(1,1:end-1);

j = find(A>0,1);

W=T(2:end,j);
if W<=0                      %checking for unboundedness. expression evaluates the whole vector.
    u =1;index=[1,1];return; % u, to check if unbounded.
end

R = T(2:end,end)./T(2:end,j);
R(R<=0) = inf;
k = find(R==min(R(:)));
if length(k)>1
    i = find(B==min(B(k)));

else
    i=k;
end
u =0;
index = [i+1,j];
end
    

     

function o = chkopt(T) %checks if the tableau is optimal or not. 
c = T(1,1:end-1);     
if all(c(:)<=0)
    o =1;              %returns 1 if optimal
else
    o =0;              %returns 0 if not optimal
end
end






