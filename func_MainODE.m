function FC=func_MainODE()
global npi np N K max_t Ct T phi Phivec beta
for np=npi:npi % # of prophages     TO BE GIVEN
    fprintf('Number of Prophage =  %d\n',np);
N=2^np;% Permutations with repitition
if np>1; sims=func_sims(np); else sims=1; end
% Parameters_M4();
d2b=(dec2bin(0:N-1,np));% decimal to binary conversion of each possible permutation of number of prophage
arr_mx_t=max_t;%[10/100*max_t;20/100*max_t;40/100*max_t;60/100*max_t;100/100*max_t];%;350;400;500;600;700];
for tm=1:length(arr_mx_t)
    tmax=arr_mx_t(tm);
    fprintf('Maximum time: %d \n',tmax);
for num_sim=1:2^(np-1) % Simulations/Prophage sequence Number
T=0; Ct=0;% Current Time
seq_prophage=0;
seq_prophage=sims(num_sim,1:find(sims(num_sim,:),1,'last'));%[3 1 1];
    sprintf('Prophage sequence %s \n',int2str(seq_prophage));
if sum(seq_prophage)~=np; error('The sum of the sequence of injected prophage must be equal to the total number of prophage'); end;
num_pphge=length(seq_prophage);% Number of therapies
Arrphi=func_Arrphi(seq_prophage);% PHI Array
%% Initial conditions and ode solution
L0(1:N,1)=0;
C0(1,1)=K;  C0(2:N,1)=0;
V0(1:np,1)=0;
X0=[L0;C0;V0];
t=0;
t_pphge=0;
for k=1:num_pphge
    if k==1
        t_pphge(2*k+1)=t_pphge(2*k-1)+tmax*seq_prophage(k)/np;
        t_pphge(2*k)=t_pphge(2*k-1)+t_pphge(3)/2;%[0 tmax/np 2*tmax/np 3*tmax/np 4*tmax/np tmax];
    else
        t_pphge(2*k+1)=t_pphge(2*k-1)+tmax*seq_prophage(k)/np;%/2;%[0 tmax/np 2*tmax/np 3*tmax/np 4*tmax/np tmax];
        t_pphge(2*k)=t_pphge(2*k-1)+(t_pphge(2*k+1)-t_pphge(2*k-1))/2;%/k;%[0 tmax/np 2*tmax/np 3*tmax/np 4*tmax/np tmax];
    end
end
% t_pphge(2*k+1)=tmax;
j=1;
options=odeset('RelTol',1e-8,'AbsTol',1e-8);
for k=1:length(t_pphge)-1
    if rem(k,2)==1
    fprintf('Insertion of %s Lysogen with %d prophage \n',num2ordinal(j), seq_prophage(j));
        Phivec(1:N,1)=Arrphi(1:N,j);
        j=j+1;
    else
        fprintf('Insertion of Lysogen stopped \n');
        Phivec(1:N,1)=0;
    end
    [t_i,SOL_i]=ode15s('ODE_sol3',[t_pphge(k) t_pphge(k+1)],X0,options);
    for i=1:length(X0); 
        if SOL_i(end,i)<1/K;  r=rand(1);%%%CHECK IT
            if r<10^-7; SOL_i(end,i)=1;   else  SOL_i(end,i)=0;  end;
        end;
    end;
    tA=length(t);tB=length(t_i);
    if k==1; 
        t(tA:tA+tB-1)=t_i;SOL(np,tm,num_sim,tA:tA+tB-1,1:2*N+np)=SOL_i(1:length(tA:tA+tB-1),1:2*N+np); 
    else
        t(tA+1:tA+tB)=t_i;SOL(np,tm,num_sim,tA+1:tA+tB,1:2*N+np)=SOL_i(1:length(tA+1:tA+tB),1:2*N+np);
    end;
    X0=SOL_i(end,:);
    SOL_i=0;
end
l_t(np,tm,num_sim)=length(t);%SOL(num_sim,:,1));
FC(np,tm,num_sim,1:2*N+np)=SOL(np,tm,num_sim,l_t(np,tm,num_sim),1:2*N+np);
Phivec=0;
% l(tm,num_sim)=length(SOL(tm,num_sim,1,:));
t_all(tm,num_sim,1:l_t(np,tm,num_sim))=t;
end
end
end
toc