%% This program is meant to simulate the model related to phage bacteria interaction.
% This program requires input of number of prophage (np), and maximum time array
% (arr_mx_t). The lysogens with different combinations of prophage are set
% in the array (sims) made by calling function func_sims(np). The parameter
% seq_prophage take one row of each "sims" array and inject one by one the 
% lysogen type with number of prophage as mentioned in seq_prophage.

clear all
%MCV rechecked; MLI rechecked; MLV rechecked; MVL rechecked;MVC
tic
global np K Ct T N phi Phivec %MVL MVC L C MLV MCV V
max_t=3000000; % Maximum time  TO BE GIVEN
dt=max_t/5;
n=5;
for np=1:n % # of prophages     TO BE GIVEN
    fprintf('Number of Prophage =  %d\n',np);
N=2^np;% Permutations with repitition
if np>1; sims=func_sims(np); else sims=1; end
Parameters();
d2b=(dec2bin(0:N-1,np));% decimal to binary conversion of each possible permutation of number of prophage
arr_mx_t=[max_t*5/5];%;350;400;500;600;700];
for tm=1:length(arr_mx_t)
    tmax=arr_mx_t(tm);
    fprintf('Maximum time: %d \n',tmax);
for num_seq=2^(np-1):2^(np-1) % Simulations/Prophage sequence Number
T=0; Ct=0;% Current Time
seq_prophage=0;
seq_prophage=sims(num_seq,1:find(sims(num_seq,:),1,'last'));%[3 1 1];
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
% Lysogen insertion and non-insertion period
for k=1:num_pphge
%     if k==1
%         t_pphge(2*k+1)=t_pphge(2*k-1)+tmax*seq_prophage(k)/np;
%         t_pphge(2*k)=t_pphge(2*k-1)+t_pphge(3)/2;%[0 tmax/np 2*tmax/np 3*tmax/np 4*tmax/np tmax];
%     else
        t_pphge(2*k+1)=t_pphge(2*k-1)+tmax*seq_prophage(k)/np;%/2;%[0 tmax/np 2*tmax/np 3*tmax/np 4*tmax/np tmax];
        t_pphge(2*k)=15*np;%t_pphge(2*k-1)+tmax*seq_prophage(k)/np/4;%/k;%[0 tmax/np 2*tmax/np 3*tmax/np 4*tmax/np tmax];
%     end
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
    if t_pphge(k) == t_pphge(k+1)
        t_i=0;SOL_i=0;
    else
    [t_i,SOL_i]=ode15s('ODE_sol3',[t_pphge(k) t_pphge(k+1)],X0,options);
    for i=1:length(X0)
        if SOL_i(end,i)<1/K;  r=rand(1);%%%CHECK IT
            if r<SOL_i(end,i); SOL_i(end,i)=1;   else  SOL_i(end,i)=0;  end;
        end;
    end;
    tA=length(t);tB=length(t_i);
    if k==1
        t(tA:tA+tB-1)=t_i;SOL(np,tm,num_seq,tA:tA+tB-1,1:2*N+np)=SOL_i(1:length(tA:tA+tB-1),1:2*N+np); 
    else
        t(tA+1:tA+tB)=t_i;SOL(np,tm,num_seq,tA+1:tA+tB,1:2*N+np)=SOL_i(1:length(tA+1:tA+tB),1:2*N+np);
    end;
    X0=SOL_i(end,:);
    SOL_i=0;
    end
end
l_t(np,tm,num_seq)=length(t);%SOL(num_seq,:,1));
FC(np,tm,num_seq,1:2*N+np)=SOL(np,tm,num_seq,l_t(np,tm,num_seq),1:2*N+np);
Phivec=0;
% l(tm,num_seq)=length(SOL(tm,num_seq,1,:));
t_all(tm,num_seq,1:l_t(np,tm,num_seq))=t;
%% Figures
% Fig(tm,num_seq)=
% figure;
% hold on;
% axis([0 tmax 10^-7 max(max(SOL(np,tm,num_seq,:,:)))]);
% xlabel('Time t');ylabel('Population densities');
% pop=N-sum(2.^seq_prophage)+1;
% title(sprintf('Bacteriophage: %d, Lysogen insertion sequence [%s]',np,int2str(nonzeros(seq_prophage)')));box on;
% set(gca,'XScale','linear','XMinorTick','on','Xdir','normal','FontSize',14)
% set(gca,'YScale','log','YMinorTick','on','Ydir','normal','FontSize',14)
% plot(0,0,'k-','LineWidth',2);  plot(0,0,'k-.','LineWidth',2);   plot(0,0,'k--','LineWidth',2);
% legend('Lysogens','CRISPR Bacteria','Bacteriophage');
% for i=1:N+np
%     c=rand(1,3);
%     if i<=N
%         A=0;B=0;
%         A(1:l_t(np,tm,num_seq))=t_all(tm,num_seq,1:l_t(np,tm,num_seq));
%         B(1:l_t(np,tm,num_seq))=SOL(np,tm,num_seq,1:l_t(np,tm,num_seq),i);
%         plot(A,B,'-','Color',[c(1) c(2) c(3)],'LineWidth'...
%             ,2,'DisplayName',sprintf('Lysogens {%s}', d2b(i,:)));
%         A=0;B=0;
%         A(1:l_t(np,tm,num_seq))=t_all(tm,num_seq,1:l_t(np,tm,num_seq));
%         B(1:l_t(np,tm,num_seq))=SOL(np,tm,num_seq,1:l_t(np,tm,num_seq),N+i);
%         plot(A,B,'-.','Color',[c(1) c(2) c(3)],'LineWidth'...
%             ,2,'DisplayName',sprintf('CRISPR Bacteria {%s}', d2b(i,:)));
%     else
%         A(1:l_t(np,tm,num_seq))=t_all(tm,num_seq,1:l_t(np,tm,num_seq));
%         B(1:l_t(np,tm,num_seq))=SOL(np,tm,num_seq,1:l_t(np,tm,num_seq),N+i);
%         plot(A,B,'--','Color',[c(1) c(2) c(3)],'LineWidth'...
%                 ,2,'DisplayName',sprintf('Prophage %d',i-N));
%     end
% end
% SOL_all(tm,1:2^(np-1),1:m(tm,num_seq),1:2*N+np)=SOL(1:2^(np-1),1:m(tm,num_seq),1:2*N+np);
end
end
end
%% PLOTS ASSOCIATED WITH .m files
% for tm=1:length(arr_mx_t)
%     for num_seq=1:2^(np-1)
%         clear X T;
%         X(:,:)=SOL(3,tm,num_seq,:,:);T(:)=t_all(tm,num_seq,:);
%         Fig_popdyn(T,X(:,1:2*N,np));
%     end
% end
% 
% % BAR PLOT WITH ONLY END POINTS
% F(1:length(arr_mx_t),1:2^(np-1))=sum(FC(np,1:length(arr_mx_t),1:2^(np-1),N+1:2*N),4);
% Fig_endpoint([1 2 3], F)
G1(1:2*2^np+np)=FC(1,1,1,1:2*2^np+np);
G2(1:2*2^np+np)=FC(2,1,2,1:2*2^np+np);
G3(1:2*2^np+np)=FC(3,1,4,1:2*2^np+np);
G4(1:2*2^np+np)=FC(4,1,8,1:2*2^np+np);
G5(1:2*2^np+np)=FC(5,1,16,1:2*2^np+np);
G= [sum(G1(1:2)) 0 sum(G1(3:4)) 0 G1(5);...
sum(G2(1:4)) 0 sum(G2(5:8)) 0 sum(G2(9:10));...
sum(G3(1:8)) 0 sum(G3(9:16)) 0 sum(G3(17:19));...
sum(G4(1:16)) 0  sum(G4(17:32)) 0 sum(G4(33:36));...
sum(G5(1:32)) 0 sum(G5(33:64)) 0 sum(G5(65:69))];
figure
bar3(G);
toc