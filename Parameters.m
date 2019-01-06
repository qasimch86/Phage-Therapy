function P=Parameters()
global K K2 pF
global b p r phi alpha beta gamma np ps pL eta pI ALPHA N
K=5*10^6; % Carrying Capacity - TABLE 17.1 STEPHEN D. ABEDON
K2=10^8;
p(1:np)=np:-1:1./np;%[1:np]./np;
r=1; % Growth rate constant for Bacteria with Biofilm
beta=6*10^-9; % Adsorption rate Biofilm Bacteria (levin and abedon)
b=150;       % Phage increase rate - Burst Size
alpha=10^-5;% Biofilm bacteria with prophage lyse (Moises Santillan)
ALPHA(1,1)=0;%alpha;
ALPHA(2:N,1)=alpha;% prophage induction rate
etaL_tild=5.56*10^-15*K2; % Formation of Biofilm from Planktonic; from HSP to HBP
phi=etaL_tild*K;
pF=10^-5; % Probability of CRISPR Failure
pI=1;% Probability of infection from previously existing phage to recently inserted lysogen
ps(1:np)=10^-5;% 10^-5; % Probability of gaining CRISPR Spacer
eta=0.01;%0.1;%0.05; % Formation of Planktonic from biofilm; from HB to HS
% ps(2)=ps(1);% 10^-5; % Probability of gaining CRISPR Spacer
% ps(3)=ps(2);% 10^-5; % Probability of gaining CRISPR Spacer
pL=0.01;
gamma=0.1; % Loss of biofilm Phage
end