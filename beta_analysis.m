clear all
tic
global npi np max_t beta
npi=2;
np=npi;
max_t=300; % Maximum time  TO BE GIVEN
beta_min=10^-10;
beta_max=10^-8;
freq=20;
beta_stepsize=(beta_max-beta_min)/(freq-1);
beta_vec=beta_min:beta_stepsize:beta_max;
k=1;
for beta=beta_vec
clearvars -except beta np npi max_t beta_vec k SOL_beta
Parameters_M4();
FC=func_MainODE();
NP=length(FC(:,1,1,1));TM=length(FC(1,:,1,1));NUM_SIM=length(FC(1,1,:,1));
TOT_POP=length(FC(1,1,1,:));
SOL_beta(k,1:NP,1:TM,1:NUM_SIM,1:TOT_POP)=FC;
X(k,1:TOT_POP)=FC(1,1,1,1:TOT_POP);
k=k+1;
end