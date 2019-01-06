clear all
tic
np=3;
N=2^np;
L=zeros(N,1);
pI=0.2;
pL=0.1;
for i=1:N
    L(i,1)=i;
end
L=[4;2;5;3;7;4;2;4];%4.5*10^6;30.1482];
Mlv=zeros(np,np);
for i = 1:np
    b=0;
    for j = 1:N
        a=dec2bin(j-1);  la=length(a);
        for k=1:la
            b(np-k+1)=str2num(a(la-k+1));
        end
%         lb=length(b);
        if b(np-i+1)==0
            if i<la
                Mlv(i,i)=Mlv(i,i)+(1-pL)*pI*L(j);
            else
                Mlv(i,i)=Mlv(i,i)+(1-pL)*L(j);
            end
        end
    end
end
toc