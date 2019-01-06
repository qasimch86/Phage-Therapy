clear all
tic
np=3;
pF=0.5;
N=2^np;
C=zeros(N,1);
for i=1:np
    ps(i)=10^-5;
end
for i=1:N
    C(i,1)=i;
end
MCV=zeros(np,np);
for i = 1:np
    b=0;
    for j = 1:N
        a=dec2bin(j-1);  la=length(a);
        for k=1:la
            b(np-k+1)=str2num(a(la-k+1));
        end
        lb=length(b);
        if b(lb-i+1)==0
            MCV(i,i)=MCV(i,i)+(1-ps(i))*C(j);
        else
            MCV(i,i)=MCV(i,i)+pF*C(j);
        end
    end
end
toc