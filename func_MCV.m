function MCV=func_MCV(C)
global pF ps np N
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