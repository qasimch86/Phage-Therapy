function Mlv=func_MLV(L)
global pL np N pI
Mlv=zeros(np,np);
for i=1:np
    pl(i)=pL;
end
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
                Mlv(i,i)=Mlv(i,i)+(1-pl(i))*pI*L(j);
            else
                Mlv(i,i)=Mlv(i,i)+(1-pl(i))*L(j);
            end
        end
    end
end