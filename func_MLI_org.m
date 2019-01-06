clear all
tic
np=3;
N=2^np;
p=[1 0.5 0.33];
Mli=zeros(np,N);
Lv=zeros(np,1);
for i = 1:np
    flag=0;
    for j = 2^(i-1)+1:N
        if flag~=0
            flag=flag-1;
            continue;
        else
            Mli(i,j:j+2^(i-1)-1)=1;
            flag=2^i-1;
        end
    end
end
for i=2:N
    Mi(i,i)=p(sum(Mli(:,i)));
end    
for i=1:np
    for k=1:N
            LI(i,k)=sum(Mli(i,:)'.*Mi(:,k));
    end
end
toc