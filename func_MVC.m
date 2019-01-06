function MVC=func_MVC(V)
global pF ps np
N=2^np;
for i=1:np
    ai(i)=2^(i-1);
end
MVC=zeros(N);
MVC(1,1)=sum(V);
for i = 2:N
    for j = 1:i
        if i==j
           a=dec2bin(i-1);la=length(a);
           b=zeros(np,1);
           for k=1:la
            b(np-k+1)=str2num(a(la-k+1));
           end
           lb=length(b);
           for k=1:lb
               if b(lb-k+1)==0
                    MVC(i,j)=MVC(i,j)+V(k);
               else
                    MVC(i,j)=MVC(i,j)+pF*V(k);
               end
           end
        elseif i>j
          for k=1:np
              if i==j+ai(k)
                  for l=1:2^(k-1)
                      if or(rem(i,2^k)==0,rem(i,2^(k))==2^(k-1)+l)
                          MVC(i,j)=-ps(k)*V(k);
                      end
                  end
              end
          end
        end
    end
end