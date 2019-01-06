function MVL=func_MVL(V)
global pL np N pI MVL
MVL=zeros(N,N);
for i=1:np
    ai(i)=2^(i-1);
end
MVL(1,1)=sum(V);
for i = 2:N
    for j = 1:i
        if i==j
           a=dec2bin(i-1);la=length(a);
           b=zeros(np,1);
           for k=1:la
            b(np-k+1)=str2num(a(la-k+1));
           end
           for k=1:np
               if b(np-k+1)==0
                   if and(k<la,not(and(sum(b(1:np-k+1))>=1,sum(b(np-k+1:la))>=1)))
                        MVL(i,j)=MVL(i,j)+pI*V(k);
                   else
                        MVL(i,j)=MVL(i,j)+V(k);
                   end
               end
           end
        elseif i>j
          for k=1:np
              if i==j+ai(k)
                  for l=1:2^(k-1)
                      if or(rem(i,2^k)==0,rem(i,2^(k))==2^(k-1)+l)
                          MVL(i,j)=-pL*V(k);
                      end
                  end
              end
          end
        end
    end
end
% ml=Mvl