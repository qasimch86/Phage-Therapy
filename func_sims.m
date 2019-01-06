function Z=func_sims(np)
[y,x] = nsumk(np,np);
X=fliplr(x);
k=1;kx=1;
for i=1:length(X)
    % all possibilities without 0 as the initial value
        if not(X(i,1)==0)
            Xx(kx,:)=X(i,:);
            kx=kx+1;
        end
        % all possibilities without 0's in the injection
        for j=1:np
        if and(X(i,1)~=0,sum(X(i,:))==np)
            if j<np
            if and(X(i,j)==0,X(i,j+1)>0)
                k=k-1;
                break;
            end
            end
              Y(k,j)=X(i,j);
              flag=1;
        else
            flag=0;
            break;
        end
    end
    if flag==1
        k=k+1;
    end
    flag=1;
end
%         Z(1:2^(np-1),:)=Y(1:2^(np-1),:);
    k=1;
for j=1:np
    for i=1:2^(np-1)
        if Y(i,1)==j
            Z(k,1:np)=Y(i,1:np);
            k=k+1;
        end
    end
end
end