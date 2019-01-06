clear all
tic
V=[0 1 2];
x=perms(V);
[m,X]=nsumk(5,length(V));
i=1;
for j=1:m
if X(j,1)~=0%and(X(j,1)~=0,sum(X(j,:))==sum(V))
    Y(i,:)=X(j,:);
    i=i+1;
end
end
toc