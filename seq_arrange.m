clear all;

np=3;
n=5;
c='';
for i=0:np
c=strcat(c,int2str(i));
end
A=(nchoosek(repmat(c, 1,n), n));
j=1;
for i=1:length(A(:,1))
    if str2num(A(i,1))>0
        B(j,:)=A(i,:);
        j=j+1;
    end
end
C=sort(unique(str2num(B)));
% for i=1:length(B(:,1))
%     if C(i)>=1000,C(i)<1