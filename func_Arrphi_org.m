clear all
tic
np=5;
N=2^np;
C=zeros(N,1);
phi=2.78;
seq_prophage=[1 4];
freq=length(seq_prophage);
d2b=(dec2bin(0:N-1));% decimal to binary conversion of each possible permutation of number of prophage
Arrphi=zeros(N,freq);
for k=1:freq
if k>1;   add_xprophage=sum(seq_prophage(1:k-1)); else add_xprophage=0; end
    a=np-add_xprophage-seq_prophage(k)+1;
    b=np-add_xprophage;
    for i=1:N
    d2bsum=0;
    for j=a:b
        d2bsum=d2bsum+str2num(d2b(i,j));
    end
    if d2bsum==seq_prophage(k)
        new_prophage_val(k)=i;
        break;
    end
    end
    Arrphi(new_prophage_val(k),k)=phi;
end
toc