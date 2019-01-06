% user defined function for Standard Euler Method
function f=ODE_sol3(t,X0)
global b K r alpha ALPHA beta gamma np Ct T N eta Phivec
% global MVL MLI MCV MLV MVC VVL L C V VVC
T=T+1;
L(1:N,1)=X0(1:N);C(1:N,1)=X0(N+1:2*N);V(1:np,1)=X0(2*N+1:end);
G=1-(sum(L)+sum(C))/K;
    MVL=func_MVL(V);
    MLI=func_MLI(L);
    MCV=func_MCV(C);
    MLV=func_MLV(L);
    MVC=func_MVC(V);
% for i=1:N;  VVL(i,1)=(MVL(i,:)*L); end;
% for i=1:np; VLI(i,1)=(MLI(i,:)*L); end;
% for i=1:N;  VVC(i,1)=(MVC(i,:)*C); end;
% for i=1:np; VCV(i,1)=(MCV(i,:)*V); end;
% for i=1:np; VLV(i,1)=(MLV(i,:)*V); end;
VVL(1:N,1) =MVL(1:N,:)*L;
VLI(1:np,1)=MLI(1:np,:)*L;
VVC(1:N,1) =MVC(1:N,:)*C;
VCV(1:np,1)=MCV(1:np,:)*V;
VLV(1:np,1)=MLV(1:np,:)*V;
L1=G.*(r.*L+Phivec)-ALPHA.*L-beta.*VVL-eta.*L;
C1=G.*r.*C-beta.*VVC-eta.*C;
V1=b.*alpha.*VLI+b*beta*VCV-gamma*V+b*beta*VLV;
%   if floor(t)>floor(Ct)
%       	fprintf('Current time:%d. \r',floor(t))
%       	Ct=floor(t);
% %         MLV
% %         MCV
% %         pause;
%   end
f=[L1;C1;V1];