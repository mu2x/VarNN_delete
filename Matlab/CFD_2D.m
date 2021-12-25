%------------------------
 pin='5a847a9dcffd0'; uid='vkumar'; fn=mfilename();
 websave('EA2.m','https://classes.mu2com.com/EA2.m'); EA2
 %------------------------
 clear; 
 %Input = jsondecode(fileread('InputCFD.json')); 
 %k = Input.k; Ar = Input.Area;  BCs = Input.BCs; Ta = BCs{1}.T; Tb = BCs{2}.T;

 Lx=0.5; Ly=0.5; n = 20; m = 5; nn=n*m; % Geom descretization/mesh
 Ta=100; Tb=200; Tc=200; Td=400; % BCs
 k=0.5;  % Properties
 
 dx = Lx/n;  dy = Ly/m;   
 x = linspace(dx/2,Lx-dx/2,n); y = linspace(dy/2,Ly-dy/2,m); [X Y]=meshgrid(x,y); 
 A=zeros(nn,nn); b=zeros(1,nn); 
 
 dz=1;  G=k; q=0; 
 for j=1:m
     for i=1:n
         P = (j-1)*n+i; W = P-1; E = P+1; N=P+n; S=P-n; 
         Aw = dy*dz; Ae=Aw; An=dx*dz; As = An; Su=0; Sp=0; 
         aW = G*Aw/dx; aE = G*Ae/dx; aN = G*An/dy; aS = G*As/dy;
         if(i>1) A(P,W)=-aW; else 
             %aW=0; Sp=Sp-2*G*Aw/dx; Su=Su + 2*G*Aw/dx*Ta;
             aW=0; Sp=Sp; Su=Su + q*Aw*dx;
         end
         if(i<n) A(P,E)=-aE; else aE=0; Sp=Sp-2*G*Ae/dx; Su=Su + 2*G*Ae/dx*Tb;end
         if(j>1) A(P,S)=-aS; else aS=0; Sp=Sp-2*G*As/dy; Su=Su + 2*G*As/dy*Tc;end
         if(j<m) A(P,N)=-aN; else aN=0; Sp=Sp-2*G*An/dy; Su=Su + 2*G*An/dy*Td;end
         aP = aW + aE + aS + aN - Sp; 
         A(P,P) = aP; b(P)=Su; 
     end
 end
 
 d=A\b';  
 
 for j=1:m; for i=1:n; in=(j-1)*n+i; d2(j,i)=d(in); end; end
 contourf(X,Y,d2); grid on
 
 