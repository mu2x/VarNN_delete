%------------------------
 pin='5ba542d0e4d33'; uid='vkotteda'; fn=mfilename();
 websave('EA2.m','https://classes.mu2com.com/EA2.m'); EA2
 %------------------------
  clear; 
 %Input = jsondecode(fileread('InputCFD.json')); 

 Phi1=1; Phi2=0; u=0.1; G=0.1; rho=1;  Ar=1; L=1; n=5; 
 u=3.5; n=20;
 dx = L/n; x = linspace(dx/2,L-dx/2,n);  A=zeros(n,n); b=zeros(n,1);
 for i = 1:n
     Dw = G/dx; De = G/dx; Fw = rho*u; Fe=rho*u; 
     %aW = Dw + Fw/2; aE = De - Fe/2; Sp=0; Su = 0; 
     aW = Dw + max(Fw,0); aE = De + max(0,-Fe); Sp=0; Su = 0; 
     if(i==1) 
       A(i,i+1) = -aE; aW = 0; Su = Su + (2*Dw+Fw)*Phi1; Sp = Sp -(2*Dw+Fw); 
     elseif(i==n)
       A(i,i-1) = -aW; aE = 0; Su = Su + (2*De-Fe)*Phi2; Sp = Sp -(2*De-Fe);
     else
       A(i,i-1) = -aW; A(i,i+1) = -aE; 
     end
     aP = aW + aE + (Fe-Fw) - Sp;  A(i,i) = aP; b(i) = Su; 
     %[aW aE Su Sp aP]
 end
 Phi = A \ b;
 Phiexact = Phi1 + (Phi2-Phi1)*(exp(u*rho*x/G)-1)/(exp(u*rho*L/G)-1); 
 plot(x,Phi, x, Phiexact)
 
 