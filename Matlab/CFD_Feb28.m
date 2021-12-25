%------------------------
 pin='5a847a9dcffd0'; uid='vkumar'; fn=mfilename();
 websave('EA2.m','https://classes.mu2com.com/EA2.m'); EA2
 %------------------------
 clear; 
 %Input = jsondecode(fileread('InputCFD.json')); 
 %k = Input.k; Ar = Input.Area;  BCs = Input.BCs; Ta = BCs{1}.T; Tb = BCs{2}.T;
 k=1000; Ar=0.01; Ta=100; Tb=500; L=0.5; 

 k=0.5; Ar=1; Ta=100; Tb=200; L=2E-2; q=1000E3; 
 n = 5; dx = L/n; x = linspace(dx/2,L-dx/2,n);  c = k*Ar/dx; 
 for i = 1:n
     aW = c; aE = c; Su = q*Ar*dx; Sp = 0; 
     if(i==1) 
       A(i,i+1) = -aE; aW = 0; Sp=Sp-2*c; Su=Su + 2*c*Ta; 
     elseif(i==n)
       A(i,i-1) = -aW; aE = 0; Sp=Sp-2*c; Su=Su + 2*c*Tb;
     else
       A(i,i-1) = -aW; A(i,i+1) = -aE; 
     end
     aP = aW + aE - Sp;  A(i,i) = aP; b(i) = Su; 
 end
 T = A \ b'
 Texact = (Tb - Ta)/L * x + Ta; 
 plot(x,T, x, Texact, '*')
 
 
 
 
 
 
 
 
 
 
 
 
 
 