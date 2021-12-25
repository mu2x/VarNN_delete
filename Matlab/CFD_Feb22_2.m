%------------------------
 pin='5a847a9dcffd0'; uid='vkumar'; fn=mfilename();
 websave('EA2.m','https://classes.mu2com.com/EA2.m'); EA2
 %------------------------
 clear; 
 Input=struct('Area',10E-3, 'k',1000,'Ta',100,'Tb',500,'L',0.5);
 n = 5; dx = Input.L / n; x = dx/2:dx:Input.L-dx/2; 
 k = Input.k; Ar = Input.Area; 
 for i = 1:n
     aW = k * Ar / dx; aE = k * Ar / dx; Su = 0; Sp = 0; 
     if(i==1) 
       A(i,i+1) = -aE; aW = 0; Sp = -2*k*Ar/dx; Su = 2*k*Ar* Input.Ta /dx; 
     elseif(i==n)
       A(i,i-1) = -aW; aE = 0; Sp = -2*k*Ar/dx; Su = 2*k*Ar* Input.Tb /dx;
     else
       A(i,i-1) = -aW; A(i,i+1) = -aE; 
     end
     aP = aW + aE - Sp;  A(i,i) = aP; b(i) = Su; 
 end
 T = A \ b'; Texact = (Input.Tb - Input.Ta)/Input.L * x + Input.Ta; 
 plot(x,T, x, Texact, '*')
 
 
 
 
 
 
 
 
 
 
 
 
 
 