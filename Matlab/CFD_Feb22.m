 %------------------------
 pin='5a847a9dcffd0'; uid='vkumar'; fn=mfilename();
 websave('EA2.m','https://classes.mu2com.com/EA2.m'); EA2
 %------------------------
 
 
 P=struct('A',10e-3,'L',0.5, 'k',1000,'Ta',100,'Tb',500); 
 L=P.L; k=P.k; Area=P.A; Ta = P.Ta; Tb=P.Tb;  
  n=5;  dx=L/n; x1=dx/2; x2=L-dx/2; 
  x = linspace(x1,x2,n); 
  
   
  % aP TP = aW TW + aE PE + Su; 
  Su = 0; Sp=0; 
  for i=1:n
        aW = k*Area/dx; aE = k*Area/dx; Sp=0; Su=0; 
        if(i==1) aW = 0; Sp = -2*k*Area/dx; Su = 2*k*Area*Ta/dx; end
        if(i==n) aE = 0; Sp = -2*k*Area/dx; Su = 2*k*Area*Tb/dx; end
        aP = aW + aE - Sp;
        A(i,i)=aP; 
        if(i==1) A(i,i+1)=-aE; elseif(i==n) A(i,i-1)=-aW;else A(i,i-1)=-aW; A(i,i+1)=-aE; end
        b(i) = Su; 
  end

  T = A\b; 
  plot(x,T)
  return;
 Ta = 100; Tb = 500; 
 A = [
     300  -100  0    0    0 
     -100  200 -100  0    0 
      0   -100  200 -100  0 
      0     0  -100  200 -100 
      0     0   0   -100   300
     ];
 b = [200*Ta 0 0 0 200*Tb]'; 
 
 T = A\b
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 