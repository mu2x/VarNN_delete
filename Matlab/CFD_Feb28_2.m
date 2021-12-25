clear all
close all
clc
commandwindow

input = struct('Area', 1, 'k', 0.5,'Ta', 100, 'Tb', 200, 'L', 2E-2, 'q', 1000E3);

n = 5; Ar = input.Area; L=input.L; k = input.k; Ta=input.Ta; Tb=input.Tb; q=input.q; 
dx = input.L/n;
x = dx/2:dx:input.L - dx/2;

%set sp = 0
%set su = q*A*L

%Su = (2*input.k*input.Area*input.Ta/dx) + (input.q*input.Area*dx)

for i = 1:n
    aW = input.k * input.Area / dx;
    aE = input.k * input.Area / dx;
    Su = input.q*input.Area*dx;
    Sp = 0;
    
    if i == 1
        A(i,i+1) = -aE;
        aW = 0;
        Sp = -2*input.k*input.Area/dx;
        %Su = 2*input.k*input.Area*input.Ta/dx;
        Su = (2*input.k*input.Area*input.Ta/dx) + (input.q*input.Area*dx)
    elseif i == n
        aE = 0; 
        A(i,i-1) = -aW;
        Sp = -2*input.k*input.Area/dx;
        %Su = 2*input.k*input.Area*input.Tb/dx;
        Su = (2*input.k*input.Area*input.Tb/dx) + (input.q*input.Area*dx)
    else
        A(i,i-1) = -aW;
        A(i,i+1) = -aE;
    end
    aP = aW + aE - Sp; 
    b(i) = Su;
    A(i,i) = aP;
    
end
B = b'; %transpose top get the column vector
A   %calculated A matric from for loop
Solution = A\B
%Texact = (input.Tb - input.Ta)/input.L*x+input.Ta;
c2=Ta; c1 = (Tb - c2)*k/L + q*L^2/(2*k) * k/L; 
Texact = -q*x.^2 / (2*k) + c1*x/(k) + c2; 


figure()
plot(x, Solution, x, Texact, 'o')
title('Temperature distribution')
axis on
grid on