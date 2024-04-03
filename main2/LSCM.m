function [ x, y ] = LSCM( a, x1, y1, n ) 
%two dimension Sine Logistic modulation map 
%a ,b are control parameters,a in [0,1] and b in [0,3] 
%n is the length of x or y 
%x1,y1 is the initial value,x and y in [0,1] 
 
x(1) = x1; 
y(1) = y1; 
 
for i=1:n-1 
x(i+1)=sin(pi*(4*a*x(i)*(1-x(i))+(1-a)*sin(pi*y(i)))); 
y(i+1)=sin(pi*(4*a*y(i)*(1-y(i))+(1-a)*sin(pi*x(i+1)))); 
end 
 
end 

