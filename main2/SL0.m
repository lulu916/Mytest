function s=SL0(A, x, sigma_min ) 
 
    sigma_decrease_factor = 0.5; 
    A_pinv = pinv(A);%Î±Äæ¾ØÕó 
    mu_0 = 2; 
    L = 7; 
s = A_pinv*x; 
sigma = 2*max(max(abs(s))); 
 
% Main Loop 
while sigma>sigma_min 
    for i=1:L 
        delta = OurDelta(s,sigma); 
        s = s - mu_0*delta; 
        s = s - A_pinv*(A*s-x);   % Projection 
    end 
    sigma = sigma * sigma_decrease_factor; 
end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function delta=OurDelta(s,sigma) 
 
delta = s.*exp((-abs(s).^2)/sigma^2);
end 
