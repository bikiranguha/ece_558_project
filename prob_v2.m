% recursive function to calculate individual probability for a given list of units
function pr = prob_v2(X,u,c) 
% X = CO, u: list of FORs for each unit, c: List of capacity values

if (length(c)==1) && (X == 0)
    pr= 1-u(1);
elseif   (X<0)
        pr = 0;
elseif (length(c)==1) && (X==c(1))
        pr = u(1);

 elseif (X>sum(c))
         pr = 0;
else
    
   up = u(2:length(c));
   cp = c(2:length(c));
    
    pr= prob_v2(X,up,cp)*(1-u(1)) + prob_v2(X-c(1),up,cp)*u(1);
end


end
