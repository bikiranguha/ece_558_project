% recursive function to calculate individual positive departure rate for a given list of units
function lmp = lambdap_v2(X,u,c,mu) 
% X = CO, u: list of FORs for each unit, c: List of capacity values,
% mu: repair rate in same units of time
if (length(c)==1) && (X == 0)
    lmp= 0;
elseif   (X<0)
        lmp = 0;
elseif (length(c)==1) && (X==c(1))
        lmp = mu;
 elseif (X>sum(c))
         lmp = 0;
else
    
    if (prob_v2(X,u,c)~=0)
        
    up = u(2:length(u));
    cp = c(2:length(c));
    
    lmp= (prob_v2(X,up,cp)*(1-u(1))*lambdap_v2(X,up,cp,mu) + prob_v2(X-c(1),up,cp)*u(1)*(lambdap_v2(X-c(1),up,cp,mu)+mu))/(prob_v2(X,u,c));
    
    else
        lmp = 0;
    end
    
end 

end