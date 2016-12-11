% recursive function to calculate individual negative departure rate for a given list of units
function lmn=lambdan_v2(X,u,c,lambda) 
% X = CO, u: list of FORs for each unit, c: List of capacity values,
% lambda: failure rate in same units of time
if (length(c)==1) && (X == 0)
    lmn = lambda;
elseif   (X<0)
        lmn = 0;
elseif (length(c)==1) && (X==c(1))
        lmn = 0;

 elseif (X>sum(c))
         lmn = 0;
else
    
    if (prob_v2(X,u,c)~=0)
        
    up = u(2:length(u));
    cp = c(2:length(c));
    
    lmn= (prob_v2(X,up,cp)*(1-u(1))*(lambdan_v2(X,up,cp,lambda)+lambda) + prob_v2(X-c(1),up,cp)*u(1)*lambdan_v2(X-c(1),up,cp,lambda))/prob_v2(X,u,c);
    
    else 
        lmn=0;
    end

end


end