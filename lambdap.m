% recursive function to calculate repair rate for a number 
% of identical units
function lmp=lambdap(X,u,c,num,mu) 
% X = CO, u:FOR, c: Capacity of unit, num: No. of units
% mu: repair rate in same units of time
if (num == 1) && (X == 0)
    lmp = 0;
elseif   (X<0)
        lmp = 0;
elseif (num == 1) && (X==c)
        lmp = mu;
 elseif (X>c*num)
         lmp = 0;
else
    lmp= (prob(X,u,c,num-1)*(1-u)*lambdap(X,u,c,num-1,mu) + prob(X-c,u,c,num-1)*u*(lambdap(X-c,u,c,num-1,mu)+mu))/prob(X,u,c,num);
end


end