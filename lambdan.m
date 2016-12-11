% recursive function to calculate failure rate for a number 
% of identical units
function lmn=lambdan(X,u,c,num,lambda) 
% X = CO, u:FOR, c: Capacity of unit, num: No. of units in system
% lambda: failure rate in same units of time
if (num == 1) && (X == 0)
    lmn = lambda;
elseif   (X<0)
        lmn = 0;
elseif (num == 1) && (X==c)
        lmn = 0;
 elseif (X>c*num)
         lmn = 0;
else
    lmn= (prob(X,u,c,num-1)*(1-u)*(lambdan(X,u,c,num-1,lambda)+lambda) + prob(X-c,u,c,num-1)*u*lambdan(X-c,u,c,num-1,lambda))/prob(X,u,c,num);
end


end