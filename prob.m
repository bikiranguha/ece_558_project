% recursive function to calculate individual probability for a number 
% of identical units
function pr = prob(X,u,c,num) 
% X = CO, u:FOR, c: Capacity of unit, num: No. of units in system after
% addition of unit
if (num == 1) && (X == 0)
    pr= 1-u;
elseif   (X<0)
        pr = 0;
elseif ( num == 1) && (X==c)
        pr = u;
 elseif (X>c*num)
         pr = 0;
else
    pr= prob(X,u,c,num-1)*(1-u) + prob(X-c,u,c,num-1)*u;
end


end
