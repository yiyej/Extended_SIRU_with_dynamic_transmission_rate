function c = bisectionMethod(f,a,b,err)
%f is a vector which contains the coefficients of polynomial
c=(a+b)/2;
while abs(ppval(f, c))>err
    if (ppval(f, c)*ppval(f, a))<0
        b=c;
    elseif (ppval(f, c)*ppval(f, b))<0
        a=c;
    else
        disp('No root in the given interval or multiple roots.')
        c = NaN;
        break
    end
    c=(a+b)/2;
end
end 
