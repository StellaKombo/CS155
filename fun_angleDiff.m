function angleDiff_out = fun_angleDiff(a,b)
    dif = mod(b-a+pi,2*pi);
    if(dif < 0)
        dif = dif + 2*pi;
    end
    angleDiff_out = dif - pi;
end
