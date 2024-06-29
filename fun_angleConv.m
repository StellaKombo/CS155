function angleConv_out = fun_angleConv(angle)
    angleConv_out = mod(fun_ConstrainAngle(angle), 2*pi);
end
