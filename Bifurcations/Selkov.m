function dzdt = Selkov(t, z, parameters)
    x = real(z); y = imag(z);
    ayxxy = parameters.a*y + x*x*y;
    dzdt = -x + ayxxy + 1i*(parameters.b - ayxxy);
end

