function a = assortativity(k_in, k_out, k_accent_in, k_accent_out, N, k_mean, c)
    a = max(0, min(1, (k_accent_out.*k_in/(N*k_mean) + c*((k_accent_in - k_mean)./k_accent_out))*((k_out - k_mean)./k_in)));
end

