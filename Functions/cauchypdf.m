function pdf = cauchypdf(x, mu, gamma)
    pdf = 1 ./ (pi*gamma*(1 + power((x - mu)/gamma, 2)));
end

