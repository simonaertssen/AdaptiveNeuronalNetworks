function [tout,xout] = DOPRI45(originalfunc,ta,tb,x0,h,p)
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    dim = size(x0);
    xout = zeros(dim(1),npts); xout(:,1) = x0;    

    tout = linspace(ta,tb,npts);
    % Make new function handle to improve speed!
    func = @(t, x) originalfunc(t, x, p.e, p.K/p.N, p.a_n);
    
    for i = 1:(npts-1)
        xout(:,i+1) = DOPRIstep(func,tout(i),xout(:,i),h);
    end
end