function handle = make_GPUhandle()
if gpuDeviceCount > 0
    handle = @gpuArray;
else 
    handle = @(x) x;
end
end

