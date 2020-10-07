function handle = makeGPUinitHandle()
if gpuDeviceCount > 0
    handle = @gpuArray;
else 
    handle = @(x) x;
end
end

