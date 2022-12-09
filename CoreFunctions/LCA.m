function [choice, rt, x] = LCA(rho, k, beta, sgm, T0, thresh, dur, dt)
tau = .1; %s arbitrary defined. Once fixed, the time constant obtained from the free parameter k is relative to it. 
sizeVinput = size(rho);
if sizeVinput(1) > 1 error('Error: the size of input values has to be 1xN'); end
x = zeros(sizeVinput)';
rt = NaN;
choice = NaN;
%%
for ti = 1:(dur/dt)
    dx = (rho' - k*x(:,ti) - beta*x(:,ti))*dt/tau + randn(sizeVinput)'*sgm*sqrt(dt/tau);
    x(:,ti+1) = x(:,ti) + dx;
    x(x(:,ti+1) < 0,ti+1) = 0;
    % threshold detecting
    if max(x(:,ti+1)) >= thresh
        rt = ti*dt + T0;
        choice = find(x(:,ti+1) == max(x(:,ti+1)));
        break;
    end
end
x = x';
end
