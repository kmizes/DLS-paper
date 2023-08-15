function [vel,speed] = getVelocityFromTrajectory(x,dt)

if nargin < 2
    dt = 1:size(x,1); % assume even spacing
    dt = (1/40):(1/40):(1/40)*size(x,1);
end

tol = .1/ (1/40); % smooth every .1 seconds, 1/40 is frame rate

vel = zeros(size(x));
for dim=1:size(x,2)
    v=diff(smooth(double(x(:,dim)),tol)) ./ diff(dt');
    vdt = mean([dt(1:end-1);dt(2:end)]);
    % smooth
    v = smooth(v,tol);
    % interp to original dimension
    vel(:,dim) = interp1(vdt, v, dt,'linear','extrap');
end

speed = sqrt(sum(vel.^2, 2));



end