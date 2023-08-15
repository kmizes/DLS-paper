function [xinterp] = interpTrajectories(x, landmarks, targets)

assert(isequal(size(x), size(y)));
assert(numel(landmarks)==numel(targets),'Interp points must be same size');
assert(length(targets)>1, 'Need multiple points to warp to');

xinterp = x(1:(landmarks(1)-1));

for i = 1:(length(landmarks)-1)

xwarp = linspace(landmarks(i), landmarks(i+2), targets(i+2)-targets(i));
xinterp = [xinterp , interp1(xwarp, x(landmarks(i):landmarks(i+1)), targets(i):targets(i+1))];

end


end