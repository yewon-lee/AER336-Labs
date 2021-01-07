function plotfield2d(U)
% Plot 2d solution field associated with the wave equation problem.
% U is the vector of field values associated with equispaced interior grid 
% points of a unit square domain.
n1 = sqrt(length(U));
np = n1+2;
xp1 = linspace(0,1,np)';
x1 = xp1(2:end-1);

xp = [repmat(xp1,np,1),reshape(repmat(xp1',np,1),[],1)];
ibnd = find(any(abs(xp-0.5) > 0.5-1e-8,2));
ii = setdiff((1:np*np)',ibnd);

clf,
xmat = reshape(xp(:,1),[np,np]);
ymat = reshape(xp(:,2),[np,np]);
up = zeros(np*np,1);
up(ii) = U;
umat = reshape(up,[np,np]);
h = surf(xmat,ymat,umat);
set(h,'edgecolor','none');
colormap('jet');
shading interp;
zlim([-0.5,1]);
drawnow;
end