%solve in the interior
function [u] = in_solver(src, N, tar, M, nu_y, mu)

iprec = 4;
nsource = N;
source = src';
ifcharge = 0;
charge = [];
ifdipole = 1;
dipstr = -1 / N * mu;
dipvec = nu_y';
ifpot = 0;
ifgrad = 0;
ifhess = 0;
ntarget = M;
targets = tar';
ifpottarg = 1;
ifgradtarg = 0;
ifhesstarg = 0;

[U] = rfmm2dpart(iprec,nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,...
        ntarget,targets,ifpottarg,ifgradtarg,ifhesstarg);
if U.ier ~= 0
    error('Error in FMM')
end

u = U.pottarg';

end