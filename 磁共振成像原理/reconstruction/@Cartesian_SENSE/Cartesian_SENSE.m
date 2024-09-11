function  res = Cartesian_SENSE(U,C)

res.adjoint = 0;
res.U = U;
res.C = C;
res.coil_rss = sqrt(sum(C.*conj(C),3));
res = class(res,'Cartesian_SENSE');

