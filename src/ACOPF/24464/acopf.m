function [out,co,co2,obj,status,times] = acopf(reg)
% ACOPF local solve (IPOPT), with generalized ADMM consistency terms.
% Only the consistency residuals + their derivatives are modified for mode='primal'.
% Everything else is unchanged.

filenamesp = ['mpc', reg];
mpc_struct = evalin('base', filenamesp);
mpc = preprocess(mpc_struct);

t0 = cputime; %#ok<NASGU>
co2 = 1;
ngens = size(mpc.gen,1);
options.auxdata = {mpc};
if ngens~=0
    co = mpc.gencost(:, [5 6 7]);
else
    co = [];
    co2 = 0;
end

options.ipopt.print_level              = 0;
options.ipopt.tol                      = 1e-6;
options.ipopt.max_iter                 = 500;
options.ipopt.dual_inf_tol             = 1e-1;
options.ipopt.compl_inf_tol            = 1e-5;
options.ipopt.acceptable_tol           = 1e-8;
options.ipopt.acceptable_compl_inf_tol = 1e-3;
options.ipopt.mu_strategy              = 'adaptive';
options.ipopt.nlp_scaling_method       = 'none';

x1 = initialx0(options.auxdata);
[options.lb, options.ub] = bounds(options.auxdata);
[options.cl, options.cu] = constraintbounds(options.auxdata);

funcs.objective            = @objective;
funcs.constraints          = @constraints;
funcs.gradient             = @gradient;
funcs.jacobian             = @jacobian;
funcs.jacobianstructure    = @jacobianstructure;
funcs.hessian              = @hessian;
funcs.hessianstructure     = @hessianstructure;

[out, info_opf] = ipopt_auxdata(x1, funcs, options);
obj    = objective(out, options.auxdata);
status = info_opf.status;

times= info_opf.cpu;
