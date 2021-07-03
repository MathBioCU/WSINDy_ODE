function [t_dd,x_dd,rhs_dd] = generate_ddd(weights_dd,tags,tol_ode,x0,tspan,thresh)

if isempty(tol_ode)
    tol_ode = 10^-10;
end
if isempty(thresh)
    thresh = 5;
end
  
rhs_dd = build_vector_field(weights_dd, tags);
options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)),'Events',@(T,Y)myEvent(T,Y,thresh));
[t_dd,x_dd] =ode45(@(t,x) rhs_dd(x),tspan,x0,options);  
% 
% if and(norm(x_dd,inf)>thresh,thresh>0)
%     x_dd(abs(x_dd)>thresh) = 0;
% end

end

function [value, isterminal, direction] = myEvent(T, Y,thresh)
    value      = norm(Y) >= thresh;
    isterminal = 1;   % Stop the integration
    direction  = 0;
end