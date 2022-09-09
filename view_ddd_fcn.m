function [t_test,x_test,t_dd,x_dd,F_dd] = view_ddd_fcn(thresh,mult,w,ode_tol,x0_dd,tspan,tol_ode,ode_name,ode_params,tags,toggle_plot)

    n=length(x0_dd);
    
    x0_new=x0_dd+mult*(rand(size(x0_dd))-0.5);
    
    [~,t_test,x_test,~] = sim_ode(x0_new,tspan,tol_ode,ode_name,ode_params,'s');
    [t_dd,x_dd,F_dd] = generate_ddd(w,tags,ode_tol,x0_new,tspan,norm(x_test,thresh{1})*thresh{2});
    
    if toggle_plot
        if n>1
            for d=1:n
                legs = {};
                subplot(2*n,1,d) 
                plot(t_test,x_test(:,d),'k','linewidth',2);  legs{1} = 'true';
                hold on 
                plot(t_dd,x_dd(:,d),'b--','linewidth',2); legs{2} = 'learned';
                hold off
                legend(legs)
            end
    
            if n==2
                legs = {};
                subplot(2*n,1,n+1:2*n)
                plot(x_test(:,1),x_test(:,2),'k','linewidth',2);  legs{1} = 'true';
                hold on 
                plot(x_dd(:,1),x_dd(:,2),'b--','linewidth',2); legs{2} = 'learned';
                hold off
            elseif n==3
                legs = {};
                subplot(2*n,1,n+1:2*n)
                plot3(x_test(:,1),x_test(:,2),x_test(:,3),'k','linewidth',2);  legs{1} = 'true';
                hold on 
                plot3(x_dd(:,1),x_dd(:,2),x_dd(:,3),'g--','linewidth',2); legs{2} = 'learned';
                hold off
            end
        else
            plot(t_dd,x_dd,'g-'); legs{1} = 'DD';
            hold on 
            plot(t_test,x_test,'k');  legs{2} = 'true';
            hold off
            legend(legs)
        end
    end

end