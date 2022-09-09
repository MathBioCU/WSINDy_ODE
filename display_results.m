function [err_wsindy,err_sindy,tp_w,tp_s,t_dd,x_dd,F_dd,lambda_hat_w,lambda_hat_s,resid,resid_true] = ...
    display_results(...
    w_sparse,true_nz_weights,w_sparse_sindy,loss_wsindy,...
    loss_sindy,lambda,noise_ratio,noise_ratio_obs,sigma,ET_wsindy,ET_sindy,xobs,x,tobs,t,grids,Gs,RTs,bs,Theta_0,...
    mts,pts,toggle_print_w,toggle_plot,toggle_plot_ddd,thresh,mult,toggle_plot_resid,...
    toggle_plot_loss,toggle_plot_derivs,toggle_plot_approx_sys,toggle_plot_fft,bweaks,useFD,dxobs_0,...
    tol_ode,x0,ode_name,ode_params,tags,vs)

    n = length(x0);
    
    err_wsindy = [norm(w_sparse(:)-true_nz_weights(:));norm(w_sparse(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights(:));
    err_sindy = [norm(w_sparse_sindy(:)-true_nz_weights(:));norm(w_sparse_sindy(true_nz_weights~=0)-true_nz_weights(true_nz_weights~=0))]/norm(true_nz_weights(:));        
    
    lambda_loss_wsindy = []; lambda_loss_sindy = [];
    for nn=1:n
        if length(lambda)>1
            indtemp = find(loss_wsindy{nn}(4,:)==0,1)-1; if isempty(indtemp);indtemp=length(lambda);end
            lambda_loss_wsindy(nn,:) = [indtemp lambda(indtemp)];
            indtemp = find(loss_sindy{nn}(4,:)==0,1)-1; if isempty(indtemp);indtemp=length(lambda);end
            lambda_loss_sindy(nn,:) = [indtemp lambda(indtemp)];
        else
            lambda_loss_wsindy(nn,:) = [1 lambda];
            lambda_loss_sindy(nn,:) = [1 lambda];
        end
    end
    
    disp(['log10 2norm err (all weights) (WSINDy)=',num2str(log10(err_wsindy(1)))])
    disp(['log10 2norm err (all weights) (SINDy)=',num2str(log10(err_sindy(1)))])
    disp(['log10 2norm err (true nz weights) (WSINDy)=',num2str(log10(err_wsindy(2)))])
    disp(['log10 2norm err (true nz weights) (SINDy)=',num2str(log10(err_sindy(2)))])
    tp_w = tpscore(w_sparse,true_nz_weights);
    tp_s = tpscore(w_sparse_sindy,true_nz_weights);
    disp(['TPR (WSINDy, SINDy)=',num2str([tp_w tp_s])])
    disp(' ')
    disp(['Noise_ratio, sigma =',num2str([noise_ratio_obs sigma])])
    disp(['Run time (WSINDy) =',num2str(ET_wsindy)])
    disp(['Run time (SINDy)  =',num2str(ET_sindy)])
    lambda_hat_w  = lambda_loss_wsindy(:,2)';
    lambda_hat_s  = lambda_loss_sindy(:,2)';
    disp(['lambda_hat_w  =', num2str(lambda_hat_w)])
    disp(['lambda_hat_s  =', num2str(lambda_hat_s)])
    disp(' ')
    disp(['Num timepoints =',num2str(size(xobs,1))])
    disp(['Num Trial Fcns =',num2str(size(w_sparse,1))])
    disp(['Num Basis Fcns =',num2str(cellfun(@(x) length(x),grids)')])
    disp(['Basis degrees  =',num2str(pts')])
    disp(['Basis supports =',num2str(2*mts'+1)])
    disp(['log10(cond(G)) =',num2str(cellfun(@(x) log10(cond(x)),[Gs{:} {Theta_0}]))])
    disp(' ')
    
    if toggle_print_w == 1
        disp(['true, wsindy, sindy'])
        disp([true_nz_weights w_sparse w_sparse_sindy])
    elseif toggle_print_w == 2
        disp(['true, wsindy, sindy, lib term'])
        nzs = find(true_nz_weights);
        nzsw = find(w_sparse);
        nzss = find(w_sparse_sindy);
        nzs=unique([nzs;nzsw;nzss]);
        rems = mod(nzs,size(w_sparse,1));
        rems(rems==0) = size(w_sparse,1);
        wws = [true_nz_weights(nzs) w_sparse(nzs) w_sparse_sindy(nzs) tags(rems,:)];
        if n==1
            disp(wws);
        elseif n==2
            cut = find(nzs>size(w_sparse,1),1)-1;
            disp(wws(1:cut,:))
            disp('------------')
            disp(wws(cut+1:end,:))
        elseif n==3
            cuts = [find(nzs>size(w_sparse,1),1)-1;find(nzs>2*size(w_sparse,1),1)-1];
            disp(wws(1:cuts(1),:))
            disp('------------')
            disp(wws(cuts(1)+1:cuts(2),:))
            disp('------------')
            disp(wws(cuts(2)+1:end,:))
        end     
    elseif toggle_print_w == 3
        nzs = [];
        for k=1:n
            nzs = [nzs;find(true_nz_weights(:,k))];
            nzs = [nzs;find(w_sparse(:,k))];
            nzs = [nzs;find(w_sparse_sindy(:,k))];
        end
        nzs = unique(nzs);
        disp(['true, wsindy, sindy'])
        disp([true_nz_weights(nzs,:) w_sparse(nzs,:) w_sparse_sindy(nzs,:)])    
    end
    
    if toggle_plot>0
        figure(1); clf
%         set(gcf, 'units','normalized','outerposition',[0 0.5 0.5 0.5])
        for nn=1:n
        subplot(2,n,nn)
        plot(tobs,xobs(:,nn),'r-',tobs(grids{nn}+mts(nn)),mean(xobs(:,nn))*ones(length(grids{nn}),1),'.k')
        subplot(2,n,n+nn)
        spy(RTs{nn})
        end
    end
    
    if or(toggle_plot_ddd==1,toggle_plot_ddd==2)
        figure(2); clf
%         set(gcf, 'units','normalized','outerposition',[0.5 0.5 0.5 0.5])
        if toggle_plot_ddd == 1
            w = w_sparse;
        elseif toggle_plot_ddd == 2
            w = w_sparse_sindy;
        end
        [~,~,t_dd,x_dd,F_dd] = view_ddd_fcn(thresh,mult,w,tol_ode,x0,tobs(1:5:end),tol_ode,ode_name,ode_params,tags,toggle_plot_ddd);
        title(num2str(tobs(end)))
    elseif toggle_plot_ddd>2
        figure(2); clf
%         set(gcf, 'units','normalized','outerposition',[0.5 0.5 0.5 0.5])
        np = toggle_plot_ddd;
        F_dd = build_vector_field(w_sparse,tags);
        rhs = build_vector_field(true_nz_weights,tags);
        xmin = min(xobs);
        xmax = max(xobs);
        vec_grids = {};
        for nn=1:n
            vec_grids{nn} = linspace(xmin(nn),xmax(nn),np);
        end
        if n==1
            errs = abs(F_dd(vec_grids{1})-rhs(vec_grids{1}))/rms(rhs(vec_grids{1}));
            plot(vec_grids{1},errs)
        elseif n==2
            [xx,yy] = meshgrid(vec_grids{:});        
            vec_points = [xx(:)';yy(:)'];
            errs = vecnorm(F_dd(vec_points)-rhs(vec_points),2,1)/rms(vecnorm(rhs(vec_points),2,1));
            errs = reshape(errs,np,np);
            imagesc(vec_grids{:},errs)
            hold on
            plot(xobs(:,1),xobs(:,2),'r.')
            colorbar
        elseif n==3
             errs = vecnorm(F_dd(x')-rhs(x'),2,1)/rms(rhs(x'));
             scatter3(x(:,1),x(:,2),x(:,3),1.5,errs)
             colorbar
        end
        t_dd = [];
        x_dd = [];
    else
        t_dd = [];
        x_dd = [];
        F_dd = [];
    end
    
    if toggle_plot_resid
        figure(3); clf
%         set(gcf, 'units','normalized','outerposition',[0.5 0 0.5 0.5])
        grids_r = grids;
        for nn=1:n
            grids_r{nn} = grids{nn}+mts(nn);
        end
        for nn=1:n
            subplot(2,n,nn)
            resid{nn} = (bs{nn}-Gs{nn}*w_sparse(:,nn))/norm(bs{nn});
            plot(t,x(:,nn)*(noise_ratio/max(x(:,nn))),'r--',t(grids_r{nn}),resid{nn},'k.')
            ylim([-noise_ratio noise_ratio])
            legend({'true sol.','WSINDy res'})
            subplot(2,n,n+nn)
            resid_true{nn} = (bs{nn}-Gs{nn}*true_nz_weights(:,nn))/norm(bs{nn});
            plot(t,x(:,nn)*(noise_ratio/max(x(:,nn))),'r--',t(grids_r{nn}),resid_true{nn},'k.')
            ylim([-noise_ratio noise_ratio])
            legend({'true sol.','True res'})
        end
    else
        resid = [];
        resid_true = [];
    end
    
    if and(length(lambda)>1,toggle_plot_loss)
        figure(4); clf
%         set(gcf, 'units','normalized','outerposition',[0 0 0.5 0.5])
        for nn=1:n
            subplot(n,1,nn)
            semilogx(lambda,loss_wsindy{nn}(1,:),'ob-',lambda_loss_wsindy(nn,2),loss_wsindy{nn}(3,lambda_loss_wsindy(nn,1)),'rx'); hold on;
            if ~isempty(dxobs_0)
                semilogx(lambda,loss_sindy{nn}(1,:),'ok-',lambda_loss_sindy(nn,2),loss_sindy{nn}(3,lambda_loss_sindy(nn,1)),'rx');hold off;
                legend({'wsindy','','sindy',''}, 'location', 'best')
            else
                legend({'wsindy',''}, 'location', 'best')
            end
        end
    end
    
    if and(toggle_plot_derivs>0,~isempty(bweaks{1}))
        dt = mean(diff(t));
        figure(5)
        for nn=1:n
            subplot(n,toggle_plot_derivs,toggle_plot_derivs*nn-(toggle_plot_derivs-1))
            plot(t,bweaks{nn},t(1:end-1),diff(x(:,nn))/dt)
            legend({'WSINDy',''})
            if and(~isempty(dxobs_0),toggle_plot_derivs==2)
                subplot(n,2,2*nn)
                plot(t(max(useFD,1)+1:end-max(useFD,1)),dxobs_0(:,nn),t(1:end-1),diff(x(:,nn))/dt)
                legend({'SINDy',''})
            end
        end
    end
    
    if toggle_plot_approx_sys
        figure(6); clf;
        x0obs = mean(xobs-cumtrapz(tobs,Theta_0*w_sparse));
        disp(['Approx Error =',num2str(mean(cumtrapz(tobs,Theta_0*w_sparse)+x0obs-x))])
        x_approx = cumtrapz(tobs,Theta_0*w_sparse)+x0obs;
        legs = {};
        hold on 
        for nn=1:n
            plot(t,x(:,nn),'b-',tobs,xobs(:,nn),'r.',tobs,x_approx(:,nn), 'g--')
            legs = {legs{:},'x','x_{obs}','x_{approx}'};
        end
        legend(legs,'location','best')
    end

    if toggle_plot_fft>0
        figure(7);clf
        for j=1:n
            m = (length(vs{j,1})-1)/2;
            Ufft = abs(fft(xobs(:,j)));
            Ufft = Ufft(floor(end/2):end);
            L = length(Ufft)-1;
            ks = -L:L;
            Ufft = [Ufft; flipud(Ufft(1:end-1))]/max(Ufft);
            subplot(n,1,j)
                semilogy(ks,Ufft)
                hold on
                Cfs_ffts = fft([zeros(1,length(tobs)-2*m-1) vs{j,1}]);
                Cfs_ffts = abs(Cfs_ffts(floor(end/2):end));
                Cfs_ffts = [Cfs_ffts fliplr(Cfs_ffts(1:end-1))];
                Cfs_ffts = Cfs_ffts/max(Cfs_ffts);
                semilogy(ks,Cfs_ffts)
                [corner,~] = findcornerpts(xobs(:,j),tobs);
                k = corner(2);
                semilogy([-k k],Ufft(max(L+1-k,1))*[1 1],'o','markersize',12)
                hold off     
                ylim([min(Ufft)*0.1 max(Ufft)])
                xlim([min(ks) max(ks)])
                legend({'$\mathcal{F}(y)$','$\mathcal{F}(\phi)$','$k^*$'},'interpreter','latex','fontsize',14)
                title(['coord ',num2str(j)])
        end
        xlabel('k')
    end



end

