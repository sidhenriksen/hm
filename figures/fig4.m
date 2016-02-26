function fig4()
    run_parallel = 1;
    bootstrap_mode = 1;

    silent = 1;
    bem = BEMunit('x0',0,'y0',0,'silent',silent);
    bem.deg_per_pixel=0.025;
    bem.outputNL = @(x)(x.^2);
    
    % Correlated random dot stereogram
    crds = RDS();
    crds.Nx = bem.Nx;
    crds.Ny = bem.Ny;    
    crds.dotsize=round(0.1./bem.deg_per_pixel);%0.1 deg dot radius
    
    % Half-matched random dot stereogram
    hmrds = crds;
    hmrds.toggle_match = 1;
    hmrds.correlation = 0;
    
    % Uncorrelated random dot stereogram
    urds = crds;
    urds.correlation = 0;
    
    densities = [0.05,0.24];
    
    sxs = linspace(0.01,0.3,21);
    
    N = 0;
    
    
    
    model_ratio = zeros(length(sxs),length(densities));
    
    for k = 1:length(densities);
        crds.density = densities(k);
        hmrds.density = densities(k);
        urds.density = densities(k);
        
        for j = 1:length(sxs);
            current_sx = bem.subunits(1).rf_params.left.sx;
            bem = bem.rescale(sxs(j)/current_sx);             
            
            if bootstrap_mode && ~N
                bem = bem.load_bootstrap(crds);
                C = mean(bem.simulate_spatial(crds,N,bootstrap_mode,run_parallel));
                bem = bem.load_bootstrap(hmrds);
                HM = mean(bem.simulate_spatial(hmrds,N,bootstrap_mode,run_parallel));
                bem = bem.load_bootstrap(urds);
                U = mean(bem.simulate_spatial(urds,N,bootstrap_mode,run_parallel));
            else
                
            C = mean(bem.simulate_spatial(crds,N,bootstrap_mode,run_parallel));
            HM = mean(bem.simulate_spatial(hmrds,N,bootstrap_mode,run_parallel));
            U = mean(bem.simulate_spatial(urds,N,bootstrap_mode,run_parallel));
            end

            model_ratio(j,k) = (HM-U)./(C-U);        
        end
    end
    
    
    
    load('../CuratedCells.mat');
    reg_low = cat(1,Base(1).Cells.regHm);
    reg_high = cat(1,Base(2).Cells.regHm);
    
    fnames = {Base(1).Cells.filename};
    f = @(x)(~isempty(strfind(x,'jbe')));
    isjbe = cellfun(f,fnames);
    
    slopes_low = reg_low(:,2);
    slopes_high = reg_high(:,2);
    
    coeffs = mypolyfit(model_ratio(:,1),model_ratio(:,2),3);
    coeffs_obs = mypolyfit(slopes_low,slopes_high,3);
    
    x = linspace(0,0.4,501);
    modelfit = get_poly(x,coeffs);
    %datafit = get_poly(x,coeffs_obs);
        
    figure(); hold on;
    plot(slopes_low(isjbe),slopes_high(isjbe),'k o','markerfacecolor','m','markersize',10);
    plot(slopes_low(~isjbe),slopes_high(~isjbe),'k ^','markerfacecolor',[0.1,0.1,0.8],'markersize',10);    
    plot(x,modelfit,'r -','linewidth',3);
    legend('Jbe','Lem','Model prediction','location','northwest');
    %plot(x,datafit,'k -','linewidth',3);
    
    set(gca,'xtick',0:0.1:0.4,'ytick',0:0.1:4,'ylim',[-0.1,0.4],'xlim',[-0.1,0.4]);
    
    xlabel('5% density half-matched slope');
    ylabel('24% density half-matched slope')
    
    plot([0,0],[-0.1,0.4],'k --','linewidth',2);
    plot([-0.1,0.4],[0,0],'k --','linewidth',2);
    
    set_plot_params(gcf);
    set(gcf,'position',[100,100,800,600])
    savefig('fig4.fig');
end

function y = get_poly(x,coeffs)
    y = zeros(size(x));
    coeffs = coeffs(end:-1:1);
    
    for j = 1:length(coeffs);
        k = (j-1);
        y = y + coeffs(j) * x.^k;
    end
    
end

function coeffs = mypolyfit(x,y,N);
    
    
    
    pcoeffs = polyfit(x,y,N);
    init = pcoeffs(2:end);
    
    coeffs = fminsearch(@(coeffs)(poly_cost(coeffs,x,y)),init);
    coeffs = [coeffs,0];
end

function cost = poly_cost(coeffs,x,y)
    coeffs = [coeffs,0];
    yhat = get_poly(x,coeffs);
    cost = sum((y-yhat).^2);
end