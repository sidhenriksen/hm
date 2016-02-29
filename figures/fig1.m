function fig1(N)
    % This recreates Figure 1 from Henriksen, Read, & Cumming (2016)
    %
    % This simulation requires BEMtoolbox; available from
    % https://github.com/sidh0/BEMtoolbox
    
    % Run-time options
    run_parallel = 1;
    bootstrap_mode = 1; % will save responses to disk
    
    % Create the BEM object
    bem = BEMunit('x0',0,'y0',0,'silent',1);
    bem.outputNL = @(x)(x.^2);
    bem.Nx = 292; bem.Ny = 292;
    bem.deg_per_pixel=0.02;
    bem = bem.rescale(1.5);
    
    % Create stimulus object
    rds = pairedRDS(); rds.dotsize=4;
    rds.Nx = bem.Nx; rds.Ny = bem.Nx;
    
    % Set stimulus parameters
    densities = [0.05,0.24];
    corrs = [-1,0,1];
    dxs = -24:2:24;
    
    % Pre-allocate space for low and high dot density
    low_tcs = zeros(length(dxs),length(corrs));
    high_tcs = zeros(length(dxs),length(corrs));
    
    if nargin <1
        N = 2e4;
    end
    
    low_rds_ids = zeros(length(dxs),length(corrs));
    high_rds_ids = zeros(length(dxs),length(corrs));
    
    % Loop over correlations and disparities; compute responses to each
    for k = 1:length(corrs);
        rds.correlation = corrs(k);        
        for j = 1:length(dxs);
            rds.dx = dxs(j);
                                    
            rds.density = densities(1);
            if ~N
                bem = bem.load_bootstrap(rds);
            end
            
            low_rds_ids(j,k) = rds.get_identifier();
            
            
            
            low_tcs(j,k) = mean(bem.simulate_spatial(rds,N,bootstrap_mode,run_parallel));
            
            rds.density = densities(2);
            
            if ~N
                bem = bem.load_bootstrap(rds);
            end
            
            high_rds_ids(j,k) = rds.get_identifier();
            
            high_tcs(j,k) = mean(bem.simulate_spatial(rds,N,bootstrap_mode,run_parallel));
            
        end
    end
    
    
    
    
    
    % Finally, plot the data
    fig=figure();
    cols = {'k','b','r'};
    buffer = 0.05;
    for d = 1:length(densities);
        subplot(1,length(densities),d); hold on;
        switch d
            case 1
                tcs = low_tcs;
                ylabel('Normalized response')
            case 2
                tcs = high_tcs;
        end
                                    
        
        markers = {'s','^','o'};
        for k = 1:length(corrs);
            plot(dxs*bem.deg_per_pixel,tcs(:,k)./max(tcs(:)),'-','marker',markers{k},'markersize',8,'markerfacecolor',cols{k},...
                'color',cols{k},'linewidth',3);            
        end
        xlabel('Disparity (deg)')
        set(gca,'ytick',0:0.25:1);
        xlim([min(dxs)*bem.deg_per_pixel - buffer,max(dxs)*bem.deg_per_pixel + buffer]);
    end
    
    set_plot_params(fig);
    b = subplot(1,length(densities),length(densities));
    a=figure();
    hold on;
    plot(dxs*bem.deg_per_pixel,high_tcs(:,2)./max(high_tcs(:)),'color','b','linewidth',2,'markersize',5,'linestyle','-',...
        'markerfacecolor','b','markeredgecolor','b','marker','^');

    
    
    ylim([0.225,0.35]); xlim([min(dxs)*bem.deg_per_pixel - buffer,max(dxs)*bem.deg_per_pixel + buffer]);
    set(gca,'xtick',[-0.5,-0,0.5],'ytick',[0.25,0.35],'fontsize',16);
    
    
    inset_size = 0.15;
    inset_handle = gca;
    new_handle = copyobj(inset_handle,fig);
    close(a);
    ax=get(b,'Position');
    set(new_handle,'Position',[1.1*ax(1)+ax(3)-inset_size 1*ax(2)+ax(4)-inset_size inset_size inset_size])
    
    
    savefig(fig,'fig1.fig');

end