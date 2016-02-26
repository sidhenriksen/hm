

function estimate_rf_sizes();
   run_viz = 0;
    
    
   data=load('curated_RF_data.mat');
   Xs = data.Xs; Ys = data.Ys; fxs = data.fxs;
    
    R2s = zeros(size(Xs));    
    Ps = cell(size(R2s));
    if run_viz
        figure();
    end
    %%% And now just compute everything in a new loop.. Much easier.                
    for j = 1:size(Xs,2);
        
        for i = 1:size(Xs,1);
            
            x = Xs{i,j};
            y = Ys{i,j};
            params = fit_gaussian(x,y);
            
            
            yhat = get_gaussian(params,x);
            
            R2s(i,j) = corr(y',yhat').^2;
            Ps{i,j} = params;
            
            if run_viz
                subplot(2,2,i); cla;
                hold on;
                plot(x,y,'k o','markerfacecolor','k'); 
                plot(x,yhat,'- r','linewidth',3);
                title(sprintf('R^2=%.2f',R2s(i,j)));

                subplot(2,2,i+2); cla;
                hold on;
                plot(y,yhat,'k o','markerfacecolor','k');
                plot(y,y,'r -','linewidth',3);
            end
        end
        if run_viz
            drawnow;
            pause(3.0);
        end                        
    end
    
    
    
    
    R2sm = mean(R2s);
    theta = 0.8;
    include = R2sm > theta;
    

    Xs_include = Xs(:,include);
    Ys_include = Ys(:,include);
    Ps_include = Ps(:,include);
    isfoveal = abs(fxs(include)) < 5;
    R2sm_include = R2sm(include);

    %sxs = Ps_include(
    all_params1 = cat(1,Ps_include{1,:});
    all_params2 = cat(1,Ps_include{2,:});
    sxs = sqrt(all_params1(:,4).*all_params2(:,4));
    
    figure();
    hold on;
    bins = linspace(0,0.5,31);
    histogram(sxs(isfoveal),bins,'facecolor','blue','facealpha',0.3);
    histogram(sxs(~isfoveal),bins,'facecolor','red','facealpha',0.3);
    legend('Fov','Ecc');
    
    
    k = 0;
    for j = 1:length(Ps_include);
        
        if ~isfoveal(j);           
           if ~mod(k,9)
               figure();
               k = 0;
           end
           k = k+1;
           subplot(3,3,k); hold on;
           X = Xs_include{1,j}; Y = Ys_include{1,j};
           plot(X,Y,'k o','markerfacecolor','k');
           x = linspace(min(X),max(X),301);
           y = get_gaussian(Ps_include{1,j},x);
           plot(x,y,'r -','linewidth',2);
           title(sprintf('s=%.2f, R^2=%.2f,',sxs(j),R2sm_include(j)));
        end
        
    end
    
end
