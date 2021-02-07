function [] = timeGlobeTwin()
    
    c = 1;
    
    WL{1} = [
        0 0 nan
        0 20 nan];
    WL{2} = [
        0 0 nan
        sqrt(10^2-9^2) 10 nan
        0 20 nan];
    
    % Insert Local Times
    maxTime = -inf;
    for idxWL = 1:length(WL)
        
        thisWL = WL{idxWL};
        
        T0_ = 0;
        for idxSections = 2:size(thisWL,1)
            
            WLpointA = thisWL(idxSections-1,:);
            WLpointB = thisWL(idxSections,:);
            
            X = [WLpointA(1) WLpointB(1)];
            T = [WLpointA(2) WLpointB(2)];
            v = diff(X)/diff(T);
            gamma = 1/sqrt(1-v^2/c^2);
            X_ = gamma*(X-v*T);
            T_ = gamma*(T-v*X/c^2);
            T_ = T_ - T_(1) + T0_;
            T0_ = T_(end);
            
            thisWL(idxSections-1,3) = T_(1);
            thisWL(idxSections,3) = T_(2);
            
            maxTime = max(maxTime,WLpointB(2));

        end
        
        WL{idxWL} = thisWL;
        
    end
    
    % Plot World Lines
    plot(maxTime*[-1 0 1],maxTime*[1 0 1],'k','LineWidth',2);
    set(gcf,'MenuBar','none','Toolbar','figure','name','Time Globe','numbertitle','off');
    grid on;
    hTitle = title('0');
    axis equal;
    ax = axis;
    set(gca,'XAxisLocation','origin','YAxisLocation','origin', ...
        'XTick',ceil(ax(1)):floor(ax(2)),'YTick',ceil(ax(3)):floor(ax(4)), ...
        'XTickLabel',{},'YTickLabel',{});
    xlabel('X (Light Years)');
    ylabel('T (Years)');
    hold on;
    hPlots = [];
    set(gca,'ColorOrderIndex',1);
    vFrame = 0;
    for idxWL = 1:length(WL)
        
        thisWL = WL{idxWL};
        
        xl_ = [];
        tl_ = [];
        xm_ = [];
        tm_ = [];
        txm_ = [];
        for idxSections = 2:size(thisWL,1)
            
            WLpointA = thisWL(idxSections-1,:);
            WLpointB = thisWL(idxSections,:);
            
            X = [WLpointA(1) WLpointB(1)];
            T = [WLpointA(2) WLpointB(2)];
            T_ = [WLpointA(3) WLpointB(3)];
            v = vFrame;
            gamma = 1/sqrt(1-v^2/c^2);
            x_ = gamma*(X-v*T);
            t_ = gamma*(T-v*X/c^2);
            xl_ = [xl_ x_];
            tl_ = [tl_ t_];
            tx_ = ceil(T_(1)):floor(T_(2));
            xm_ = [xm_ interp1(T_,x_,tx_)];
            tm_ = [tm_ interp1(T_,t_,tx_)];
            txm_ = [txm_ tx_];
            
        end
        cOi = get(gca,'ColorOrderIndex');
        hPlots(idxWL).lineHandle = plot(xl_,tl_,'LineWidth',2);
        set(gca,'ColorOrderIndex',cOi);
        hPlots(idxWL).markerHandle = plot(xm_,tm_,'o');
        hPlots(idxWL).markerHandle.MarkerFaceColor = hPlots(idxWL).markerHandle.Color;
        hPlots(idxWL).stringHandle = text(xm_,tm_,cellstr(num2str(txm_'))');
        for idx = 1:length(txm_)
            hPlots(idxWL).stringHandle(idx).Position = [xm_(idx) tm_(idx) 0];
            hPlots(idxWL).stringHandle(idx).String = sprintf(' %d',txm_(idx));
        end
        set(hPlots(idxWL).stringHandle,'Color',hPlots(idxWL).markerHandle.Color);
    end
    hold off;
    
    hControl = uicontrol(gcf,'Style','slider','Units','normal','Position',[.1 .02 .8 .05], ...
        'Value',0,'Min',-1,'Max',1,'SliderStep',[.005 .05],'Callback',@UpdateGlobe);

    function [source] = UpdateGlobe(source,eventdata)

        source.Value = round(100*source.Value)/100;
        source.Value = min(max(source.Value,-1),1);
        disp(source.Value);
        v = source.Value;
        hTitle.String = num2str(v);
        
        gamma = 1/sqrt(1-v^2/c^2);
        
        for idxWL = 1:length(WL)
            
            thisWL = WL{idxWL};
            
            xl_ = [];
            tl_ = [];
            xm_ = [];
            tm_ = [];
            txm_ = [];
            for idxSections = 2:size(thisWL,1)
                
                WLpointA = thisWL(idxSections-1,:);
                WLpointB = thisWL(idxSections,:);
                
                X = [WLpointA(1) WLpointB(1)];
                T = [WLpointA(2) WLpointB(2)];
                T_ = [WLpointA(3) WLpointB(3)];
                gamma = 1/sqrt(1-v^2/c^2);
                x_ = gamma*(X-v*T);
                t_ = gamma*(T-v*X/c^2);
                xl_ = [xl_ x_];
                tl_ = [tl_ t_];
                tx_ = ceil(T_(1)):floor(T_(2));
                xm_ = [xm_ interp1(T_,x_,tx_)];
                tm_ = [tm_ interp1(T_,t_,tx_)];
                txm_ = [txm_ tx_];
                
            end
            hPlots(idxWL).lineHandle.XData = xl_;
            hPlots(idxWL).lineHandle.YData = tl_;
            hPlots(idxWL).markerHandle.XData = xm_;
            hPlots(idxWL).markerHandle.YData = tm_;
            for idx = 1:length(txm_)
                hPlots(idxWL).stringHandle(idx).Position = [xm_(idx) tm_(idx) 0];
                hPlots(idxWL).stringHandle(idx).String = sprintf(' %d',txm_(idx));
            end
            
        end
        
        ax = axis;
        set(gca,'XTick',ceil(ax(1)):floor(ax(2)),'YTick',ceil(ax(3)):floor(ax(4)));
    
    end

end
