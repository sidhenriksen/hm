function smile(varargin)
    x = linspace(0,1,500);
    y = x;
    
    figure(); hold on;
    % Head
    plot(sin(2*pi*x),cos(2*pi*y),'k -','linewidth',3);

    % Eyes
    plot(0.1*sin(2*pi*x)+0.5,0.1*cos(2*pi*y)+0.25,'k -',0.1*sin(2*pi*x)-0.5,0.1*cos(2*pi*y)+0.25,'k -','linewidth',3)
    % Mouth
    plot(0.5*sin(pi*x + pi/2),0.5*cos(pi*x + pi/2)-0.35, 'k -', 'linewidth',3)
    
end