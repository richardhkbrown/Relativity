close all;
clear all;

const = [];
const.G = 1;
const.M = 1;
const.m = 1;
const.L = 0.5;
const.E = -1;
const.c = 10;
const.rDir = -1;

rPrev_newton = inf;
rSlope_newton = -inf;
rPrev_gr = inf;
rSlope_gr = -inf;

r = 0.001:.001:3;
Vref_newton = p_newton(r,const);
Vref_gr = p_gr(r,const);

idx= r>0.4;
r0 = interp1(Vref_newton(idx),r(idx),const.E);
r0 = r0 - 1e-6;
phi0 = 0;
t0 = 0;
y_newton = [r0 phi0 t0];
y_gr = [r0 phi0 t0];

r_newton = r0;
phi_newton = phi0;
r_gr = r0;
phi_gr = phi0;

figure(1);
ax1 = axes;
hE = plot(ax1,r_newton,const.E,'go',r_gr,const.E,'mo','MarkerSize',20,'LineWidth',2);
hold on;
set(gca,'ColorOrderIndex',1);
plot(r,Vref_newton,'g',r,Vref_gr,'m','LineWidth',2);
hold off;
axis([0 1.5 -3 0]);
set(figure(1),'Position',[100 100 1280 720]);
set(ax1,'Position',[0.563 0 1-.563 1],'XColor','none','YColor','none','Color','k');

ax2 = axes;
G = const.G;
M = const.M;
c = const.c;
r_s = 2*G*M/c^2;
angs = 0:1:360;
hTop = plot(ax2,r_newton.*cosd(phi_newton),r_newton.*sin(phi_newton),'g', ...
    r_gr.*cosd(phi_gr),r_gr.*sin(phi_gr),'m','LineWidth',2);
hold on;
hTop = [hTop;plot(ax2,r_newton.*cosd(phi_newton),r_newton.*sin(phi_newton),'go', ...
    r_gr.*cosd(phi_gr),r_gr.*sin(phi_gr),'mo', ...
    r_s*cosd(angs),r_s*sind(angs),'y','LineWidth',2,'MarkerSize',15)];
axis equal;
axis([-1 1 -1 1]);
set(ax2,'Position',[0 0 .563 1],'XColor','none','YColor','none','Color','k');

vidObj = VideoWriter('eom.mp4','MPEG-4');
open(vidObj);

tau = 0;
h = 0.0001;
iPlot = 0;
while ( tau<12.1 )
    
    y = y_newton;
    
    const.rDir = 1;
    k1 = f_newton(y,const);
    k2 = f_newton(y+0.5*h*k1,const);
    k3 = f_newton(y+0.5*h*k2,const);
    k4 = f_newton(y+h*k3,const);
    y1 = y + h/6*(k1+2*k2+2*k3+k4);
    const.rDir = -1;
    k1 = f_newton(y,const);
    k2 = f_newton(y+0.5*h*k1,const);
    k3 = f_newton(y+0.5*h*k2,const);
    k4 = f_newton(y+h*k3,const);
    y2 = y + h/6*(k1+2*k2+2*k3+k4);

    if ( any(imag(y1)) && ~any(imag(y2)) )
        y_newton = y2;
    elseif ( any(imag(y2)) && ~any(imag(y1)) )
        y_newton = y1;
    elseif ( isinf(rPrev_newton) || isinf(rPrev_newton) )
        y_newton = y2;
    else
        rSlope1 = y1(1)-y_newton(1);
        rSlope2 = y2(1)-y_newton(1);
        if ( abs(rSlope1-rSlope_newton) < abs(rSlope2-rSlope_newton) )
            y_newton = y1;
        else
            y_newton = y2;
        end
    end
    r = y_newton(1);
    rSlope_newton = r-rPrev_newton;
    rPrev_newton = r;
    phi = y_newton(2);
    
    if ( mod(iPlot,100)==0 )
        r_newton(end+1) = r;
        phi_newton(end+1) = phi;
    end
    
    y = y_gr;
    
    const.rDir = 1;
    k1 = f_gr(y,const);
    k2 = f_gr(y+0.5*h*k1,const);
    k3 = f_gr(y+0.5*h*k2,const);
    k4 = f_gr(y+h*k3,const);
    y1 = y + h/6*(k1+2*k2+2*k3+k4);
    const.rDir = -1;
    k1 = f_gr(y,const);
    k2 = f_gr(y+0.5*h*k1,const);
    k3 = f_gr(y+0.5*h*k2,const);
    k4 = f_gr(y+h*k3,const);
    y2 = y + h/6*(k1+2*k2+2*k3+k4);

    if ( any(imag(y1)) && ~any(imag(y2)) )
        y_gr = y2;
    elseif ( any(imag(y2)) && ~any(imag(y1)) )
        y_gr = y1;
    elseif ( isinf(rPrev_newton) || isinf(rPrev_newton) )
        y_gr = y2;
    else
        rSlope1 = y1(1)-y_gr(1);
        rSlope2 = y2(1)-y_gr(1);
        if ( abs(rSlope1-rSlope_gr) < abs(rSlope2-rSlope_gr) )
            y_gr = y1;
        else
            y_gr = y2;
        end
    end
    r = y_gr(1);
    rSlope_gr = r-rPrev_gr;
    rPrev_gr = r;
    phi = y_gr(2);
    
    tau = tau + h;
    
    if ( mod(iPlot,100)==0 )
        r_gr(end+1) = r;
        phi_gr(end+1) = phi;
    end
    
    if ( mod(iPlot,100)==0 )
        hTop(1).XData = r_newton.*cos(phi_newton);
        hTop(1).YData = r_newton.*sin(phi_newton);
        hTop(2).XData = r_gr.*cos(phi_gr);
        hTop(2).YData = r_gr.*sin(phi_gr);
        hTop(3).XData = r_newton(end).*cos(phi_newton(end));
        hTop(3).YData = r_newton(end).*sin(phi_newton(end));
        hTop(4).XData = r_gr(end).*cos(phi_gr(end));
        hTop(4).YData = r_gr(end).*sin(phi_gr(end));
        m = const.m;
        E = const.E;
        hE(1).XData = r_newton(end);
        dY = f_newton(y_newton,const);
        hE(1).YData = E-0.5*m*dY(1)^2;
        hE(2).XData = r_gr(end);
        dY = f_gr(y_gr,const);
        hE(2).YData = E-0.5*m*dY(1)^2;
        drawnow;
        writeVideo(vidObj, getframe(gcf));
    end
    iPlot = iPlot+1;
    
end

close(vidObj);

function Vref = p_newton(r,const)
    
    m = const.m;
    G = const.G;
    M = const.M;
    L = const.L;

    Vref = m*(-G*M./r+L^2./(2*m^2*r.^2));

end

function Vref = p_gr(r,const)
    
    m = const.m;
    G = const.G;
    M = const.M;
    L = const.L;
    c = const.c;

    Vref = m*(-G*M./r+L^2./(2*m^2*r.^2)-G*M/c^2*L.^2./(m^2*r.^3));

end

function dY = f_newton(y,const)

    m = const.m;
    G = const.G;
    M = const.M;
    L = const.L;
    c = const.c;
    E = const.E;
    r = y(1);
    Vref = p_newton(r,const);
    r_s = 2*G*M/c^2;
    
    dr_dtau = const.rDir*sqrt(2*(E-Vref)/m);
    dphi_dtau = L/(m*r^2);
    dt_dtau = E/m/(1-r_s/r);
    
    dY = [dr_dtau dphi_dtau dt_dtau];
    
end

function dY = f_gr(y,const)

    m = const.m;
    G = const.G;
    M = const.M;
    L = const.L;
    c = const.c;
    E = const.E;
    r = y(1);
    Vref = p_gr(r,const);
    r_s = 2*G*M/c^2;
    
    dr_dtau = const.rDir*sqrt(2*(E-Vref)/m);
    dphi_dtau = L/(m*r^2);
    dt_dtau = E/m/(1-r_s/r);
    
    dY = [dr_dtau dphi_dtau dt_dtau];
    
end