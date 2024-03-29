close all;
clear all;

% Constants
const = [];
const.G = 1;
const.M = 1;
const.c = 10;

% Initial Conditions
v_o = 0.5;
r = 1;
t = 0;
phi = 0;
theta = pi/2;
v_r = 0;
v_phi = v_o / r;
v_theta = 0;
[g__t_t,g__t_r,g__t_theta,g__t_phi,g__r_r,g__r_theta,g__r_phi,g__theta_theta,g__theta_phi,g__phi_phi] = ...
    metric(r,theta,const);
c = const.c;
v_t = sqrt( (c^2 - g__r_r*v_r*v_r - g__phi_phi*v_phi*v_phi - g__theta_theta*v_theta*v_theta)/g__t_t );

% Trails
r_ = r;
t_ = t;
phi_ = phi;
theta_ = theta;

% Plot
ax = axes;
G = const.G;
M = const.M;
r_s = 2*G*M/c^2;
angs = 0:1:360;
hold off;
hTop = plot(ax,r.*cos(phi),r.*sin(phi),'go');
hold on;
xTrail = [r_.*cos(phi_) r_.*cos(fliplr(phi_))];
yTrail = [r_.*sin(phi_) r_.*sin(fliplr(phi_))];
hTop = [hTop patch(xTrail,yTrail,'g')];
hTop(end).EdgeColor = 'none';
hTop(end).FaceColor = 'g';
hTop(end).FaceAlpha = 'flat';
hTop = [hTop patch(r_s.*cosd(angs),r_s.*sind(angs),'y')];
hTop(end).EdgeColor = 'none'; 
axis equal;
axis([-1.1 1.1 -1.1 1.1]);
set(ax,'Position',[0 0 1 1],'XColor','none','YColor','none','Color','k');
set(ax.Parent,'Color',[.5 .5 .5]);

% Initial States
y = [v_r v_t v_phi v_theta r t phi theta];

tau = 0;
h = 0.0001;
iPlot = 0;
while ( tau<20 )
      
    k1 = f(y,const);
    k2 = f(y+0.5*h*k1,const);
    k3 = f(y+0.5*h*k2,const);
    k4 = f(y+h*k3,const);
    y = y + h/6*(k1+2*k2+2*k3+k4);
    r = y(5);
    t = y(6);
    phi = y(7);
    theta = y(8);
    tau = tau + h;

    if ( mod(iPlot,100)==0 )
        r_(end+1) = r;
        t_(end+1) = t;
        phi_(end+1) = phi;
        theta_(end+1) = theta;
    end
    
    if ( mod(iPlot,100)==0 )
        hTop(1).XData = r.*cos(phi);
        hTop(1).YData = r.*sin(phi);
        idx = t_>=t_(end)-1;
        xTrail = [r_(idx).*cos(phi_(idx))];
        yTrail = [r_(idx).*sin(phi_(idx))];
        N = length(xTrail);
        preVertices = [xTrail;yTrail;0*xTrail]';
        preVelocities = diff(preVertices,1);
        preVelocities(end+1,:) = preVelocities(end,:);
        z = ones(N,1)*[0 0 1];
        preNormals = cross(preVelocities,z);
        vertices = [preVertices+preNormals ...
                    preVertices-preNormals];
        verticesX = [vertices(:,1) vertices(:,4)];
        verticesY = [vertices(:,2) vertices(:,5)];
        verticesX = reshape(verticesX',prod(size(verticesX)),1);
        verticesY = reshape(verticesY',prod(size(verticesY)),1);
        alphas = (1 - [0:N-1 N-1:-1:0]/(N-1))';
        hTop(2).Vertices = [verticesX verticesY];
        N = size(hTop(2).Vertices,1);
        hTop(2).Faces = ((1:(N-3)) + [0;2;3;1])';
        hTop(2).FaceVertexAlphaData = ((0:(N-4)))'/(N-4);
        drawnow;
    end
    iPlot = iPlot+1;
    
end

function dY = f(y,const)

    G = const.G;
    M = const.M;
    c = const.c;
    r = y(1);
    r_s = 2*G*M/c^2;

    % y = [v_r, v_t, v_phi, v_theta, r, t, phi, theta]

    v_r = y(1);
    v_t = y(2);
    v_phi = y(3);
    v_theta = y(4);
    r = y(5);
    t = y(6);
    phi = y(7);
    theta = y(8);

    [Gamma_t,Gamma_r,Gamma_phi,Gamma_theta] = ...
        christoffel(r,theta,const);

    v_v = [    v_t*v_t     v_t*v_r     v_t*v_theta     v_t*v_phi ...
               v_r*v_t     v_r*v_r     v_r*v_theta     v_r*v_phi ...
           v_theta*v_t v_theta*v_r v_theta*v_theta v_theta*v_phi ...
             v_phi*v_t   v_phi*v_r   v_phi*v_theta   v_phi*v_phi]';

    dvt_dtau = -Gamma_t*v_v;
    dvr_dtau = -Gamma_r*v_v;
    dvtheta_dtau = -Gamma_theta*v_v;
    dvphi_dtau = -Gamma_phi*v_v;
    
    dY = [dvr_dtau dvt_dtau dvphi_dtau dvtheta_dtau v_r v_t v_phi v_theta];
    
end

function [g__t_t,g__t_r,g__t_theta,g__t_phi, ...
        g__r_r,g__r_theta,g__r_phi,...
        g__theta_theta,g__theta_phi,g__phi_phi] = ...
    metric(r,theta,const);

    % Schwarzschild 
    
    G = const.G;
    M = const.M;
    c = const.c;
    r_s = 2*G*M/c^2;

    g__t_t = (1-r_s/r)*c^2;
    g__t_r = 0;
    g__t_theta = 0;
    g__t_phi = 0;
    g__r_r = -1/(1-r_s/r);
    g__r_theta = 0;
    g__r_phi = 0;
    g__theta_theta = -r^2;
    g__theta_phi = 0;
    g__phi_phi = -r^2*sin(theta)^2;

end

function [Gamma_t,Gamma_r,Gamma_phi,Gamma_theta] = ...
    christoffel(r,theta,const)

    % Schwatzchild
    
    G = const.G;
    M = const.M;
    c = const.c;
    r_s = 2*G*M/c^2;
    
    % t
    Gamma_t__t_t = 0;
    Gamma_t__t_r = r_s/(2*r*(r-r_s));
    Gamma_t__t_theta = 0;
    Gamma_t__t_phi = 0;
    Gamma_t__r_r = 0;
    Gamma_t__r_theta = 0;
    Gamma_t__r_phi = 0;
    Gamma_t__theta_theta = 0;
    Gamma_t__theta_phi = 0;
    Gamma_t__phi_phi = 0;

    % r
    Gamma_r__t_t = c^2*r_s*(r-r_s)/(2*r^3);
    Gamma_r__t_r = 0;
    Gamma_r__t_theta = 0;
    Gamma_r__t_phi = 0;
    Gamma_r__r_r = -r_s/(2*r*(r-r_s));
    Gamma_r__r_theta = 0;
    Gamma_r__r_phi = 0;
    Gamma_r__theta_theta = r_s-r;
    Gamma_r__theta_phi = 0;
    Gamma_r__phi_phi = (r_s-r)*sin(theta)^2;

    % theta
    Gamma_theta__t_t = 0;
    Gamma_theta__t_r = 0;
    Gamma_theta__t_theta = 0;
    Gamma_theta__t_phi = 0;
    Gamma_theta__r_r = 0;
    Gamma_theta__r_theta = 1/r;
    Gamma_theta__r_phi = 0;
    Gamma_theta__theta_theta = 0;
    Gamma_theta__theta_phi = 0;
    Gamma_theta__phi_phi = -sin(theta)*cos(theta);

    % phi
    Gamma_phi__t_t = 0;
    Gamma_phi__t_r = 0;
    Gamma_phi__t_theta = 0;
    Gamma_phi__t_phi = 0;
    Gamma_phi__r_r = 0;
    Gamma_phi__r_theta = 0;
    Gamma_phi__r_phi = 1/r;
    Gamma_phi__theta_theta = 0;
    Gamma_phi__theta_phi = cot(theta);
    Gamma_phi__phi_phi = 0;

    Gamma_t = [Gamma_t__t_t     Gamma_t__t_r     Gamma_t__t_theta     Gamma_t__t_phi ...
               Gamma_t__t_r     Gamma_t__r_r     Gamma_t__r_theta     Gamma_t__r_phi ...
               Gamma_t__t_theta Gamma_t__r_theta Gamma_t__theta_theta Gamma_t__theta_phi ...
               Gamma_t__t_phi   Gamma_t__r_phi   Gamma_t__theta_phi   Gamma_t__phi_phi];

    Gamma_r = [Gamma_r__t_t     Gamma_r__t_r     Gamma_r__t_theta     Gamma_r__t_phi ...
               Gamma_r__t_r     Gamma_r__r_r     Gamma_r__r_theta     Gamma_r__r_phi ...
               Gamma_r__t_theta Gamma_r__r_theta Gamma_r__theta_theta Gamma_r__theta_phi ...
               Gamma_r__t_phi   Gamma_r__r_phi   Gamma_r__theta_phi   Gamma_r__phi_phi];

    Gamma_theta = [Gamma_theta__t_t     Gamma_theta__t_r     Gamma_theta__t_theta     Gamma_theta__t_phi ...
                   Gamma_theta__t_r     Gamma_theta__r_r     Gamma_theta__r_theta     Gamma_theta__r_phi ...
                   Gamma_theta__t_theta Gamma_theta__r_theta Gamma_theta__theta_theta Gamma_theta__theta_phi ...
                   Gamma_theta__t_phi   Gamma_theta__r_phi   Gamma_theta__theta_phi   Gamma_theta__phi_phi];

    Gamma_phi = [Gamma_phi__t_t     Gamma_phi__t_r     Gamma_phi__t_theta     Gamma_phi__t_phi ...
                 Gamma_phi__t_r     Gamma_phi__r_r     Gamma_phi__r_theta     Gamma_phi__r_phi ...
                 Gamma_phi__t_theta Gamma_phi__r_theta Gamma_phi__theta_theta Gamma_phi__theta_phi ...
                 Gamma_phi__t_phi   Gamma_phi__r_phi   Gamma_phi__theta_phi   Gamma_phi__phi_phi];

end