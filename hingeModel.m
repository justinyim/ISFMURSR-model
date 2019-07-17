% hingeModel.m
% Justin Yim - November 2018

% Robot parameters:
% dimensions: 25mm x 15mm x 70um
% total mass: 50mg
% front foot dimensions: 10mg 3-5mm
% cantilever maximum deflection at DC: 1.67mm
% free-free resontant frequency: 199.9Hz (FEM)
% maximum speed freqency: 200Hz
% maximum speed: 8.7cm/s
% running displacement at 200Hz (resonance): 3mm

% I: Analytical Dynamics ==================================================
% Variables ---------------------------------------------------------------
syms t % time

% States
syms x y th1 th2
syms vx vy w1 w2
syms ax ay alpha1 alpha2
p = {x;y;th1;th2}; % front link CG x (m), CGy (m), angle (rad), rear link angle (rad)
v = {vx;vy;w1;w2}; % time derivatives of p
a = {ax;ay;alpha1;alpha2}; % time deriviatives of v
X = [x, y, th1, th2, vx, vy, w1, w2]; % state vector

syms xd yd th1d th2d
syms vxd vyd w1d w2d
pd = [xd;yd;th1d;th2d]; % temporary time derivatives of p
vd = [vxd;vyd;w1d;w2d]; % temporary time derivatives of v

% Parameters 
syms l10x l10y l12 l21 l20x l20y % kinematic parameters
syms m1 m2 I1 I2 % mass parameters

% Forces
syms F1x F1y F2x F2y tau g % force parameters

% Kinematics --------------------------------------------------------------
% Robot pose: locations of pin joint and foot points
cg1 = [x;y]; % front link CG point
pin = cg1 + [-l12*cos(th1); % pin joint point
             -l12*sin(th1)];
cg2 = pin + [-l21*cos(th2); % rear link CG point
             -l21*sin(th2)];
foot1 = cg1 + [l10x*cos(th1) + l10y*sin(th1); % front foot point
               l10x*sin(th1) - l10y*cos(th1)];
foot2 = cg2 + [-l20x*cos(th2) + l20y*sin(th2); % rear foot point
               -l20x*sin(th2) - l20y*cos(th2)];

% Velocities
v2 = diffT(cg2,p,1); % rear link CG velocity
v2 = subs(v2,pd,v);

vfoot1 = subs(diffT(foot1,p,1),pd,v); % front foot velocity
vfoot2 = subs(diffT(foot2,p,1),pd,v); % rear foot velocity

% Accelerations
a2 = diffT(v2,[p;v],1); % rear link CG acceleration
a2 = subs(a2,[pd;vd],[v;a]);

coeffs = vpa(zeros(6,length(a)));
const = [ax; ay; alpha1; a2; alpha2];
for ii = 1:length(a)
    coeffs(:,ii) = diff([ax; ay; alpha1; a2; alpha2],a(ii));
    const = const - coeffs(:,ii)*a(ii);
end

% Dynamics ----------------------------------------------------------------
% Force & torque balance matrices: M*x = F
% where x = [ax1; ay1; alpha1; ax2; ay2; alpha2]
M = diag([m1; m1; I1; m2; m2; I2]);
F = [F1x;
    F1y - m1*g;
    tau - F1x*(foot1(2)-cg1(2)) + F1y*(foot1(1)-cg1(1));
    F2x;
    F2y - m2*g;
    -tau - F2x*(foot2(2)-cg2(2)) + F2y*(foot2(1)-cg2(1))] - M*const;

% Transforming into minimal coordinates M*z = F
% where z = [ax1; ay1; alpha1; alpha2; Fpx; Fpy]
M = M*coeffs;

% Pin joint x and y reaction forces respectively
M = [M, ...
    [1, 0;
     0, 1;
     -(pin(2)-cg1(2)), (pin(1)-cg1(1));
     -1, 0;
     0, -1;
     (pin(2)-cg2(2)), -(pin(1)-cg2(1))]];

%% II: Simulate with numerical parameters =================================

% Front foot is 1.5mm back from the tip.

% Positions of points relative to the hinge:
%   CG should be an additional 1/3 of the way towards the front foot
%   (20mg * 6.25mm + 10mg * 10mm)/30mg = 7.5mm
% unimorph CG 6.25mm
% front body CG 7.5mm (1.25mm further than unimorph)
% foot 10mm (2.5mm further than front body CG)


% Numerical Parameters ----------------------------------------------------
% Inertia and dimensions
ms = 100;                   % mass scaling to avoid numerical issues
g_ = 9.81;                  % Acceleration due to gravity (m/s^2)
m1_ = 30E-6*ms;             % front link mass (kg)
m2_ = 20E-6*ms;             % rear link mass (kg)
l12_ = 0.75E-2;             % front link to pin length (m)
l10x_ = 0.25E-2;            % front link to foot x distance (m)
l10y_ = 0.35E-2;            % front link to foot y distance (m)
l21_ = 0.625E-2;            % rear link to pin length (m)
l20x_ = 0.625E-2;           % rear link to foot x distance (m)
l20y_ = 0.05E-2;            % rear link to foot y distance (m)

% Ground parameters
mu_ = 0.3;          % coefficient of friction (ratio)
ky = 8.5E0;         % ground stiffness (N/m)
dy = 5E-3;          % ground damping [N/(m/s)]
dx = 1E0;           % x ground damping for smoothing [N/(m/s)]

% Actuators parameters
freqs = logspace(2.02,2.5,8);   % Acutator frequencies to simulate (Hz)

ang = 0.2;                      % resting angle (rad)
taumag = 5E-6;          % actuator torque amplitude (N m)
kspring = 2.5E-4;       % spring stiffness (N m/rad)
dspring = 1E-8;         % spring damping [N m/(rad/s)]

% Parameters derived from those specified above
% front link moment of inertia (kg m^2)
I1_ = 1/12*m2_*(2*0.625E-2)^2 + m2_*0.125E-2^2 + 10E-6*ms*0.25E-2^2;  
I2_ = 1/12*m2_*(2*0.625E-2)^2;  % [-] rear link moment of inertia (kg m^2)

params1 =  [m1, m2, I1, I2, g, l12, l10x, l10y, l21, l20x, l20y];
params1_ = [m1_,m2_,I1_,I2_,g_,l12_,l10x_,l10y_,l21_,l20x_,l20y_];


% Numerical expressions ---------------------------------------------------
foot1_ = subs(foot1,params1,params1_);
pin_ = subs(pin,params1,params1_);
cg2_ = subs(cg2,params1,params1_);
foot2_ = subs(foot2,params1,params1_);

vfoot1_ = subs(vfoot1,params1,params1_);
vcg2_ = subs(v2,params1,params1_);
vfoot2_ = subs(vfoot2,params1,params1_);

% Functions
ffoot1 = dynamicsFunction(foot1_,X);
fpin = dynamicsFunction(pin_,X);
fcg2 = dynamicsFunction(cg2_,X);
ffoot2 = dynamicsFunction(foot2_,X);
ffoot1y = dynamicsFunction(foot1_(2),X);
ffoot2y = dynamicsFunction(foot2_(2),X);

fvcg2 = dynamicsFunction(vcg2_,X);
fvfoot1y = dynamicsFunction(vfoot1_(2),X);
fvfoot2y = dynamicsFunction(vfoot2_(2),X);
fvfoot1x = dynamicsFunction(vfoot1_(1),X);
fvfoot2x = dynamicsFunction(vfoot2_(1),X);

% Forces
min_s = @(i1,i2) heaviside(i1-i2)*i2 + heaviside(i2-i1)*i1;

F1y_ = ms*heaviside(-foot1_(2))*heaviside(-ky*foot1_(2) - dy*vfoot1_(2))*(-ky*foot1_(2) - dy*vfoot1_(2));
F1x_ = ms*-sign(vfoot1_(1))*(min_s(abs(dx*vfoot1_(1)),abs(mu_*F1y_/ms)));
F2y_ = ms*heaviside(-foot2_(2))*heaviside(-ky*foot2_(2) - dy*vfoot2_(2))*(-ky*foot2_(2) - dy*vfoot2_(2));
F2x_ = ms*-sign(vfoot2_(1))*(min_s(abs(dx*vfoot2_(1)),abs(mu_*F2y_/ms)));

fF1x = dynamicsFunction(F1x_,X);
fF1y = dynamicsFunction(F1y_,X);
fF2x = dynamicsFunction(F2x_,X);
fF2y = dynamicsFunction(F2y_,X);

trajectories = cell(length(freqs),2);

%% III: Simulation ========================================================
% Solve for initial conditions
objectiveF = matlabFunction([foot2_(2)-1E-6; foot1_(2)-1E-6; th2 - th1 - ang; foot1_(1)], ...
    'Vars',{x,y,th1,th2});
% initial pose
pos0 = lsqnonlin(@(x)objectiveF(x(1),x(2),x(3),x(4)), [0;0;0;0]);

for ii = 1:length(freqs)
    fprintf([num2str(ii),': ']);
    tic

    % Motor torque
    tau_ = ms*(kspring*(th2-th1-ang) + taumag*(cos(2*pi*freqs(ii)*t)) - dspring*(w1-w2));

    params =  [params1, F1x, F1y, F2x, F2y, tau];
    params_ = [params1_,F1x_,F1y_,F2x_,F2y_,tau_];

    % Numerical dynamics and simulation ---------------------------------------
    M_ = subs(M,params,params_); % numerical M matrix
    F_ = subs(F,params,params_); % numerical F vector

    fM = dynamicsFunction(M_,X); % M matrix function
    fF = dynamicsFunction(F_,X); % F vector function
    fv = @(t,X) X(5:8); % function for velocities

    fprintf(['functions ',num2str(toc,3),' s. '])
    
    X0 = [pos0 + [0;0.001;0;0];0;0;0;0];    % initial state 1mm above the ground
    options = odeset('reltol',1e-4);        % ODE solver options
    ts = [0,1];                             % simulation duration (s)
    
    % Run the simulation
    tic
    [tout,xout] = ode45(@(t,X)MFdynamicsOrdered(t,X,fM,fF,fv), ts, X0, options);
    disp(['Simulated ', num2str(length(tout)),' steps in ',num2str(toc),' s'])

    trajectories(ii,:) = {tout, xout};

end

%% IV: Post-processing ----------------------------------------------------
% Plotting parameters
framerate = 5000;   % Framerate (hz)
xScale = 1000;      % Display scale (ratio)
Fscale = 1E-2;      % Force to distance display scale (ratio)
videoName = '200Hz_5000FPS.avi';     % Save video to videoName if nonemtpy.

toPlot = 5;%1:length(freqs); % which simulation index to plot

for n = length(freqs):-1:1
    tout = trajectories{n,1};
    xout = trajectories{n,2};
    
    % Initialize arrays of uninterpolated raw timing output values
    cgs{n} = zeros(length(tout),2); % total system CG location (m)
    cg2s{n} = zeros(length(tout),2); % back body CG location (m)
    v2s{n} = zeros(length(tout),2); % back body CG velocity (m/s)
    fys{n} = zeros(length(tout),2); % front and back foot heights (m)
    cgvs{n} = zeros(length(tout),2); % total system CG velocity (m/s)
    vfys{n} = zeros(length(tout),2); % front and back foot y velocities (m/s)
    vfxs{n} = zeros(length(tout),2); % front and back foot x velocities (m/s)
    KEs{n} = zeros(length(tout),2); % front and back body kinetic energy (J)
    PEs{n} = zeros(length(tout),2); % front and back body potential energy (J)
    SPEs{n} = 0.5*kspring*(xout(:,3)-xout(:,4)-ang).^2; % hinge spring potential energy (J)

    F1sx{n} = zeros(length(tout),1); % Front foot x force (N)
    F1sy{n} = zeros(length(tout),1); % Front foot y force (N)
    F2sx{n} = zeros(length(tout),1); % Back foot x force (N)
    F2sy{n} = zeros(length(tout),1); % Back foot y force (N)
    
    ptsx{n} = zeros(length(tout),5); % [front foot, front CG, joint, back CG, back foot] for plotting
    ptsy{n} = zeros(length(tout),5); % [front foot, front CG, joint, back CG, back foot] for plotting

    % Populate uninterpolated raw timing output values
    for ii = 1:length(tout)
        cg2s{n}(ii,:) = fcg2(tout(ii),xout(ii,:))';
        v2s{n}(ii,:) = fvcg2(tout(ii),xout(ii,:))';
        cgs{n}(ii,:) = (m1_*xout(ii,1:2) + m2_*cg2s{n}(ii,:))/(m1_+m2_);
        fys{n}(ii,:) = [ffoot1y(tout(ii),xout(ii,:)), ffoot2y(tout(ii),xout(ii,:))];
        vfys{n}(ii,:) = [fvfoot1y(tout(ii),xout(ii,:)), fvfoot2y(tout(ii),xout(ii,:))];
        vfxs{n}(ii,:) = [fvfoot1x(tout(ii),xout(ii,:)), fvfoot2x(tout(ii),xout(ii,:))];
        KEs{n}(ii,:) = [0.5*m1_*sum(xout(ii,5:6).^2) + 0.5*I1_*xout(ii,7).^2, ...
            0.5*m2_*sum(v2s{n}(ii,:).^2) + 0.5*I2_*xout(ii,8).^2];
        PEs{n}(ii,:) = [m1_*g_*xout(ii,2), m2_*g_*cg2s{n}(ii,2)];

        pt = [ffoot1(tout(ii),xout(ii,:)), xout(ii,1:2)', fpin(tout(ii),xout(ii,:)), ...
            fcg2(tout(ii),xout(ii,:)), ffoot2(tout(ii),xout(ii,:))];
        ptsx{n}(ii,:) = pt(1,:);
        ptsy{n}(ii,:) = pt(2,:);
        F1sx{n}(ii) = fF1x(tout(ii),xout(ii,:));
        F1sy{n}(ii) = fF1y(tout(ii),xout(ii,:));
        F2sx{n}(ii) = fF2x(tout(ii),xout(ii,:));
        F2sy{n}(ii) = fF2y(tout(ii),xout(ii,:));
    end
    cgvs{n} = (m2_*v2s{n} + m1_*xout(:,[5,6]))/(m1_+m2_);

    % Interpolated timings at constant framerate for plotting
    tint{n} = [tout(1):(1/framerate):tout(end), tout(end)];
    xint{n} = interp1(tout,xout,tint{n});

    ptisx{n} = zeros(length(tint{n}),5);
    ptisy{n} = zeros(length(tint{n}),5);
    cgis{n} = zeros(length(tint{n}),2);
    F1isx{n} = zeros(length(tint{n}),1);
    F1isy{n} = zeros(length(tint{n}),1);
    F2isx{n} = zeros(length(tint{n}),1);
    F2isy{n} = zeros(length(tint{n}),1);

    % Populate interpolated timings
    for ii = 1:length(tint{n})
        pt = [ffoot1(tint{n}(ii),xint{n}(ii,:)), xint{n}(ii,1:2)', fpin(tint{n}(ii),xint{n}(ii,:)), ...
            fcg2(tint{n}(ii),xint{n}(ii,:)), ffoot2(tint{n}(ii),xint{n}(ii,:))];
        ptisx{n}(ii,:) = pt(1,:);
        ptisy{n}(ii,:) = pt(2,:);
        cgis{n}(ii,:) = sum(pt(:,[2,4]).*[m1_,m2_; m1_,m2_],2)/(m1_+m2_);

        F1isx{n}(ii) = fF1x(tint{n}(ii),xint{n}(ii,:));
        F1isy{n}(ii) = fF1y(tint{n}(ii),xint{n}(ii,:));
        F2isx{n}(ii) = fF2x(tint{n}(ii),xint{n}(ii,:));
        F2isy{n}(ii) = fF2y(tint{n}(ii),xint{n}(ii,:));
    end

    % Gestures
    t_window = (tout > 0.2); % Times at which to extract gestures
    contracting = xout(:,8) > xout(:,7); % Contraction-expansion phase
    front_touch = F1sy{n} > 0; % Front foot on the ground
    back_touch = F2sy{n} > 0; % Back foot on the ground

    dt = tout(2:end) - tout(1:(end-1)); % time step (s)

    % Gesture arrays
    aci{n} = contracting & ~front_touch & ~back_touch & t_window; % Aerial contracting
    aei{n} = ~contracting & ~front_touch & ~back_touch & t_window; % Aerial expanding
    fci{n} = contracting & front_touch & ~back_touch & t_window; % Front-touching contracting
    fei{n} = ~contracting & front_touch & ~back_touch & t_window; % Front-touching expanding
    rci{n} = contracting & ~front_touch & back_touch & t_window; % Rear-touching contracting
    rei{n} = ~contracting & ~front_touch & back_touch & t_window; % Rear-touching expanding
    bci{n} = contracting & front_touch & back_touch & t_window; % Both-touching contracting
    bei{n} = ~contracting & front_touch & back_touch & t_window; % Both-touching expanding
    
    % Direction of foot movement
    fr{n} = vfxs{n}(:,1) < 0 & front_touch & ~back_touch & t_window; % Front foot drags back
    ff{n} = vfxs{n}(:,1) > 0 & front_touch & ~back_touch & t_window; % Front foot drags forward
    br{n} = vfxs{n}(:,2) < 0 & ~front_touch & back_touch & t_window; % Rear foot drags back
    bf{n} = vfxs{n}(:,2) > 0 & ~front_touch & back_touch & t_window; % Rear foot drags forward

    % Count gestures and transitions between gestures
    gestureVector{n} = [aci{n}, aei{n}, fci{n}, fei{n}, rci{n}, rei{n}, bci{n}, bei{n}];
    gestureTransitions{n} = find(gestureVector{n}(1,:));
    gestureMatrix{n} = zeros(8,8);
    for ii = 2:length(tout)
        prev = find(gestureVector{n}(ii-1,:));
        next = find(gestureVector{n}(ii,:));
        if prev ~= next
            gestureMatrix{n}(next,prev) = gestureMatrix{n}(next,prev)+1;
            gestureTransitions{n} = [gestureTransitions{n}; next];
        end
    end
    
    % Duty cycles for each gesture
    ac{n} = sum(dt.*aci{n}(1:(end-1)));
    ae{n} = sum(dt.*aei{n}(1:(end-1)));
    fc{n} = sum(dt.*fci{n}(1:(end-1)));
    fe{n} = sum(dt.*fei{n}(1:(end-1)));
    rc{n} = sum(dt.*rci{n}(1:(end-1)));
    re{n} = sum(dt.*rei{n}(1:(end-1)));
    bc{n} = sum(dt.*bci{n}(1:(end-1)));
    be{n} = sum(dt.*bei{n}(1:(end-1)));
    
    
    % Additional processing -----------------------------------------------
    % The following are additional values calculated from the simulation
    % for interest.
    
    % Impulses in stance phase
    td1{n} = find(front_touch(2:end) & ~front_touch(1:(end-1)));
    to1{n} = find(~front_touch(2:end) & front_touch(1:(end-1)));
    td2{n} = find(back_touch(2:end)& ~back_touch(1:(end-1)));
    to2{n} = find(~back_touch(2:end)& back_touch(1:(end-1)));
    
    totalImpulse1x = cumtrapz(tout, F1sx{n});
    totalImpulse1y = cumtrapz(tout, F1sy{n});
    totalImpulse2x = cumtrapz(tout, F2sx{n});
    totalImpulse2y = cumtrapz(tout, F2sy{n});
    
    impulse1x{n} = zeros(length(td1{n})-1,1);
    impulse1y{n} = impulse1x{n};
    force1x{n} = impulse1x{n};
    force1y{n} = impulse1x{n};
    velocity1y{n} = impulse1x{n};
    stance1{n} = impulse1x{n};
    KE{n} = impulse1x{n};
    PE{n} = impulse1x{n};
    for ii = 1:(length(td1{n})-1)
        impulse1x{n}(ii) = totalImpulse1x(to1{n}(ii)) - totalImpulse1x(td1{n}(ii)); % stance impulse
        impulse1y{n}(ii) = totalImpulse1y(to1{n}(ii)) - totalImpulse1y(td1{n}(ii));
        stance1{n}(ii) = tout(to1{n}(ii))-tout(td1{n}(ii));
        force1x{n}(ii) = impulse1x{n}(ii)/(tout(to1{n}(ii))-tout(td1{n}(ii))); % average force
        force1y{n}(ii) = impulse1y{n}(ii)/(tout(to1{n}(ii))-tout(td1{n}(ii)));
        velocity1y{n}(ii) = vfys{n}(td1{n}(ii),1); % touchdown velocity
    end
    impulse2x{n} = zeros(length(td2{n})-1,1);
    impulse2y{n} = impulse2x{n};
    force2x{n} = impulse2x{n};
    force2y{n} = impulse2x{n};
    velocity2y{n} = impulse2x{n};
    stance2{n} = impulse2x{n};
    for ii = 1:(length(td2{n})-1)
        impulse2x{n}(ii) = totalImpulse2x(to2{n}(ii)) - totalImpulse2x(td2{n}(ii));
        impulse2y{n}(ii) = totalImpulse2y(to2{n}(ii)) - totalImpulse2y(td2{n}(ii));
        stance2{n}(ii) = tout(to2{n}(ii))-tout(td2{n}(ii));
        force2x{n}(ii) = impulse2x{n}(ii)/(tout(to2{n}(ii))-tout(td2{n}(ii)));
        force2y{n}(ii) = impulse2y{n}(ii)/(tout(to2{n}(ii))-tout(td2{n}(ii)));
        velocity2y{n}(ii) = vfys{n}(td2{n}(ii),2);
    end
end
   
%% V: Plotting ============================================================

figure(1)
clf
hold all
for n = toPlot
    stem((n-1)/(2*length(freqs))+(1:8),100*[ac{n},ae{n},fc{n},fe{n},rc{n},re{n},bc{n},be{n}]/ ...
        sum(dt.*t_window(1:(end-1))))
end
hold off
set(gca,'xticklabel',{'','ac','ae','fc','fe','rc','re','bc','be',''})
set(gca,'xlim',[0,9],'ylim',[0,50])
ylabel('Duty cycle')
xlabel('Gesture')
set(gca,'fontsize',14)


speed = zeros(size(freqs));
for ii = 1:size(trajectories,1)
    startInd = find(trajectories{ii,1} > 0.2, 1, 'first');
    speed(ii) = (trajectories{ii,2}(end,1) - trajectories{ii,2}(startInd,1)) / ...
        (trajectories{ii,1}(end,1) - trajectories{ii,1}(startInd,1));
end

figure(2)
clf
plot(freqs,speed,'.-')
hold on
plot(freqs([1,end]),[0,0],'k')
hold off
xlabel('Frequency (Hz)')
ylabel('Speed (m/s)')
set(gca,'fontsize',14)

for n = toPlot
    tout = trajectories{n,1};
    xout = trajectories{n,2};

    figure(100*n+1)
    clf
    subplot(2,2,1)
    plot(tout,xout(:,1:2))
    hold on
    plot(tout,cgs{n}(:,1),'k')
    plot(tout,cgs{n}(:,2),':k')
    hold off
    legend('x','y','cgx','cgy')
    ylabel('Distance (m)')
    set(gca,'fontsize',14)
    subplot(2,2,3)
    plot(tout,xout(:,3:4))
    legend('th1','th2')
    ylabel('Angle (rad)')
    xlabel('Time (s)')
    set(gca,'fontsize',14)
    subplot(2,2,2)
    plot(tout,xout(:,5:6))
    ylabel('Velocity (m/s)')
    set(gca,'fontsize',14)
    subplot(2,2,4)
    plot(tout,xout(:,7:8))
    ylabel('Angular velocity (m/s)')
    xlabel('Time (s)')
    set(gca,'fontsize',14)

    figure(100*n+2)
    clf
    h4_1 = subplot(2,1,1);
    plot(tout,F1sx{n},tout,F2sx{n})
    ylabel('x Force (N)')
    set(gca,'fontsize',14)
    h4_2 = subplot(2,1,2);
    plot(tout,F1sy{n},tout,F2sy{n})
    ylabel('y Force (N)')
    xlabel('Time (s)')
    set(gca,'fontsize',14)
    linkaxes([h4_1, h4_2],'x')
    
    figure(100*n+3)
    clf
    plot(tout,cgvs{n})
    xlabel('Time (s)')
    ylabel('Total CG velocity (m/s)')
    set(gca,'fontsize',14)
    
end
% Animate

% Write a video file
saveVideo = ~isempty(videoName);
if saveVideo
    vid = VideoWriter(videoName);
    vid.FrameRate = 30;
    vid.Quality = 98;
    open(vid);
end

n = toPlot(end);

tout = trajectories{n,1};
xout = trajectories{n,2};
f_animate = figure(101);
clf
a = axes;
hold on
h = plot(0,0,'.-k','markersize',10);
hF1 = quiver(0,0,0,0,'r','linewidth',1);
hF2 = quiver(0,0,0,0,'r','linewidth',1);
set(hF1,'autoscale','off','showarrowhead','off');
set(hF2,'autoscale','off','showarrowhead','off');
h_gnd = plot([-1,1],[0,0],'k');
hold off
set(gca,'ylim',[-0.01,0.05])
axis equal
xlabel('x (mm)')
ylabel('y (mm)')
h_title = title('Initializing ...');
set(gca,'fontsize',14)
for ii = 1:2500%find(tint{n}>0.2,1,'first'):min(length(tint{n}),find(tint{n}>0.2,1,'first')+300)%1:2500%
    tic
    set(h,'xdata',xScale*ptisx{n}(ii,:),'ydata',xScale*ptisy{n}(ii,:))
    set(hF1,'xdata',xScale*ptisx{n}(ii,1),'ydata',xScale*ptisy{n}(ii,1),...
        'udata',F1isx{n}(ii)*Fscale*xScale,'vdata',F1isy{n}(ii)*Fscale*xScale)
    set(hF2,'xdata',xScale*ptisx{n}(ii,5),'ydata',xScale*ptisy{n}(ii,5),...
        'udata',F2isx{n}(ii)*Fscale*xScale,'vdata',F2isy{n}(ii)*Fscale*xScale)
    set(a,'xlim',xScale*(cgis{n}(ii,1)+[2*(-l20x_-l21_),2*(l10x_+l12_)]), ...
        'ylim',xScale*[-0.01,l10x_+l12_+l20x_+l21_])
    set(h_gnd,'xdata',xScale*(cgis{n}(ii,1)+[2*(-l20x_-l21_),2*(l10x_+l12_)]))
    set(h_title,'string',[num2str(floor(tint{n}(ii)*1000),3),' ms'])
    drawnow
    
    if saveVideo
        frame = getframe(f_animate);
        frame = frame.cdata;
        writeVideo(vid, frame);
    else
        pause(toc-(1/30))
    end
    
end

if saveVideo
    close(vid);
end