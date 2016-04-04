clear all
close all;clc

lambda = 1; %633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);

k_air = 2*pi/lambda; % propagation constant of air
k_silver = omega * sqrt(mu_0*epsilon_0*eps_silver); % propagation constant of silver
%%
len = 1e4; % Vector Length
%%
kxx = linspace(-1e1*k_air,1e4*k_air,len);
kxy = linspace(-1e4*k_air,1e1*k_air,len);

% Find the branch point location on the kx x-axis
dif = abs(kxx - k_air);
bp1_loc = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

x = linspace(1e-2*lambda,1e4*lambda,len);
H_c1 = zeros(length(x),1); % Initialize the magnetic field vector
H_c2 = zeros(length(x),1); % Initialize the magnetic field vector
H_c3 = zeros(length(x),1); % Initialize the magnetic field vector
%% Define Contour
%                 C2
%               --->--
%               | k1 |
%               |    |
%             C1|    |C3
%               |    |
%       ------> |    V <--------
%       bottom  |    | top sheet
%               ^    |
% Define real and imaginary start points for c1
%
c1_loc = bp1_loc - 1; % Path location in terms of re_kx in kx plane
c1_start_real = kxx(bp1_loc - 1);
c1_start_imag = -kxx(bp1_loc - 1)*1e3;
diff = abs(kxy - c1_start_imag);
c1_sty_loc = find(diff == min(diff)); % Starting point on the Im_kx axis

%
% Define real and imaginary end points for c1
%
c1_end_real = kxx(bp1_loc - 1);
diff = abs(kxy - 0);
c1_eny_loc = find(diff == min(diff)); % Ending point on the Im_kx axis
c1_end_imag = kxy(c1_eny_loc - 1)/len;
%
% Make C1 Contour
c1_real = c1_start_real*ones(len,1);
c1_imag = linspace(c1_start_imag, c1_end_imag,len);
%
c1 = horzcat(c1_real, c1_imag'); % concatenate real and imaginary parts
kx_c1 = c1_real(:,1) + 1i*c1_imag(1,:)';

%
% Define real and imaginary start points for c2
%
c2_start_real = kxx(bp1_loc - 1);
c2_start_imag = kxy(c1_eny_loc + 1)/len;
diff = abs(kxy - c2_start_imag);
c2_y_loc = find(diff == min(diff)); % Starting point on the Im_kx axis
%
% Define real and imaginary end points for c2
%
c2_end_real = kxx(bp1_loc + 1);
c2_end_imag = kxy(c1_eny_loc + 1)/len;
%
%
% Make C3 Contour
%
c2_real = linspace(c2_start_real, c2_end_real,3);
c2_imag = c2_end_imag*ones(3,1);

c2 = horzcat(c2_real', c2_imag);
kx_c2 = c2_real(1,:)' + 1i*c2_imag(1,:);
% Define real and imaginary start points for c3
%
c3_start_real = kxx(bp1_loc + 1);
c3_start_imag = kxy(c1_eny_loc - 1)/len;
%
% Define real and imaginary end points for c3
%
c3_end_real = kxx(bp1_loc + 1);
c3_end_imag = -kxx(bp1_loc - 1)*1e3;
%
% Make C3 Contour
c3_real = c3_start_real*ones(len,1);
c3_imag = linspace(c3_start_imag, c3_end_imag,len);

c3 = horzcat(c3_real, c3_imag');
kx_c3 = c3_real(:,1) + 1i*c3_imag(1,:)';



%% Define Green's function


kz_1 = @(kx) sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
D = @(kz_1, kz_2) kz_2/eps_silver + kz_1/eps_0;
G = @(kz_1, kz_2) 1./D;

% On Contour C1
kz1_c1 = kz_1(kx_c1);
kz2_c1 = kz_2(kx_c1);

% On Contour C2
kz1_c2 = kz_1(kx_c2);
kz2_c2 = kz_2(kx_c2);

% On Contour C3
kz1_c3 = kz_1(kx_c3);
kz2_c3 = kz_2(kx_c3);

%%
% Branch Cut Curve
hyp_silver = imag(k_silver^2)./(2*kxx); % Hyperbolic cruve for silver

% Intersection on C1
flip_c1_imag = flip(c1_imag);
diff1 = abs(hyp_silver - flip_c1_imag);
c1_hyp_silver_int_loc = find(diff1 == min(diff1));
Silver_branch_cut_loc_c1 = flip_c1_imag(c1_hyp_silver_int_loc);

% Intersection on C3
diff2 = abs(hyp_silver - c3_imag);
c3_hyp_silver_int_loc = find(diff2 == min(diff2));
Silver_branch_cut_loc_c3 = c3_imag(c3_hyp_silver_int_loc);

%% Integrate
for i = 1 : length (x)
    % Integrate on left edge
    % This is a bottom (improper) Riemann sheet
    % Im(kz) > 0
    % Due to the path being on the left of the cut
    % Re(kz) > 0 [1]
    for j = 1 : length(kx_c1)
        % Enforce real parts to be wavevectors to be negative
        %
        if imag(kx_c1) < Silver_branch_cut_loc_c1
            % This is bottom (improper) Riemann sheet of both kz1 and kz2
            % Im(kz) > 0
            % Due to the path being on the left of the cut
            % Re(kz) > 0 [1]
            if real(kz1_c1(j)) < 0
                kz1_c1(j) = -real(kz1_c1(j)) + 1i*imag(kz1_c1(j));
            end
            if real(kz2_c1(j)) < 0
                kz2_c1(j) = -real(kz2_c1(j)) + 1i*imag(kz2_c1(j));
            end
            %
            % Satisfy Imaginary parts
            if imag(kz1_c1(j)) < 0
                kz1_c1(j) = conj(kz1_c1(j));
            end
            if imag(kz2_c1(j)) < 0
                kz2_c1(j) = conj(kz2_c1(j));
            end
        else
            % This is bottom sheet of kz1, but top sheet ofkz2
            % Im(kz1) > 0, Re(kz1) > 0
            % Im(kz2) < 0, Re(kz2) < 0
            if real(kz1_c1(j)) < 0
                kz1_c1(j) = -real(kz1_c1(j)) + 1i*imag(kz1_c1(j));
            end
            if real(kz2_c1(j)) > 0
                kz2_c1(j) = -real(kz2_c1(j)) + 1i*imag(kz2_c1(j));
            end
            %
            % Satisfy Imaginary parts
            if imag(kz1_c1(j)) < 0
                kz1_c1(j) = conj(kz1_c1(j));
            end
            if imag(kz2_c1(j)) > 0
                kz2_c1(j) = conj(kz2_c1(j));
            end
        end
        
        
    end
    D_c1 = D(kz1_c1, kz2_c1);
    G_c1 = 1./D_c1;
    integrand_1 = G_c1.*exp(-1i*kx_c1(j)*x(i));
    H_c1(i) = sum(integrand_1(1: length(integrand_1 - 3)));
    
    % Integrate on top edge
    % This is a top (improper) Riemann sheet for both
    % Im(kz) > 0
    % Due to the path being on the left of the cut
    % Re(kz) < 0 [1]
    for j = 1 : length(kx_c2)
        % Enforce real parts to be wavevectors to be negative
        %
        if real(kz1_c2(j)) < 0
            kz1_c2(j) = -real(kz1_c2(j)) + 1i*imag(kz1_c2(j));
        end
        if real(kz2_c2(j)) > 0
            kz2_c2(j) = -real(kz2_c2(j)) + 1i*imag(kz2_c2(j));
        end
        %
        % Satisfy Imaginary parts
        if imag(kz1_c2(j)) > 0
            kz1_c2(j) = conj(kz1_c2(j));
        end
        if imag(kz2_c2(j)) < 0
            kz2_c2(j) = conj(kz2_c2(j));
        end
    end
    D_c2 = D(kz1_c2, kz2_c2);
    G_c2 = 1./D_c2;
    integrand_2 = G_c2.*exp(-1i*kx_c2*x(i));
    H_c2(i) = sum(integrand_2);
    
    % Integrate on right edge
    % This is a top (proper) Riemann sheet
    % Im(kz) < 0
    % Due to the path being on the left of the cut
    % Re(kz) < 0 [1]
    for j = 1 : length(kx_c3)
        
        
        if imag(kx_c1) > Silver_branch_cut_loc_c3
            % This is the top sheet of both kz1 and kz2
            % Im(kz1) < 0, Re(kz1) < 0
            % Im(kz2) < 0, Re(kz2) < 0
            if real(kz1_c3(j)) > 0
                kz1_c3(j) = -real(kz1_c3(j)) + 1i*imag(kz1_c3(j));
            end
            if real(kz2_c3(j)) < 0
                kz2_c3(j) = -real(kz2_c3(j)) + 1i*imag(kz2_c3(j));
            end
            %
            % Satisfy Imaginary parts
            if imag(kz1_c1(j)) < 0
                kz1_c1(j) = conj(kz1_c1(j));
            end
            if imag(kz2_c1(j)) > 0
                kz2_c1(j) = conj(kz2_c1(j));
            end
        else
            % This is bottom sheet of of both kz1 and kz2
            % Im(kz1) < 0, Re(kz1) > 0
            % Im(kz2) > 0, Re(kz2) > 0
            if real(kz1_c1(j)) < 0
                kz1_c1(j) = -real(kz1_c1(j)) + 1i*imag(kz1_c1(j));
            end
            if real(kz2_c1(j)) < 0
                kz2_c1(j) = -real(kz2_c1(j)) + 1i*imag(kz2_c1(j));
            end
            %
            % Satisfy Imaginary parts
            if imag(kz1_c1(j)) < 0
                kz1_c1(j) = conj(kz1_c1(j));
            end
            if imag(kz2_c1(j)) < 0
                kz2_c1(j) = conj(kz2_c1(j));
            end
        end        
    end
    D_c3 = D(kz1_c3, kz2_c3);
    G_c3 = 1./D_c3;
    integrand_3 = G_c3.*exp(-1i*kx_c3*x(i));
    H_c3(i) = sum(integrand_3);
    
end
%% Define H by EQ. 1 in [1]
H = -k_air/(2*pi*eta_0)*(H_c1 +  H_c3); % C2 Contour must be zero anyway, therefore taken out
%% Plot Figure
figure('Name','Decay of the Creeping Wave');
loglog(x/lambda, abs(H)/abs(max(H)),'LineWidth',1.4,'Color','black')
set(gcf,'Color','white');

ylabel('$\vert Creeping Wave\vert$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\frac{x}{\lambda}$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

ylim([10e-10 10e1])

title('Decay Plot of Creeping Wave part');

cleanfigure();
matlab2tikz('filename',sprintf('nevels_michalski_decay_plot.tex'))
%%
save nevels_michalski_creeping_wave.mat % Save data to plot the branch cuts
save('test.mat', 'H', 'x')              % Save H and x for plotting along SPP
%% Transpose All variables for all variables
% kz_1 = kz_1.';
% kz_2 = kz_2.';
% D = D.';
% G = G.';
% H = H.';
% kx = kx.';
