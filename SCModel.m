%% SCModel.m
%
% Authors:
% Pasquale Marzaioli
% Paolo Iaccarino
% Federica Pirozzi
% Silvia Preziosi
%
% Description:
% This function defines the spacecraft (SC) structural model, including geometry, mass
% properties, and surface characteristics. It computes the inertia matrix, surface normals,
% force radii, optical properties (diffuse and specular coefficients), and center of mass
% shift. The model represents a large satellite with a cubic body and deployable solar
% panels, used in the attitude dynamics simulation.
%
% Outputs:
% - EXTRACTED_DATA: Structure with extracted model data (inertia, surfaces, etc.)
% - MODEL_SC: Structure array with per-surface data

function [EXTRACTED_DATA, MODEL_SC] = SCModel
    
    % -------------------------------
    % Initial parameters
    % -------------------------------
    clc;
    num_body  = 6;                     % Number of main body surfaces
    num_panel = 4;                     % Number of solar panels
    n_sur = num_body + num_panel;      % Total number of surfaces
    
    mass_body = 500;                   % Mass of the main body [kg]
    mass_panel = 3;                    % Mass of a solar panel [kg]
    
    % Dimensions of the main body (cuboid)
    length_body = 1;   % Length [m]
    width_body  = 1;   % Width [m]
    height_body = 1;   % Height [m]
    surface_value_body = width_body * height_body; % Area of one surface of the main body
    
    % Dimensions of the panels
    length_panel = 2;  % Length [m]
    width_panel  = 2;  % Width [m]
    surface_value_panel = length_panel * width_panel; % Area of a solar panel
    
    % -------------------------------
    % Surface area vector
    % -------------------------------
    surface = [ones(num_body, 1) * surface_value_body; ... % Areas of the main body
               ones(num_panel, 1) * surface_value_panel];  % Areas of the panels
    
    % -------------------------------
    % Moment of inertia calculations
    % -------------------------------
    % Moments of inertia for the main body
    I_body = diag([
        (mass_body / 12) * (width_body^2 + height_body^2), % Ixx
        (mass_body / 12) * (length_body^2 + height_body^2), % Iyy
        (mass_body / 12) * (length_body^2 + width_body^2)   % Izz
    ]);
    
    % Moments of inertia for a solar panel
    I_panel = diag([
        (mass_panel / 12) * (width_panel^2 + 0),  % Ixx (width)
        (mass_panel / 12) * (length_panel^2 + 0), % Iyy (length)
        (mass_panel / 12) * 0                     % Izz negligible
    ]);
    
    % Total moment of inertia (body + panels)
    I = I_body + num_panel * I_panel;
    InvI = inv(I);  % Inverse inertia matrix
    
    % -------------------------------
    % Preparation of vectors and matrices
    % -------------------------------
    N_B = [  1  0  0  1
             0  1  0  1
            -1  0  0  1
             0 -1  0  1
             0  0  1  1
             0  0 -1  1
             1  0  0  0
            -1  0  0  0
             1  0  0  0
            -1  0  0  0];  % Surface normals
    
    cm_shift = [0.0, -0.1, 0.0]; % Center of mass shift vector [m]
    % Initialization of vectors
    rho_S = zeros(n_sur, 1);          % Reflection coefficients
    rho_D = zeros(n_sur, 1);          % Dissipation coefficients
    r_F = zeros(n_sur, 3);            % Force vector radii
    type = cell(n_sur, 1);            % Surface type ('body' or 'panel')
    
    % -------------------------------
    % Calculation of force radii and types
    % -------------------------------
    for i = 1:n_sur
        if N_B(i, 4) == 1 % Main body
            type{i}  = 'body';
            rho_S(i) = 0.5;
            rho_D(i) = 0.3;
    
            if N_B(i, 1) == 1
                r_F(i, :) = [length_body / 2, 0, 0] + cm_shift;
            elseif N_B(i, 1) == -1
                r_F(i, :) = [-length_body / 2, 0, 0] + cm_shift;
            elseif N_B(i, 2) == 1
                r_F(i, :) = [0, width_body / 2, 0] + cm_shift;
            elseif N_B(i, 2) == -1
                r_F(i, :) = [0, -width_body / 2, 0] + cm_shift;
            elseif N_B(i, 3) == 1
                r_F(i, :) = [0, 0, height_body / 2] + cm_shift;
            elseif N_B(i, 3) == -1
                r_F(i, :) = [0, 0, -height_body / 2] + cm_shift;
            end
        else % Solar panels
            type{i}  = 'panel';
            rho_S(i) = 0.1;
            rho_D(i) = 0.05;
    
            % Calculation of relative positions of panels
            % panel_offset = width_body / 2 + width_panel / 2;
            panel_offset = 0;
            if i == 7 % First panel
                r_F(i, :) = [0, length_panel / 2, panel_offset] + cm_shift;
            elseif i == 8 % Second panel
                r_F(i, :) = [0, length_panel / 2, panel_offset] + cm_shift;
            elseif i == 9 % Third panel
                r_F(i, :) = [0, -length_panel / 2, panel_offset] + cm_shift;
            elseif i == 10 % Fourth panel
                r_F(i, :) = [0, -length_panel / 2, panel_offset] + cm_shift;
            end
        end
    end
    
    % -------------------------------
    % Updating the MODEL_SC structure
    % -------------------------------
    MODEL_SC = struct();
    for i = 1:n_sur
        MODEL_SC(i).N_B     = N_B(i, 1:3);
        MODEL_SC(i).rho_S   = rho_S(i);
        MODEL_SC(i).rho_D   = rho_D(i);
        MODEL_SC(i).surface = surface(i);
        MODEL_SC(i).r_F     = r_F(i, :);
        MODEL_SC(i).type    = type{i};
    end
    
    % -------------------------------
    % Consolidating data
    % -------------------------------
    EXTRACTED_DATA = struct();
    EXTRACTED_DATA.N_B       = N_B(:, 1:3);
    EXTRACTED_DATA.type      = type;
    EXTRACTED_DATA.rho_S     = rho_S;
    EXTRACTED_DATA.rho_D     = rho_D;
    EXTRACTED_DATA.surface   = surface;
    EXTRACTED_DATA.r_F       = r_F;
    EXTRACTED_DATA.I_total   = I;
    EXTRACTED_DATA.cm_shift  = cm_shift;
    
    % -------------------------------
    % Display
    % -------------------------------
    disp('Total Moments of Inertia [kg·m²]:');
    disp(I);
    
    disp('Center of Mass Shift [m]:');
    disp(cm_shift);
    
    disp('Force Vectors with Shift [m]:');
    disp(r_F);
    
    disp('Updated MODEL_SC structure:');
    disp(MODEL_SC);
end