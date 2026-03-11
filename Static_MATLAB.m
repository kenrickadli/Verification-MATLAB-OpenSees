% =========================================================================
% 3D PORTAL FRAME - Static Analysis (Direct Stiffness Method)
% =========================================================================

clear; clc;

% Units: Force (N), Length (m), Time (sec)
g = 9.80665;  % m/s²

fprintf('\n=== STATIC ANALYSIS STARTING ===\n');

%% ========================================================================
% 1. STRUCTURE GEOMETRY
% =========================================================================
Lx = 7.2;   % Bay width (m)
Ly = 7.2;   % Bay depth (m)
Lz = 3.6;   % Story height (m)
nx = 1;     % Number of bays in X
ny = 1;     % Number of bays in Y
nz = 1;     % Number of stories

ndt = (nx + 1) * (ny + 1) * (nz + 1);
nds = (nx + 1) * (ny + 1);

fprintf('\nStructure: %dx%d bays, %d stories\n', nx, ny, nz);
fprintf('Total nodes: %d\n', ndt);

%% ========================================================================
% 2. LOADING
% =========================================================================
SIDL = 2.0;     % kN/m²
Live = 2.4;     % kN/m²
h_slab = 180;   % mm
rho_c = 2400;   % kg/m³

floor_area = (Lx * nx) * (Ly * ny);
m_Dead = rho_c * (h_slab * 1e-3) * floor_area;
m_SIDL = (SIDL * 1e3) * floor_area / g;
m_Live = (Live * 1e3) * floor_area / g;
mass_source = (m_Dead + m_SIDL + 0.5 * m_Live);
mz = mass_source * ((Lx * nx)^2 + (Ly * ny)^2) / 12;

fprintf('Total mass per floor: %.0f kg\n', mass_source);

%% ========================================================================
% 3. MATERIAL PROPERTIES
% =========================================================================
Fy = 345;       % MPa
Ry = 1.1;
Es = 200000;    % MPa
v = 0.3;

Es = Es * 1e6;  % Pa
G = Es / (2 * (1 + v));
J = 1e-10;

fprintf('\nMaterial: E = %.2e Pa, G = %.2e Pa\n', Es, G);

%% ========================================================================
% 4. BEAM SECTION PROPERTIES
% =========================================================================
d_beam = 312.4;   % mm
bf_beam = 304.8;  % mm
tf_beam = 17;     % mm
tw_beam = 10.9;   % mm

A_beam = (bf_beam * tf_beam * 2 + (d_beam - 2 * tf_beam) * tw_beam) * 1e-6;  % m²
Ix_beam = ((bf_beam * d_beam^3 / 12) - ((bf_beam - tw_beam) * (d_beam - 2 * tf_beam)^3 / 12)) * 1e-12;  % m⁴
Iy_beam = ((d_beam - 2 * tf_beam) * tw_beam^3 / 12 + 2 * (tf_beam * bf_beam^3 / 12)) * 1e-12;  % m⁴

fprintf('\n=== BEAM PROPERTIES ===\n');
fprintf('Area: %.5e m²\n', A_beam);
fprintf('Ix: %.5e m⁴, Iy: %.5e m⁴\n', Ix_beam, Iy_beam);

%% ========================================================================
% 5. DEFINE NODES
% =========================================================================
fprintf('\n--- Building Model: Nodes ---\n');

num_joints = ndt;
joint_coords = zeros(num_joints, 3);

n = 1;
for k = 0:nz
    for j = 0:ny
        for i = 0:nx
            joint_coords(n, 1) = i * Lx;
            joint_coords(n, 2) = j * Ly;
            joint_coords(n, 3) = k * Lz;
            n = n + 1;
        end
    end
end

fprintf('Created %d nodes\n', num_joints);

%% ========================================================================
% 6. DEFINE SUPPORTS
% =========================================================================
fprintf('--- Building Model: Supports ---\n');

num_supports = nds;
support_data = zeros(num_supports, 7);

for i = 1:nds
    support_data(i, 1:7) = [i 1 1 1 1 1 1];  % Fixed
end

fprintf('Created %d fixed supports at base\n', num_supports);

%% ========================================================================
% 7. DEFINE MEMBERS
% =========================================================================
fprintf('--- Building Model: Members ---\n');

member_data = [];

% Columns
fprintf('Creating columns...\n');
nn = 1;
for k = 0:(nz-1)
    for j = 0:ny
        for i = 0:nx
            node1 = nn;
            node2 = nn + nds;
            member_data = [member_data; node1, node2, 1, 1, 0];
            nn = nn + 1;
        end
    end
end
num_columns = size(member_data, 1);

% Beams X-direction
fprintf('Creating beams (X-direction)...\n');
nn = nds + 1;
for k = 1:nz
    for j = 0:ny
        for i = 0:nx
            if i < nx
                node1 = nn;
                node2 = nn + 1;
                member_data = [member_data; node1, node2, 1, 1, 90];
            end
            nn = nn + 1;
        end
    end
end

% Beams Y-direction
fprintf('Creating beams (Y-direction)...\n');
nn = nds + 1;
for k = 1:nz
    for j = 0:ny
        for i = 0:nx
            if j < ny
                node1 = nn;
                node2 = nn + (nx + 1);
                member_data = [member_data; node1, node2, 1, 1, 90];
            end
            nn = nn + 1;
        end
    end
end

num_members = size(member_data, 1);
fprintf('Total members created: %d (%d columns, %d beams)\n', ...
    num_members, num_columns, num_members - num_columns);

material_props = [Es, G] / 1e9;  % GPa
cross_section_props = [A_beam*1e6, Ix_beam*1e12, Iy_beam*1e12, J*1e12];  % mm units

%% ========================================================================
% 8. DEFINE LOADS
% =========================================================================
fprintf('--- Building Model: Loads ---\n');

nodal_load = 1/4 * mass_source * g / (nx * ny);  % N

num_joint_loads = (ndt - nds);
joint_load_nums = (nds+1:ndt)';
joint_loads = zeros(num_joint_loads, 6);
joint_loads(:, 3) = -nodal_load;  % Z-direction

fprintf('Applied %.2f N vertical load to %d nodes\n', nodal_load, num_joint_loads);

%% ========================================================================
% 9. DOF NUMBERING
% =========================================================================
fprintf('\n--- Setting up Analysis ---\n');

DOF_PER_JOINT = 6;
total_nodes = num_joints;

num_reactions = num_supports * 6;
total_structural_dof = DOF_PER_JOINT * total_nodes;
num_dof = total_structural_dof - num_reactions;

fprintf('Total DOF: %d (Free: %d, Restrained: %d)\n', ...
    total_structural_dof, num_dof, num_reactions);

structure_coords = zeros(DOF_PER_JOINT * total_nodes, 1);
free_dof_counter = 0;
restrained_dof_counter = num_dof;

% Number original joints
for i = 1:num_joints
    is_supported = any(support_data(:, 1) == i);
    
    if is_supported
        for dof = 1:DOF_PER_JOINT
            coord_index = (i-1) * DOF_PER_JOINT + dof;
            restrained_dof_counter = restrained_dof_counter + 1;
            structure_coords(coord_index) = restrained_dof_counter;
        end
    else
        for dof = 1:DOF_PER_JOINT
            coord_index = (i-1) * DOF_PER_JOINT + dof;
            free_dof_counter = free_dof_counter + 1;
            structure_coords(coord_index) = free_dof_counter;
        end
    end
end

%% ========================================================================
% 10. ASSEMBLE GLOBAL STIFFNESS MATRIX
% =========================================================================
fprintf('--- Assembling Global Stiffness Matrix ---\n');

K_global = zeros(num_dof, num_dof);
P_global = zeros(num_dof, 1);

% Convert material props back to Pa
material_props = material_props * 1e9;
cross_section_props(1) = cross_section_props(1) * 1e-6;  % m²
cross_section_props(2:4) = cross_section_props(2:4) * 1e-12;  % m⁴

for member_id = 1:num_members
    joint_begin = member_data(member_id, 1);
    joint_end = member_data(member_id, 2);
    roll_angle = member_data(member_id, 5);
    
    E = material_props(1);
    G = material_props(2);
    Area = cross_section_props(1);
    Ix = cross_section_props(2);
    Iy = cross_section_props(3);
    J = cross_section_props(4);
    
    x_begin = joint_coords(joint_begin, 1);
    y_begin = joint_coords(joint_begin, 2);
    z_begin = joint_coords(joint_begin, 3);
    x_end = joint_coords(joint_end, 1);
    y_end = joint_coords(joint_end, 2);
    z_end = joint_coords(joint_end, 3);
    
    member_length = sqrt((x_end - x_begin)^2 + (y_end - y_begin)^2 + (z_end - z_begin)^2);
    
    cos_x = (x_end - x_begin) / member_length;
    cos_y = (y_end - y_begin) / member_length;
    cos_z = (z_end - z_begin) / member_length;
    
    if abs(cos_x) < 1e-6 && abs(cos_z) < 1e-6
        cos_yx = -cos_y * cosd(roll_angle);
        cos_yy = 0;
        cos_yz = sind(roll_angle);
        cos_zx = cos_y * sind(roll_angle);
        cos_zy = 0;
        cos_zz = cosd(roll_angle);
    else
        denom = sqrt(cos_x^2 + cos_z^2);
        cos_yx = (-cos_x * cos_y * cosd(roll_angle) - cos_z * sind(roll_angle)) / denom;
        cos_yy = denom * cosd(roll_angle);
        cos_yz = (-cos_y * cos_z * cosd(roll_angle) + cos_x * sind(roll_angle)) / denom;
        cos_zx = (cos_x * cos_y * sind(roll_angle) - cos_z * cosd(roll_angle)) / denom;
        cos_zy = -denom * sind(roll_angle);
        cos_zz = (cos_y * cos_z * sind(roll_angle) + cos_x * cosd(roll_angle)) / denom;
    end
    
    K_local = form_local_stiffness(E, G, Area, Ix, Iy, J, member_length);
    T_matrix = form_transformation_matrix(cos_x, cos_y, cos_z, ...
        cos_yx, cos_yy, cos_yz, cos_zx, cos_zy, cos_zz);
    K_member_global = T_matrix' * K_local * T_matrix;
    
    K_global = assemble_member(K_global, K_member_global, joint_begin, joint_end, ...
        structure_coords, num_dof);
end

fprintf('Global stiffness matrix assembled (%dx%d)\n', num_dof, num_dof);

%% ========================================================================
% 11. APPLY LOADS
% =========================================================================
fprintf('--- Applying Loads ---\n');

for i = 1:num_joint_loads
    joint_num = joint_load_nums(i);
    
    for dof = 1:6
        coord_index = (joint_num - 1) * 6 + dof;
        global_dof = structure_coords(coord_index);
        
        if global_dof <= num_dof
            P_global(global_dof) = P_global(global_dof) + joint_loads(i, dof);
        end
    end
end

fprintf('Loads applied to %d nodes\n', num_joint_loads);

%% ========================================================================
% 12. SOLVE SYSTEM
% =========================================================================
fprintf('\n--- Solving System ---\n');

tic;
D_global = K_global \ P_global;
solve_time = toc;

fprintf('Solution completed in %.4f seconds\n', solve_time);

%% ========================================================================
% 13. EXTRACT RESULTS
% =========================================================================
fprintf('\n--- Extracting Results ---\n');

joint_displacements = zeros(num_joints, 6);
for i = 1:num_joints
    for dof = 1:6
        coord_index = (i - 1) * 6 + dof;
        global_dof = structure_coords(coord_index);
        
        if global_dof <= num_dof
            joint_displacements(i, dof) = D_global(global_dof);
        end
    end
end

fprintf('\n=== FLOOR NODE DISPLACEMENTS (mm) ===\n');
displacements = [];
for i = (nds+1):ndt
    uz = joint_displacements(i, 3) * 1000;
    displacements = [displacements, uz];
    fprintf('Node %d: Uz = %.6f mm\n', i, uz);
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');

%% ========================================================================
% FUNCTIONS
% =========================================================================
function K = form_local_stiffness(E, G, A, Ix, Iy, J, L)
    K = zeros(12, 12);
    
    % Axial
    K([1,7], [1,7]) = (E*A/L) * [1, -1; -1, 1];
    
    % Shear Y & Bending Z
    K([2,8], [2,8]) = (12*E*Ix/L^3) * [1, -1; -1, 1];
    K([6,12], [6,12]) = (4*E*Ix/L) * [1, 0.5; 0.5, 1];
    K([2,8], [6,12]) = (6*E*Ix/L^2) * [1, 1; -1, -1];
    K([6,12], [2,8]) = K([2,8], [6,12])';
    
    % Shear Z & Bending Y
    K([3,9], [3,9]) = (12*E*Iy/L^3) * [1, -1; -1, 1];
    K([5,11], [5,11]) = (4*E*Iy/L) * [1, 0.5; 0.5, 1];
    K([3,9], [5,11]) = (-6*E*Iy/L^2) * [1, 1; -1, -1];
    K([5,11], [3,9]) = K([3,9], [5,11])';
    
    % Torsion
    K([4,10], [4,10]) = (G*J/L) * [1, -1; -1, 1];
end

function T = form_transformation_matrix(cx, cy, cz, cyx, cyy, cyz, czx, czy, czz)
    R = [cx, cy, cz;
         cyx, cyy, cyz;
         czx, czy, czz];
    
    T = zeros(12, 12);
    T(1:3, 1:3) = R;
    T(4:6, 4:6) = R;
    T(7:9, 7:9) = R;
    T(10:12, 10:12) = R;
end

function K_global = assemble_member(K_global, K_member, node_i, node_j, coords, num_dof)
    for i = 1:12
        if i <= 6
            global_i = (node_i - 1) * 6 + i;
        else
            global_i = (node_j - 1) * 6 + (i - 6);
        end
        row = coords(global_i);
        
        if row <= num_dof
            for j = 1:12
                if j <= 6
                    global_j = (node_i - 1) * 6 + j;
                else
                    global_j = (node_j - 1) * 6 + (j - 6);
                end
                col = coords(global_j);
                
                if col <= num_dof
                    K_global(row, col) = K_global(row, col) + K_member(i, j);
                end
            end
        end
    end
end
