function cost = antenna_sim(parameters, Sim_Path)

display_geom = 0;
run_sim = 1;

physical_constants;

Sim_CSX = 'uwb.xml';

fov = 90 * pi/180;

f_0 = 6.4896e9;
f_c = 0.4992e9 / 0.3925;

helix_radius = round(parameters(1)*1e6)/1e6;
helix_length = round(parameters(2)*1e6)/1e6;
helix_pitch2 = round(parameters(3)*1e6)/1e6;
helix_pitch1 = round(parameters(4)*1e6)/1e6;

helix_turns = 0.25*(4.0*helix_length - helix_pitch1 + helix_pitch2)/helix_pitch2;
taper_ratio = 1;

lambda0 = C0/f_0;

feed_height = 1e-3;
wire_radius = 5e-4;
ground_radius = 3e-2;
ground_height = 1e-3;
feed_R = 50;

ang = linspace(0,pi/2,50);
helix_points1(1,:) = helix_radius*cos(ang);
helix_points1(2,:) = helix_radius*sin(ang);
helix_points1(3,:) = feed_height + ang/(2*pi)*helix_pitch1;

ang = linspace(pi/2,2*pi*helix_turns,200*(helix_turns-.25));

helix_points2(1,:) = helix_radius*cos(ang);
helix_points2(2,:) = helix_radius*sin(ang);
helix_points2(3,:) = feed_height + .25*helix_pitch1 + (ang/(2*pi)-.25)*helix_pitch2;

mesh_res = C0 / (f_0 + f_c) / 10;

CSX = InitCSX();
CSX = AddMetal(CSX, 'Ground');
CSX = AddCylinder(CSX,'Ground',1,[0 0 -ground_height],[0 0 0],ground_radius);

CSX = AddMetal(CSX, 'Helix');
CSX = AddCurve(CSX, 'Helix', 1, helix_points1);
CSX = AddWire(CSX, 'Helix', 1, helix_points1, wire_radius);
CSX = AddCurve(CSX, 'Helix', 1, helix_points2);
CSX = AddWire(CSX, 'Helix', 1, helix_points2, wire_radius);

CSX = AddMaterial(CSX, 'Structure');
CSX = SetMaterialProperty(CSX, 'Structure', 'Epsilon', 3);
inside_radius = helix_radius-wire_radius-4e-4;
outside_radius = helix_radius+wire_radius;
CSX = AddCylindricalShell(CSX,'Structure',0,[0 0 0],[0 0 feed_height+helix_length], (inside_radius+outside_radius)/2, outside_radius-inside_radius);

[CSX port] = AddLumpedPort(CSX, 999, 1, feed_R, [helix_radius 0 0], [helix_radius 0 feed_height], [0 0 1], true);

mesh.x = [-(2*lambda0+ground_radius) SmoothMeshLines2([-helix_radius, 0, helix_radius], wire_radius) (2*lambda0+ground_radius)];
mesh.y = mesh.x;
mesh.z = [-(2*lambda0+ground_height) -ground_height -ground_height/2 0 feed_height/2 SmoothMeshLines2([feed_height, feed_height+helix_length], wire_radius) (2*lambda0+feed_height+helix_length)];

mesh = SmoothMesh(mesh, mesh_res);
CSX = DefineRectGrid(CSX, 1, mesh);

start = [mesh.x(11)     mesh.y(11)     mesh.z(11)];
stop  = [mesh.x(end-10) mesh.y(end-10) mesh.z(end-10)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

FDTD = InitFDTD( 'NrTs', 50000, 'EndCriteria', 1e-5, 'OverSampling', 2);
FDTD = SetGaussExcite(FDTD, f_0, f_c );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond(FDTD, BC);

if (run_sim==1 || display_geom==1)
    confirm_recursive_rmdir(0);
    [status, message, messageid] = rmdir( Sim_Path, 's' );
    [status, message, messageid] = mkdir( Sim_Path );

    WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

    if (display_geom==1)
        CSXGeomPlot([Sim_Path '/' Sim_CSX]);
        return;
    end

    RunOpenEMS( Sim_Path, Sim_CSX);
end

freq = linspace(f_0-0.4992e9, f_0+0.4992e9, 11);
theta = linspace(0, fov/2, 5);
phi = linspace(-pi, pi, 10);

port = calcPort(port, Sim_Path, freq);
nf2ff = CalcNF2FF(nf2ff, Sim_Path, freq, theta, phi,'Mode',1);

s11 = port.uf.ref ./ port.uf.inc;

s11_db = 20*log10(abs(s11))

vswr = (1+abs(s11)) ./ (1-abs(s11))
max_vswr = max(vswr)

rhcp_gain = 4*pi*(abs(nf2ff.E_cprh{5}).^2./abs(nf2ff.E_norm{5}).^2 .* nf2ff.P_rad{5}) ./ port.P_inc(5)

min_rhcp_gain = min(min(rhcp_gain))

vswr_target = 2.5;

cost = max([1, 1+(max_vswr-vswr_target)*10])/min_rhcp_gain
return;
