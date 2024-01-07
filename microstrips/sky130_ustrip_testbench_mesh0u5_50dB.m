% sky130 Microstrip testbench on openEMS
% Author: Leonardo Gomes
% Based on Volker Muelhaus' simulation script found on 
% https://github.com/IHP-GmbH/IHP-Open-PDK/tree/dev/ihp-sg13g2/libs.tech/openems/testcase/SG13_Octagon_L2n0/OpenEMS_Matlab

close all
clear
clc
confirm_recursive_rmdir(0);   % delete old data without asking

basename = mfilename ; % get name of current model from *.m filename
physical_constants;    % load table of physical constants

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unit = 1e-6; % specify everything in um

Boundaries   = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'};  % xmin xmax ymin ymax zmin zmax

refined_cellsize = 0.25;  % mesh.z resolution inside thick metals to account for the skin effect

f_start = 0e9;
f_stop  = 110e9;

energy_limit = -50;    % end criteria for residual energy

max_cellsize = c0/(f_stop*sqrt(11.9))/unit /20;  % max cellsize 1/20 wavelength in silicon
converg_port_length = 5; % denser mesh around ports to improve accuracy

ustrip_width = 6;

used_metals = [1 0 0 0 1];
used_vias = [0 0 0 0 ];


%% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
lim = exp(energy_limit/10 * log(10)); % cconvert energy limit from dB to number
FDTD = InitFDTD('endCriteria',lim);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
FDTD = SetBoundaryCond( FDTD, Boundaries );

CSX = InitCSX();

%%%%%%%%%%%%%%%%%%%%%%%%%  sky130 stackup  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_M = [1.376 2.006 2.786 4.021 5.371];
th_M = [0.36 0.36 0.845 0.845 1.26];

%% silicon substrate
CSX = AddMaterial( CSX, 'Sub' );
CSX = SetMaterialProperty( CSX, 'Sub', 'Epsilon', 11.9, 'Kappa', 0.3 );
Sub.thick = 200;
Sub.zmin = -283.75;
Sub.zmax = Sub.zmin + Sub.thick;

mesh.z = [linspace(Sub.zmin,Sub.zmax,20) ];


%% EPI -- even though sky130 doesn't have an EPI layer, 
%%        I decided to keep it to enforce a thinner mesh at the BEOL/FEOL interface
CSX = AddMaterial( CSX, 'EPI' );
CSX = SetMaterialProperty( CSX, 'EPI', 'Epsilon', 11.9, 'Kappa', 0.3 );
EPI.thick = 3.75;
EPI.zmin = Sub.zmax;
EPI.zmax = EPI.zmin + EPI.thick;

mesh.z = [mesh.z linspace(EPI.zmin,EPI.zmax,3)];


%% FOX
CSX = AddMaterial( CSX, 'FOX' );
CSX = SetMaterialProperty( CSX, 'FOX', 'Epsilon', 3.9 );
FOX.thick = 0.936;
FOX.zmin = EPI.zmax;
FOX.zmax = FOX.zmin + FOX.thick;

mesh.z = [mesh.z FOX.zmax];

%% Sub_LI
CSX = AddMaterial( CSX, 'Sub_LI' );
CSX = SetMaterialProperty( CSX, 'Sub_LI', 'Epsilon', 4.383 );
Sub_LI.thick = 0.44;
Sub_LI.zmin = FOX.zmax;
Sub_LI.zmax = Sub_LI.zmin + Sub_LI.thick;

mesh.z = [mesh.z Sub_LI.zmax];

%% Sub_Metal1
CSX = AddMaterial( CSX, 'Sub_Metal1' );
CSX = SetMaterialProperty( CSX, 'Sub_Metal1', 'Epsilon', 4.5 );
Sub_Metal1.thick = 0.63;
Sub_Metal1.zmin = Sub_LI.zmax;
Sub_Metal1.zmax = Sub_Metal1.zmin + Sub_Metal1.thick;

mesh.z = [mesh.z Sub_Metal1.zmax];

%% Sub_Metal2
CSX = AddMaterial( CSX, 'Sub_Metal2' );
CSX = SetMaterialProperty( CSX, 'Sub_Metal2', 'Epsilon', 4.2 );
Sub_Metal2.thick = 0.78;
Sub_Metal2.zmin = Sub_Metal1.zmax;
Sub_Metal2.zmax = Sub_Metal2.zmin + Sub_Metal2.thick;

mesh.z = [mesh.z Sub_Metal2.zmax];

%% Sub_Metal3
CSX = AddMaterial( CSX, 'Sub_Metal3' );
CSX = SetMaterialProperty( CSX, 'Sub_Metal3', 'Epsilon', 4.1 );
Sub_Metal3.thick = 1.235;
Sub_Metal3.zmin = Sub_Metal2.zmax;
Sub_Metal3.zmax = Sub_Metal3.zmin + Sub_Metal3.thick;

mesh.z = [mesh.z Sub_Metal3.zmax];

%% Sub_Metal4
CSX = AddMaterial( CSX, 'Sub_Metal4' );
CSX = SetMaterialProperty( CSX, 'Sub_Metal4', 'Epsilon', 4 );
Sub_Metal4.thick = 1.35;
Sub_Metal4.zmin = Sub_Metal3.zmax;
Sub_Metal4.zmax = Sub_Metal4.zmin + Sub_Metal4.thick;

mesh.z = [mesh.z Sub_Metal4.zmax];

%% Sub_Metal5
CSX = AddMaterial( CSX, 'Sub_Metal5' );
CSX = SetMaterialProperty( CSX, 'Sub_Metal5', 'Epsilon', 3.5 );
Sub_Metal5.thick = 1.979;
Sub_Metal5.zmin = Sub_Metal4.zmax;
Sub_Metal5.zmax = Sub_Metal5.zmin + Sub_Metal5.thick;

mesh.z = [mesh.z Sub_Metal5.zmax];

%% air above is background material, no need to place box, just add mesh line
Air.thick = 100;
Air.zmax = Sub_Metal5.zmax + Air.thick;
mesh.z = [mesh.z Air.zmax];

%% Metal5
Metal5.sigma = 27400000.0;
Metal5.thick = 1.26;
Metal5.zmin  = FOX.zmin + 5.371;
Metal5.zmax  = Metal5.zmin + Metal5.thick;
CSX = AddMaterial( CSX, 'Metal5' );
CSX = SetMaterialProperty( CSX, 'Metal5', 'Kappa', Metal5.sigma );

if(used_metals(5)==1)
	mesh.z = [mesh.z linspace(Metal5.zmin,Metal5.zmax,(ceil(Metal5.thick/refined_cellsize)))];
endif

%% Metal4
Metal4.sigma = 25200000.0;
Metal4.thick = 0.845;
Metal4.zmin  = FOX.zmin + 4.021;
Metal4.zmax  = Metal4.zmin + Metal4.thick;
CSX = AddMaterial( CSX, 'Metal4' );
CSX = SetMaterialProperty( CSX, 'Metal4', 'Kappa', Metal4.sigma );

if(used_metals(4)==1)
	mesh.z = [mesh.z linspace(Metal4.zmin,Metal4.zmax,(ceil(Metal4.thick/refined_cellsize)))];
endif

%% Metal3
Metal3.sigma = 25200000.0;
Metal3.thick = 0.845;
Metal3.zmin  = FOX.zmin + 2.786;
Metal3.zmax  = Metal3.zmin + Metal3.thick;
CSX = AddMaterial( CSX, 'Metal3' );
CSX = SetMaterialProperty( CSX, 'Metal3', 'Kappa', Metal3.sigma );

if(used_metals(3)==1)
	mesh.z = [mesh.z linspace(Metal3.zmin,Metal3.zmax,(ceil(Metal3.thick/refined_cellsize)))];
endif

%% Metal2
Metal2.sigma = 22200000.0;
Metal2.thick = 0.36;
Metal2.zmin  = FOX.zmin + 2.006;
Metal2.zmax  = Metal2.zmin + Metal2.thick;
CSX = AddMaterial( CSX, 'Metal2' );
CSX = SetMaterialProperty( CSX, 'Metal2', 'Kappa', Metal2.sigma );

if(used_metals(2)==1)
	mesh.z = [mesh.z Metal2.zmax];
endif

%% Metal1
Metal1.sigma = 22200000.0;
Metal1.thick = 0.36;
Metal1.zmin  = FOX.zmin + 1.376;
Metal1.zmax  = Metal1.zmin + Metal1.thick;
CSX = AddMaterial( CSX, 'Metal1' );
CSX = SetMaterialProperty( CSX, 'Metal1', 'Kappa', Metal1.sigma );

if(used_metals(1)==1)
	mesh.z = [mesh.z Metal1.zmax];
endif

%% Via4
Via4.sigma = 2610000.0;
Via4.thick = 0.505;
Via4.zmin  = FOX.zmin + 4.866;
Via4.zmax  = Via4.zmin + Via4.thick;
CSX = AddMaterial( CSX, 'Via4' );
CSX = SetMaterialProperty( CSX, 'Via4', 'Kappa', Via4.sigma );

%% Via3
Via3.sigma = 3760000.0;
Via3.thick = 0.39;
Via3.zmin  = FOX.zmin + 3.363;
Via3.zmax  = Via3.zmin + Via3.thick;
CSX = AddMaterial( CSX, 'Via3' );
CSX = SetMaterialProperty( CSX, 'Via3', 'Kappa', Via3.sigma );

%% Via2
Via2.sigma = 3490000.0;
Via2.thick = 0.42;
Via2.zmin  = FOX.zmin + 2.366;
Via2.zmax  = Via2.zmin + Via2.thick;
CSX = AddMaterial( CSX, 'Via2' );
CSX = SetMaterialProperty( CSX, 'Via2', 'Kappa', Via2.sigma );

%% Via1
Via2.sigma = 3620000.0;
Via2.thick = 0.27;
Via2.zmin  = FOX.zmin + 1.736;
Via2.zmax  = Via2.zmin + Via2.thick;
CSX = AddMaterial( CSX, 'Via1' );
CSX = SetMaterialProperty( CSX, 'Via2', 'Kappa', Via2.sigma );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  conductors  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cell ("test", 2 polygons, 0 paths, 0 labels, 0 references)

clear p;
p(1,1) = -15.000;
p(2,1) = 0.000;
p(1,2) = -15.000;
p(2,2) = 500.000;
p(1,3) = 15.000;
p(2,3) = 500.000;
p(1,4) = 15.000;
p(2,4) = 0.000;
CSX = AddLinPoly( CSX, 'Metal1', 200, 'z', Metal1.zmin, p, Metal1.thick); 
clear p;
p(1,1) = -ustrip_width/2;
p(2,1) = 0.000;
p(1,2) = -ustrip_width/2;
p(2,2) = 500.000;
p(1,3) = ustrip_width/2;
p(2,3) = 500.000;
p(1,4) = ustrip_width/2;
p(2,4) = 0.000;
CSX = AddLinPoly( CSX, 'Metal5', 200, 'z', Metal5.zmin, p, Metal5.thick); 

% Bounding box of geometry
geometry.xmin= -15.000;
geometry.xmax= 15.000;
geometry.ymin= 0;
geometry.ymax= 500.000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end of conductors  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%  ports created manually  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[CSX, port{1}] = AddLumpedPort( CSX, 500, 1, 50, [-ustrip_width/2 geometry.ymin Metal5.zmin], [ustrip_width/2 geometry.ymin Metal1.zmax], [0 0 -1], true);
[CSX, port{2}] = AddLumpedPort( CSX, 500, 2, 50, [-ustrip_width/2 geometry.ymax Metal5.zmin], [ustrip_width/2 geometry.ymax Metal1.zmax], [0 0 -1], false);
% [CSX, port] = AddLumpedPort( CSX, prio, portnr, R, start, stop, dir, excite, varargin )
mesh.y = [linspace(geometry.ymin,geometry.ymin+converg_port_length,ceil(converg_port_length/refined_cellsize))];
mesh.x = [linspace(-ustrip_width/2, ustrip_width/2, 12)];
mesh.y = [mesh.y linspace(geometry.ymax-converg_port_length,geometry.ymax, ceil(converg_port_length/refined_cellsize))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end of ports  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geometry.width = geometry.xmax - geometry.xmin;
geometry.height = geometry.ymax - geometry.ymin;

% add margin on each side, so that metal walls have no effect
box.xmin = geometry.xmin - 2 * geometry.width;
box.xmax = geometry.xmax + 2 * geometry.width;
box.ymin = geometry.ymin - 0.2 * geometry.height;
box.ymax = geometry.ymax + 0.2 * geometry.height;


%%%%%%%%%%%%%%%%%%%%% boxes for dielectrics and substrate %%%%%%%%%%%%%%%%%%%%%%%%%%

CSX = AddBox( CSX, 'Sub', 100, [box.xmin, box.ymin, Sub.zmin], [box.xmax, box.ymax, Sub.zmax] );
CSX = AddBox( CSX, 'EPI', 100, [box.xmin, box.ymin, EPI.zmin], [box.xmax, box.ymax, EPI.zmax] );
CSX = AddBox( CSX, 'FOX', 100, [box.xmin, box.ymin, FOX.zmin], [box.xmax, box.ymax, FOX.zmax] );
CSX = AddBox( CSX, 'Sub_LI', 100, [box.xmin, box.ymin, Sub_LI.zmin], [box.xmax, box.ymax, Sub_LI.zmax] );
CSX = AddBox( CSX, 'Sub_Metal1', 100, [box.xmin, box.ymin, Sub_Metal1.zmin], [box.xmax, box.ymax, Sub_Metal1.zmax] );
CSX = AddBox( CSX, 'Sub_Metal2', 100, [box.xmin, box.ymin, Sub_Metal2.zmin], [box.xmax, box.ymax, Sub_Metal2.zmax] );
CSX = AddBox( CSX, 'Sub_Metal3', 100, [box.xmin, box.ymin, Sub_Metal3.zmin], [box.xmax, box.ymax, Sub_Metal3.zmax] );
CSX = AddBox( CSX, 'Sub_Metal4', 100, [box.xmin, box.ymin, Sub_Metal4.zmin], [box.xmax, box.ymax, Sub_Metal4.zmax] );
CSX = AddBox( CSX, 'Sub_Metal5', 100, [box.xmin, box.ymin, Sub_Metal5.zmin], [box.xmax, box.ymax, Sub_Metal5.zmax] );


%%%%%%%%%%%%%%%%%%%%%%%%%%%  build final mesh   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mesh.x = [mesh.x box.xmin box.xmax];
mesh.y = [mesh.y box.ymin box.ymax];
mesh.z = unique(sort(mesh.z));

% refine mesh in conductor regions
mesh.x = [mesh.x geometry.xmin, geometry.xmax];
mesh.x = unique(sort(mesh.x));

mesh.x = SmoothMeshLines( mesh.x, max_cellsize, 1.3, 0);
mesh.y = SmoothMeshLines( mesh.y, max_cellsize, 1.3, 0);
mesh.z = SmoothMeshLines( mesh.z, max_cellsize, 1.3, 0);

CSX = DefineRectGrid( CSX, unit, mesh );


%% write/show/run the openEMS compatible xml-file
Sim_Path = ['data/' basename];
Sim_CSX = [basename '.xml'];

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

CSXGeomPlot( [Sim_Path '/' Sim_CSX] );  % open in viewer before start
RunOpenEMS( Sim_Path, Sim_CSX ,'');  % start solver

%% post-processing
close all
f = linspace( f_start, f_stop, 1601 );
port = calcPort(port, Sim_Path, f, 'RefImpedance', 50);

s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;

%% plot results

% S11
figure 1;
plot(f/1e9,20*log10(abs(s11)),'r-','LineWidth',2);
grid on;
ylabel('|S_11| (dB)','Interpreter','None');
xlabel('frequency (GHz)');

% S21
figure 2;
plot(f/1e9,20*log10(abs(s21)),'r-','LineWidth',2);
grid on;
ylabel('|S_21| (dB)','Interpreter','None');
xlabel('frequency (GHz)');

s22 = s11;
s12 = s21;

for i=1:length(f)
  s(1,1,i) = s11(i);
  s(1,2,i) = s12(i);
  s(2,1,i) = s21(i);
  s(2,2,i) = s22(i);
  end
  
a = s2a(s);
A(:,:) = a(1,1,:);
B(:,:) = a(1,2,:);
C(:,:) = a(2,1,:);
A = transpose(A);
B = transpose(B);
C = transpose(C);
gamma = acosh(A)/((geometry.ymax-geometry.ymin)*unit);

alpha = abs(real(gamma));
alpha_dBmm = alpha*8.68/1000;
beta = unwrap(angle(s21))./((geometry.ymax-geometry.ymin)*unit);
beta_degUm = beta*1e-6*180/pi;

Q = -beta./(2*alpha);
Zc = sqrt(B./C);
Zc_im = real(Zc);

% alpha
figure 3;
plot(f/1e9,alpha_dBmm,'r-','LineWidth',2);
grid on;
ylabel('alpha (dB/mm)','Interpreter','None');
xlabel('frequency (GHz)');

% Zc
figure 4;
plot(f/1e9,Zc_im,'r-','LineWidth',2);
grid on;
ylabel('Zc (ohm)','Interpreter','None');
ylim([0 100])
xlabel('frequency (GHz)');

% Q
figure 5;
plot(f/1e9,Q,'r-','LineWidth',2);
grid on;
ylabel('Quality factor','Interpreter','None');
xlabel('frequency (GHz)');

% beta
figure 6;
plot(f/1e9,beta_degUm,'r-','LineWidth',2);
grid on;
ylabel('Propagation constant (deg/um)','Interpreter','None');
xlabel('frequency (GHz)');

save ([Sim_Path '/results']);

%% export Touchstone *.s1p
% write_s1p('s',f,s11,[Sim_Path '/' basename '.s1p'],50);
