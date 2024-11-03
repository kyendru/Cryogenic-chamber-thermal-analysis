clc,clear,close all


%% Thermal properties of a hollow cylinder where the inner and outer surfaces are held at a constant temperature
%Defining geometry
model = createpde('thermal');
geo = multicylinder([20 25 35],20,'Void',[1 0 0]);
model.Geometry = geo;
figure
pdegplot(model, 'FaceLabels','on')
generateMesh(model);
title("Visualizing the Geometry")
figure
pdemesh(model)
title("Visualizing the Mesh")

%Defining thermal properties and boundary conditions of the solid
thermalProperties(model,'Cell',1,'ThermalConductivity',40);
thermalProperties(model,'Cell',2,'ThermalConductivity',0.15);
thermalBC(model,'Face',3,'Temperature',85);
thermalBC(model,'Face',7,'Temperature',4); %add radiation as BC

%Solving the problem
result = solve(model);
figure
pdeplot3D(model,'ColorMapData',result.Temperature)
title("Thermal Analysis Results")


%% Robot arm analysis exposed to heat through conduction (thermal)

%Importing geometry
thermalmodel = createpde('thermal','steadystate');
importGeometry(thermalmodel,'Gripper Pivot.stl'); 
figure
pdegplot(thermalmodel, 'FaceAlpha',0.5)
title("Visualizing the Geometry")

%Generating mesh
generateMesh(thermalmodel,'Hmax',0.09);
thermalmodel.Mesh
figure
pdemesh(thermalmodel)
title("Visualizing the Mesh")

%Measuring mesh dimensions
xLength = max(thermalmodel.Mesh.Nodes(1,:)) - min(thermalmodel.Mesh.Nodes(1,:));
yLength = max(thermalmodel.Mesh.Nodes(2,:)) - min(thermalmodel.Mesh.Nodes(2,:));
zLength = max(thermalmodel.Mesh.Nodes(3,:)) - min(thermalmodel.Mesh.Nodes(3,:));

%Defining material properties: magnesium alloy
u = symunit; 
tc = vpa(52) *u.Watt/(u.m*u.Kelvin);
tc = rewrite(tc,u.W/(u.in*u.Kelvin));
tc = double(separateUnits(tc));
thermalProperties(thermalmodel,'Cell',1,'ThermalConductivity',tc); 

%Setting up boundary conditions
htc = vpa(25) *u.Watt/(u.m^2*u.Kelvin);
htc = rewrite(htc,u.Watt/(u.in^2*u.Kelvin));
htc = double(separateUnits(htc));
thermalBC(thermalmodel,'Face',1:thermalmodel.Geometry.NumFaces,...
        'ConvectionCoefficient',htc,'AmbientTemperature',288.15);

figure
pdegplot(thermalmodel,'FaceLabels','on', 'FaceAlpha',0.5)

thermalBC(thermalmodel,'Face',43,'Temperature',288.15); %robotic arm is acting as a heat sink
thermalBC(thermalmodel,'Face',38,'HeatFlux',25);
thermalBC(thermalmodel,'Face',44,'HeatFlux',10);

%Solving the problem
tic
results = solve(thermalmodel);
time = toc;

%Computing max temperature
t_max = max(results.Temperature);
vpa(rewrite(t_max*u.Kelvin,u.Celsius,'Temperature','absolute'),5)
figure
vpa(rewrite(t_max*u.Kelvin,u.Fahrenheit,'Temperature','absolute'),5)
pdeplot3D(thermalmodel,'ColorMapData',results.Temperature)
title('Temperature in K')
colorbar
%Highest temperature on the part closest to ciruit board, heat is dissipated
%as it moves along the part.

%Running multiple analysis by testing different materials
thermalConductivities = [45.3 51.2 63.2 76.9 83.5]; 
runtime = time * length(thermalConductivities);

%Running a parallel for loop to speed up analysis
maxTemperatures = zeros(size(thermalConductivities)); 
thermalConductivities_inches = double(separateUnits(...
    rewrite(thermalConductivities *u.Watt/(u.m*u.Kelvin),u.W/(u.in*u.Kelvin)))); 
tic
parfor i = 1: length(thermalConductivities)
    thermalProperties(thermalmodel,'Cell',1,...
    'ThermalConductivity',thermalConductivities_inches(i)); 
    results = solve(thermalmodel);
    maxTemperatures(i) = max(results.Temperature)
end
tF = toc;

figure
plot(thermalConductivities,maxTemperatures,'o')
xlabel('Thermal Conductivity $\frac{W}{mK}$','Interpreter','latex')
ylabel('Max Temperature $^\circ$K','Interpreter','latex')
title('Max Temperature v. Thermal Conductivity')
%As thermal conductivity increases, max temperature decreases.

%Creating predictive model
linearmodel = fitlm(thermalConductivities(:),maxTemperatures(:),'VarNames',...
    {'ThermalConductivity', 'MaxTemperature'})
figure
plot(linearmodel) 
%Can be seen that the linear model is a good fit.

%Finding associated allowable heat transfer coefficient for the alloy
shutOffTemp = 311;
objective = @(x) predict(linearmodel,x) - shutOffTemp;
max_tc = fzero(objective,50);

%% Deflection of a bracket (structural)

%Importing geometry
model = createpde('structural','static-solid');
importGeometry(model,'BracketWithHole.stl');
figure
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Bracket with Face Labels')

%Defining structural properties
structuralProperties(model,'Cell',1,'YoungsModulus',200e9, ...
                                    'PoissonsRatio',0.3);

%Defining boundary conditions
structuralBC(model,'Face',4,'Constraint','fixed'); %Cantilever
distributedLoad = 1e4; % Applied load in Pascals
structuralBoundaryLoad (model,'Face',8,'SurfaceTraction',[0;0;-distributedLoad]);

%Creating a mesh
bracketThickness = 1e-2; % Thickness of horizontal plate with hole, meters
generateMesh(model,'Hmax',bracketThickness);
figure
pdeplot3D(model)
title('Visualizing the Mesh');

%Solving the problem
result = solve(model);

%Finding maximum deflection and plotting displacement components
minUz = min(result.Displacement.uz);
fprintf('Maximal deflection in the z-direction is %g meters.', minUz)

figure
pdeplot3D(model,'ColorMapData',result.Displacement.ux)
title('x-displacement')
colormap('jet')

figure
pdeplot3D(model,'ColorMapData',result.Displacement.uy)
title('y-displacement')

figure
pdeplot3D(model,'ColorMapData',result.Displacement.uz)
title('z-displacement')

%Plotting von Mises stress ( used to predict yielding of materials under
%complex loading from the results of uniaxial tensile tests)
figure
pdeplot3D(model,'ColorMapData',result.VonMisesStress,'Deformation',result.Displacement,'DeformationScaleFactor',1000)
title('von Mises stress')

%Analysing multiple materials to find material that results in least delfection

%Importing material properties
gcp;
material_parameters = readtable('Materials_ParametersList.txt','ReadVariableNames',false,'HeaderLines',1);
material_parameters.Properties.VariableNames = {'Material','YoungsModulus','PoissonsRatio'};
material_parameters.Properties.VariableUnits = {'','GPa',''} ;

idx_max = size(material_parameters,1); % Number of materials to consider
distributedLoad = 1e4; % Applied load in Pascals
bracketThickness = 1e-2; % Thickness of horizontal plate with hole, meters

Model = cellfun(@(~) createpde('structural','static-solid'),material_parameters.Material,'UniformOutput',false);
YoungMV = material_parameters.YoungsModulus*1e9; % GPa 
PoissR = material_parameters.PoissonsRatio;

%Creating structural analysis models
pb = cell(1,idx_max);
results = cell(1,idx_max);

%Solving the problem
parfor ii = 1:idx_max
    pb{ii} = setupFEA(Model{ii},YoungMV(ii),PoissR(ii),distributedLoad);
    results{ii} = solve(pb{ii});    
end

%Finding minimum deflection
minUz = cellfun(@(x) min(x.Displacement.uz),results);
[mUz,Imin] = max(minUz);
fprintf(['The material with the minimum deflection in the z-direction is '...
    material_parameters.Material{Imin} ' with ' num2str(mUz) ' m']);
material_parameters.Max_Deflection = abs(minUz)'; % Append deflection values to our table
%Steel has the least deflection, aluminum has the most

%Plotting Z-displacement since maximum deflecion is in the z direction
for jj = 1:idx_max
    figure
    subplot(2,3,jj)
    h = pdeplot3D(pb{jj},'ColorMapData',results{jj}.Displacement.uz,...
        'Deformation',results{jj}.Displacement,'DeformationScaleFactor',5e2);
    cm = colormap(gca,jet);
    colormap(gcf,flip(cm));
    clim(gca,[min(minUz) 0])
    title(material_parameters.Material(jj))
end

%Computing load deflection curves for all materialsfor 10 different loads in the range given below
F_Load = linspace(1e2,1e8,10);
F_Load_size = size(F_Load,2);
[X,Y] = meshgrid(1:idx_max,F_Load); %to iterate ona single loop of all material parameter and load combinations
lD_result = nan(F_Load_size,idx_max);

parfor ii = 1:numel(X)
   mdl = pb{X(ii)};
   lD_result(ii) = solveLD(mdl,Y(ii)); 
end

%Visualizing curves
figure
plot(lD_result,F_Load,'-.x')
xlabel('Displacement [m]')
ylabel('Load [Pa]')
legend(material_parameters.Material,'Location','southeast')
grid on;


%% Tuning Fork Analysis (structural)

%Importing geometry
model = createpde('structural','modal-solid');
importGeometry(model,'TuningFork.stl');
figure
pdegplot(model)

%Defining structural properties
E = 210E9;
nu = 0.3;
rho = 8000;
structuralProperties(model,'YoungsModulus',E, ...
                           'PoissonsRatio',nu, ...
                           'MassDensity',rho);

%Generating mesh
generateMesh(model,'Hmax',0.001);

%Solving model
RF = solve(model,'FrequencyRange',[-1,4000]*2*pi);

%Visualizing results
modeID = 1:numel(RF.NaturalFrequencies);
tmodalResults = table(modeID.',RF.NaturalFrequencies/2/pi);
tmodalResults.Properties.VariableNames = {'Mode','Frequency'};
disp(tmodalResults);
%First six modes of the tuning fork are due to rigid body motion, seventh
%mode is the first flexible mode

frames  = animateSixTuningForkModes(RF);

%Simulating fork dynamics as is it quickly struck
tmodel = createpde('structural','transient-solid');
importGeometry(tmodel,'TuningFork.stl');
mesh = generateMesh(tmodel,'Hmax',0.005);
structuralProperties(tmodel,'YoungsModulus',E, ...
                            'PoissonsRatio',nu, ...
                            'MassDensity',rho);
figure('units','normalized','outerposition',[0 0 1 1])
pdegplot(tmodel,'FaceLabels','on')
view(-50,15)
title 'Geometry with Face Labels'

%Defining boundary conditions
structuralBC(tmodel,'Face',[21,22],'Constraint','fixed'); %COnstraining face 21 and 22
T = 2*pi/RF.NaturalFrequencies(7); %pressure applied only for a fraction of the time period of the fundamental mode
structuralBoundaryLoad(tmodel,'Face',11,'Pressure',5E6,'EndTime',T/300);
structuralIC(tmodel,'Displacement',[0;0;0],'Velocity',[0;0;0]);

%Solving the problem
ncycle = 50;
samplingFrequency = 60/T;
tlist = linspace(0,ncycle*T,ncycle*T*samplingFrequency);
R = solve(tmodel,tlist);

excitedTineTipNodes = findNodes(mesh,'region','Face',12);
tipDisp = R.Displacement.uy(excitedTineTipNodes(1),:);

figure
plot(R.SolutionTimes,tipDisp)
title('Transverse Displacement at Tine Tip')
xlim([0,0.1])
xlabel('Time')
ylabel('Y-Displacement')

[fTip,PTip] = tuningForkFFT(tipDisp,samplingFrequency);
figure
TF = islocalmax(PTip,'MaxNumExtrema',6);
plot(fTip,PTip,fTip(TF),PTip(TF),'r*') 
title({'Single-sided Amplitude Spectrum', 'of Tip Vibration'})
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0,4000])

baseNodes = tmodel.Mesh.findNodes('region','Face',6);
baseDisp = R.Displacement.ux(baseNodes(1),:);
figure
plot(R.SolutionTimes,baseDisp)
title('Axial Displacement at the End of Handle')
xlim([0,0.1])
ylabel('X-Displacement')
xlabel('Time')

[fBase,PBase] = tuningForkFFT(baseDisp,samplingFrequency);
figure
TFb = islocalmax(PBase,'MaxNumExtrema',6);
plot(fBase,PBase,fBase(TFb),PBase(TFb),'r*') 
title({'Single-sided Amplitude Spectrum', 'of Base Vibration'})
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0,4000])

Function to solve load deflection models
function modelN = setupFEA(modelN,YoungM,PoissonR,distrL)
    modelN.Geometry = importGeometry(modelN,'BracketWithHole.stl');
    structuralProperties(modelN,'Cell',1,'YoungsModulus',YoungM,'PoissonsRatio',PoissonR);
    structuralBC(modelN,'Face',4,'Constraint','fixed');
    structuralBoundaryLoad (modelN,'Face',8,'SurfaceTraction',[0;0;-distrL]);
    bracketThickness = 1e-2; % Thickness of horizontal plate with hole, meters
    generateMesh(modelN,'Hmax',bracketThickness);
end

function ld_mat = solveLD(modelM,ldm)
structuralBoundaryLoad (modelM,'Face',8,'SurfaceTraction',[0;0;-ldm]);
ld_result = solve(modelM);
ld_mat = abs(min(ld_result.Displacement.uz));
end