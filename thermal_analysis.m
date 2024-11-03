clc, clear, close all
set(0,'DefaultFigureWindowStyle','docked')

%% Initialization
emissivity = 1;             %emissivity
thermal_cond = 237;         %thermal conductivity of aluminum in W/(mK)
spec_heat = 910;            %specific heat of aluminum in J/kg K
mass_den = 2710;            %density of aluminum in kg/m^3
inch_to_m = 0.0254;              %conversion factor for inch to m 
l_outer = 10*inch_to_m;          %length of cold shroud in m
thickness = 1*inch_to_m;       %thickness of cold shroud
l_inner = l_outer-(2*thickness); %inner length of cold shroud in m
T_cham = 300;                    %chamber temperature in Kelvin
T_cs = 70;                       %temperature of bottom surface of cold shroud
stefan_boltzmann_constant = 5.670373E-8; 
L = 1.992e5; %J/kg
A = (l_outer^2-l_inner^2);  %area of bottom surface in m^2
volume_cube = l_outer^3 - l_inner^3 - (thickness*l_inner^2); %volume of cube in m^3
Q = mass_den*volume_cube*spec_heat*(T_cham-T_cs);            %total heat to be removed at first time step
mass_theoretical = Q/L;          

%% Creating the thermal model, geometry and mesh
thermalmodel = createpde("thermal","transient");
gm = importGeometry(thermalmodel,"cold_shroud v4.stl");

%Visualizing faces
figure
pdegplot(thermalmodel, 'FaceLabels','on','EdgeLabels','on', 'VertexLabels', 'on')
title("Visualizing the Geometry")

outerFaces = [6 7 8 9 10];    %outer faces(except bottom)
innerFaces = [1 2 3 4 5 11];  %inner faces
botttomFace = 11;
topFace = 10;

%Generating mesh
mesh = generateMesh(thermalmodel);
figure
pdemesh(thermalmodel)

%% Solving the thermal model
%Thermal properties of aluminum
thermalProperties(thermalmodel, 'ThermalConductivity', thermal_cond, ...
                                'MassDensity', mass_den, ...
                                'SpecificHeat', spec_heat);

%Boundary conditions
thermalIC(thermalmodel, T_cham);  %initial temperature set to 300K
thermalBC(thermalmodel,"Face",11,"Temperature",T_cs);  %lower face set to 70K

%Heat flux function 
heatFluxFunc = @(region, state) -emissivity * stefan_boltzmann_constant * ...
                                (state.u.^4 - T_cham^4);

%Case 1: no heat flux
%thermalBC(thermalmodel,"Face",10,"Temperature",300);   

%Case 2: heat flux applied to one side only
%thermalBC(thermalmodel, 'Face', 6, 'HeatFlux', heatFluxFunc);

%Case 3: heat flux to two adjacent sides
%thermalBC(thermalmodel, 'Face', [6,7], 'HeatFlux', heatFluxFunc);

%Case 4: heat flux to two opposite sides
%thermalBC(thermalmodel, 'Face', [6,8], 'HeatFlux', heatFluxFunc);

%Case 5: heat flux to all side faces
%thermalBC(thermalmodel, 'Face', [6,7,8,9], 'HeatFlux', heatFluxFunc);

%Case 6: heat flux to top face only
%thermalBC(thermalmodel, 'Face', 10, 'HeatFlux', heatFluxFunc);

%Case 7: heat flux to all outer faces
thermalBC(thermalmodel, 'Face', outerFaces, 'HeatFlux', heatFluxFunc);

t = mass_den*spec_heat*l_outer^2/thermal_cond;  %time required to reach thermal equilibrium
%tspan = 0:1e2:1e4;
tspan = 0:1e2:604800;
%time steps for analysis
%tspan = 0:10:t;

%Solving the thermal problem
result = solve(thermalmodel,tspan);

% Plotting the result
figure
pdeplot3D(thermalmodel,'ColorMapData',result.Temperature(:,end));  

%% Calculating total heat 
bottomface_nodes = findNodes(mesh, 'region', 'Face', 11); %finding nodes on bottom surface
bottomface_coords = mesh.Nodes(:, bottomface_nodes)';     %finding coordinates of nodes on bottom surface

%Finding temperature at each node at each time step at the bottom surface
%and 0.01m above the bottom surface
for j = 1:length(tspan)
    for i = 1:length(bottomface_nodes)
        temp1(i,j) = interpolateTemperature(result, bottomface_coords(i,1), bottomface_coords(i,2), bottomface_coords(i,3), j);
        temp2(i,j) = interpolateTemperature(result, bottomface_coords(i,1), bottomface_coords(i,2), 0.01, j);
    end
end

%Finding mass flow rate at each time step
for in = 1:1:length(tspan)
    dT = temp2(:,in)-temp1(:,in);   %finding delta T
    dz = 0.01;                      %finding delta Z
    temp_grad = dT./dz;             %finding the temperature gradient in K/m
    heat_flux_den = -thermal_cond*temp_grad;       %finding heat flux density using Fourier's Law in W/m^2
    dA = ones(size(heat_flux_den)) * A / numel(heat_flux_den); %size of each surface area element on bottom surface
    heat_flow_rate = sum(heat_flux_den .* dA);     %heat flow rate in W (J/s)
    mdot(in) = abs(heat_flow_rate)/L;              %kg/s
end

mdot_min = mdot.*(60);  %finding heat flow rate in kg/min
 
figure()
semilogy(tspan./(60),mdot_min)
grid on
xlabel("Time in days")
ylabel("Log of Mass flow rate vs Time")

% mass_total(1) = 0;
% for i = 2:length(mdot)
%     mass_total(i) = mass_total(i-1)+ mdot(i-1).*100;
% end
% 
% figure()
% plot(tspan, mass_total)
% xlabel("Time in seconds")
% ylabel("Mass in kg")
% grid on 

% Calculate the cumulative mass using trapezoidal integration
mass = cumtrapz(tspan,mdot);

% Plot the results
figure()
plot(tspan./(24*60*60), mass)
grid on
xlabel('Time (s)')
ylabel('Mass (kg)')
title('Mass as a Function of Time')

%% References
%mass density:https://www.thyssenkrupp-materials.co.uk/density-of-aluminium.html#:~:text=The%20density%20of%20aluminium%20is,m3%20and%202%2C810kg%2Fm3.
%specific heat:https://www.engineeringtoolbox.com/specific-heat-metals-d_152.html
%thermal conductivity: https://markhammetals.com/copper-vs-aluminum-which-is-the-better-conductor-of-heat/#:~:text=Aluminum%20thermal%20conductivity%20is%20about,utensils%20and%20HVAC%20system%20production.
%STL scale changing function: https://www.mathworks.com/matlabcentral/fileexchange/164836-stl-scale-changer
%Fourier's law: https://physics.emory.edu/faculty/brody/Advanced%20Lab/phys%20222%20lecture%20notes.pdf
%Rate of heat flow: https://www.mathworks.com/matlabcentral/answers/403664-the-rate-of-heat-flow-conduction-between-two-points#:~:text=The%20rate%20of%20heat%20flow%20(conduction)%20between%20two%20points%20on,distance%20from%20the%20heated%20end.