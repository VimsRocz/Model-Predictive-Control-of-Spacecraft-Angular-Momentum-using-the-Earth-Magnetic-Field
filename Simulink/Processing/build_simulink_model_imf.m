function build_simulink_model_imf()
%BUILD_SIMULINK_MODEL_IMF Builds a Simulink model using Interpreted MATLAB Function blocks.
% xdot = tau_ext - [B]_x*m

model = "mtq_desat_mpc_model";
if bdIsLoaded(model); close_system(model,0); end
new_system(model); open_system(model);

% Basic layout positions
x0 = 40; y0 = 40; dx = 180; dy = 80;

% Pull params from base workspace
if evalin("base","exist('P','var') ~= 1")
    error("Parameter struct 'P' not found in base workspace. Run: P = params_default(); assignin('base','P',P);");
end
P = evalin("base","P");
Ts_str = num2str(P.Ts, 16);
Tend_str = num2str(P.Tend, 16);
x0_str = mat2str(P.x0(:), 16);

% Blocks
add_block("simulink/Sources/Clock", model+"/Clock", "Position",[x0 y0 x0+60 y0+30]);

add_block("simulink/User-Defined Functions/Interpreted MATLAB Function", model+"/B_field", ...
    "Position",[x0+dx y0 x0+dx+120 y0+60]);
set_param(model+"/B_field","MATLABFcn","sim_b_field","InputSignalNames","t");

add_block("simulink/User-Defined Functions/Interpreted MATLAB Function", model+"/Tau_ext", ...
    "Position",[x0+dx y0+dy x0+dx+120 y0+dy+60]);
set_param(model+"/Tau_ext","MATLABFcn","sim_tau_ext","InputSignalNames","t");

add_block("simulink/User-Defined Functions/Interpreted MATLAB Function", model+"/Controller", ...
    "Position",[x0+2*dx y0+0.5*dy x0+2*dx+140 y0+0.5*dy+90]);
set_param(model+"/Controller","MATLABFcn","sim_controller","InputSignalNames","x,B,t");

add_block("simulink/Math Operations/Gain", model+"/TsGain", ...
    "Gain",Ts_str, "Position",[x0+3*dx y0+0.5*dy x0+3*dx+60 y0+0.5*dy+50]);

add_block("simulink/Continuous/Integrator", model+"/Int_x", ...
    "Position",[x0+4*dx y0+0.5*dy x0+4*dx+60 y0+0.5*dy+60]);
set_param(model+"/Int_x","InitialCondition",x0_str);

add_block("simulink/Sinks/To Workspace", model+"/x_log", ...
    "VariableName","x_log", "SaveFormat","StructureWithTime", ...
    "Position",[x0+5*dx y0+0.5*dy x0+5*dx+100 y0+0.5*dy+60]);

add_block("simulink/Sinks/To Workspace", model+"/m_log", ...
    "VariableName","m_log", "SaveFormat","StructureWithTime", ...
    "Position",[x0+3*dx y0+1.4*dy x0+3*dx+100 y0+1.4*dy+60]);

add_block("simulink/Sinks/To Workspace", model+"/B_log", ...
    "VariableName","B_log", "SaveFormat","StructureWithTime", ...
    "Position",[x0+2*dx y0 x0+2*dx+100 y0+60]);

add_block("simulink/Sinks/To Workspace", model+"/tau_ext_log", ...
    "VariableName","tau_ext_log", "SaveFormat","StructureWithTime", ...
    "Position",[x0+2*dx y0+dy x0+2*dx+100 y0+dy+60]);

add_block("simulink/Sinks/To Workspace", model+"/tau_mtq_log", ...
    "VariableName","tau_mtq_log", "SaveFormat","StructureWithTime", ...
    "Position",[x0+3.7*dx y0+1.4*dy x0+3.7*dx+110 y0+1.4*dy+60]);

add_block("simulink/Sinks/To Workspace", model+"/tau_total_log", ...
    "VariableName","tau_total_log", "SaveFormat","StructureWithTime", ...
    "Position",[x0+4.6*dx y0+1.4*dy x0+4.6*dx+110 y0+1.4*dy+60]);

% Sum: tau_ext + tau_mtq
add_block("simulink/Math Operations/Sum", model+"/SumTau", ...
    "Inputs","++", "Position",[x0+3*dx y0+0.5*dy-40 x0+3*dx+40 y0+0.5*dy]);

% tau_mtq = cross(m, B)
add_block("simulink/Math Operations/Cross Product", model+"/MTQ_Torque", ...
    "Position",[x0+2.5*dx y0+1.4*dy x0+2.5*dx+140 y0+1.4*dy+90]);

% Wire connections
add_line(model,"Clock/1","B_field/1");
add_line(model,"Clock/1","Tau_ext/1");

add_line(model,"Int_x/1","Controller/1");     % x -> controller
add_line(model,"B_field/1","Controller/2");   % B -> controller
add_line(model,"Clock/1","Controller/3");     % t -> controller

add_line(model,"B_field/1","B_log/1");        % log B
add_line(model,"Tau_ext/1","tau_ext_log/1");  % log tau_ext

add_line(model,"Controller/1","m_log/1");     % log m
add_line(model,"Controller/1","MTQ_Torque/1");% m -> mtq torque
add_line(model,"B_field/1","MTQ_Torque/2");   % B -> mtq torque

add_line(model,"Tau_ext/1","SumTau/1");       % tau_ext
add_line(model,"MTQ_Torque/1","SumTau/2");    % tau_mtq
add_line(model,"MTQ_Torque/1","tau_mtq_log/1"); % log tau_mtq

add_line(model,"SumTau/1","TsGain/1");
add_line(model,"SumTau/1","tau_total_log/1"); % log total torque
add_line(model,"TsGain/1","Int_x/1");
add_line(model,"Int_x/1","x_log/1");

% Set model variables
set_param(model,"StopTime",Tend_str);
set_param(model,"Solver","ode4","FixedStep",Ts_str);

save_system(model);
disp("Built Simulink model: " + model + ".slx");
end
