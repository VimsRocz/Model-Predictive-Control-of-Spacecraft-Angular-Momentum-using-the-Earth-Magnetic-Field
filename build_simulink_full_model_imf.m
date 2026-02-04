function build_simulink_full_model_imf()
%BUILD_SIMULINK_FULL_MODEL_IMF Build Simulink model for full plant (IMF).
%
% Uses a MATLAB Function block that calls sim_full_step(t) via
% coder.extrinsic, so all MATLAB toolboxes (IGRF, etc.) are supported.

model = "mtq_full_model";
if bdIsLoaded(model); close_system(model,0); end
new_system(model); open_system(model);

% Layout
x0 = 40; y0 = 40; dx = 220; dy = 70;
nl = newline;

% Params from base workspace
if evalin("base","exist('P','var') ~= 1")
    error("Parameter struct 'P' not found in base workspace. Run: P = params_default(); assignin('base','P',P);");
end
P = evalin("base","P");
Ts_str = num2str(P.Ts, 16);
Tend_str = num2str(P.Tend, 16);

% Blocks
add_block("simulink/Sources/Clock", model+"/Clock", "Position",[x0 y0 x0+60 y0+30]);

add_block("simulink/User-Defined Functions/MATLAB Function", model+"/FullStep", ...
    "Position",[x0+dx y0 x0+dx+180 y0+120]);

% To Workspace blocks
outs = { ...
    "r_eci_log", "B_eci_log", "B_body_log", "h_w_log", "m_log", ...
    "tau_ext_log", "tau_mtq_log", "tau_rw_log", "q_ib_log", "q_ib_des_log", ...
    "w_body_log", "w_des_body_log" };

for i = 1:numel(outs)
    add_block("simulink/Sinks/To Workspace", model + "/" + outs{i}, ...
        "VariableName", outs{i}, "SaveFormat", "StructureWithTime", ...
        "Position",[x0+2*dx y0+(i-1)*dy x0+2*dx+120 y0+(i-1)*dy+30]);
end

% MATLAB Function block script
set_mfcn_script(model+"/FullStep", ...
"function [r_eci, B_eci, B_body, h_w, m, tau_ext, tau_mtq, tau_rw, q_ib, q_ib_des, w_body, w_des_body] = fcn(t)" + nl + ...
"%#codegen" + nl + ...
"coder.extrinsic('sim_full_step');" + nl + ...
"r_eci = zeros(3,1); B_eci = zeros(3,1); B_body = zeros(3,1);" + nl + ...
"h_w = zeros(3,1); m = zeros(3,1); tau_ext = zeros(3,1);" + nl + ...
"tau_mtq = zeros(3,1); tau_rw = zeros(3,1);" + nl + ...
"q_ib = zeros(1,4); q_ib_des = zeros(1,4);" + nl + ...
"w_body = zeros(3,1); w_des_body = zeros(3,1);" + nl + ...
"[r_eci, B_eci, B_body, h_w, m, tau_ext, tau_mtq, tau_rw, q_ib, q_ib_des, w_body, w_des_body] = sim_full_step(t);" + nl + ...
"end" + nl);

% Output sizes
set_mfcn_output_size(model+"/FullStep", "r_eci", "3");
set_mfcn_output_size(model+"/FullStep", "B_eci", "3");
set_mfcn_output_size(model+"/FullStep", "B_body", "3");
set_mfcn_output_size(model+"/FullStep", "h_w", "3");
set_mfcn_output_size(model+"/FullStep", "m", "3");
set_mfcn_output_size(model+"/FullStep", "tau_ext", "3");
set_mfcn_output_size(model+"/FullStep", "tau_mtq", "3");
set_mfcn_output_size(model+"/FullStep", "tau_rw", "3");
set_mfcn_output_size(model+"/FullStep", "q_ib", "4");
set_mfcn_output_size(model+"/FullStep", "q_ib_des", "4");
set_mfcn_output_size(model+"/FullStep", "w_body", "3");
set_mfcn_output_size(model+"/FullStep", "w_des_body", "3");
set_mfcn_input_size(model+"/FullStep", "t", "1");

% Wiring
add_line(model, "Clock/1", "FullStep/1");

for i = 1:numel(outs)
    add_line(model, "FullStep/" + string(i), outs{i} + "/1");
end

% Model solver settings
set_param(model, "Solver", "ode4", "FixedStep", Ts_str);
set_param(model, "StopTime", Tend_str);

save_system(model);
disp("Built Simulink full model: " + model + ".slx");
end

function set_mfcn_script(blockPath, scriptText)
rt = sfroot;
ch = rt.find('-isa','Stateflow.EMChart','Path', blockPath);
if isempty(ch)
    error("Could not find Stateflow chart for block: %s", blockPath);
end
ch.Script = char(scriptText);
end

function set_mfcn_output_size(blockPath, outputName, sizeSpec)
rt = sfroot;
ch = rt.find('-isa','Stateflow.EMChart','Path', blockPath);
if isempty(ch)
    error("Could not find Stateflow chart for block: %s", blockPath);
end
outData = ch.find('-isa','Stateflow.Data','Scope','Output');
for i = 1:numel(outData)
    if strcmp(outData(i).Name, outputName)
        outData(i).Props.Array.Size = sizeSpec;
    end
end
end

function set_mfcn_input_size(blockPath, inputName, sizeSpec)
rt = sfroot;
ch = rt.find('-isa','Stateflow.EMChart','Path', blockPath);
if isempty(ch)
    error("Could not find Stateflow chart for block: %s", blockPath);
end
inData = ch.find('-isa','Stateflow.Data','Scope','Input');
for i = 1:numel(inData)
    if strcmp(inData(i).Name, inputName)
        inData(i).Props.Array.Size = sizeSpec;
    end
end
end

