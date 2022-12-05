%% Clear Workspace and define pathnames

clearvars; clc; close all;
Starting_Directory=pwd;

% Add MSOT MAtlab library to path 
% This library was created by iThera
addpath( genpath('.') )
javaaddpath .\MSOTBeans\xbean.jar %%AK
javaaddpath .\MSOTBeans\msotbeans.jar %%AK

%% Inputs
prompt = {'Do you want to analyze the OE MSOT data? (answer: Y/N) ',...
    'Do you want to analyze the DCE MSOT MSP data? (answer: Y/N)'};
dlgtitle = 'Input';
dims = [1 100];
definput = {'Y','N'};
answer = inputdlg(prompt,dlgtitle,dims,definput,'on');

oe_msot_answer = answer{1,1};
msp_msot_answer = answer{2,1};

%% Identify the folder paths

if strcmpi(oe_msot_answer,'Y')
    disp('------------------------------');
    disp('OE - MSOT');
    disp('------------------------------');
    disp('Select the directory containing the msot file for Medical Air data ');
    DIR=uigetdir;
    cd(DIR)
    msot_name = dir ('*.msot');
    msot_file_medair= fullfile(DIR, msot_name.name);

    datainfo = loadMSOT(msot_file_medair);
    X = ['The selected file is: ',datainfo.Name];
    disp(X);
    prompt = 'Do you want to continue? Y/N [Y]: ';
    str = input(prompt,'s');
    if str == 'N'
        return;
    end

    disp('Select the directory containing the msot file for 100% 02 data ');
    DIR=uigetdir;
    cd(DIR)
    msot_name = dir ('*.msot');
    msot_file_02 = fullfile(DIR, msot_name.name);

    datainfo = loadMSOT(msot_file_02);
    X = ['The selected file is: ',datainfo.Name];
    disp(X);
    prompt = 'Do you want to continue? Y/N [Y]: ';
    str = input(prompt,'s');
    if str == 'N'
        return;
    end

    disp('Select the directory for saving the final results ');
    saving_DIR=uigetdir; 
    
    segm = run_OE_MSOT_AK_nonneg(msot_file_medair, msot_file_02, saving_DIR);
end

if strcmpi(msp_msot_answer,'Y')
    disp('------------------------------');
    disp('DCE - MSOT');
    disp('------------------------------');
    
    disp('Select the directory containing the msot file MSP data ');
    DIR=uigetdir;
    cd(DIR)
    msot_name = dir ('*.msot');
    msot_file= fullfile(DIR, msot_name.name);

    datainfo = loadMSOT(msot_file);
    X = ['The selected file is: ',datainfo.Name];
    disp(X);
    prompt = 'Do you want to continue? Y/N [Y]: ';
    str = input(prompt,'s');
    if str == 'N'
        return;
    end

    disp('Select the directory for saving the final results ');
    saving_DIR=uigetdir; 
    
    if oe_msot_answer == 'Y'
        run_DCE_MSOT_MSP (msot_file, saving_DIR, segm);
    else
        run_DCE_MSOT_MSP (msot_file, saving_DIR);
    end

end
    