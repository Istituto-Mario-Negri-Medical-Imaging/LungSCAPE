%% setup_matlab.m - LungScape MATLAB Environment Setup & Verification
%
% Purpose:
%   Checks that all MATLAB toolboxes and external libraries required by
%   LungScape are installed and accessible on the MATLAB path.
%   Run this script once after installation to verify your environment
%   before running LungScape.m.
%
% Usage:
%   cd /path/to/LungSCAPE_codebase
%   setup_matlab
%
% What is checked:
%   1. MATLAB toolboxes: Image Processing, Parallel Computing
%   2. External libraries: NIfTI tools, anisodiff3D, vesselness3D
%      (including MEX compilation of eig3volume.c), Image Graphs
%
% If a library is missing:
%   - The download URL and installation instructions are printed.
%   - For vesselness3D, instructions for MEX compilation are provided.
%
% After resolving all issues, re-run this script to confirm.
%
% See also: LungScape, addpath, mex

clc;
fprintf('=================================================================\n');
fprintf('  LungScape - MATLAB Environment Setup & Verification\n');
fprintf('=================================================================\n\n');

issues = {};   % accumulate problem descriptions for final summary

%% -----------------------------------------------------------------------
%  1. MATLAB Toolboxes
%% -----------------------------------------------------------------------
fprintf('--- MATLAB Toolboxes ---\n\n');

tbStatus = checkToolbox('Image Processing Toolbox',    'image_toolbox');
if ~tbStatus
    issues{end+1} = 'Image Processing Toolbox is not available (required for morphological operations and image filtering).';
end

tbStatus = checkToolbox('Parallel Computing Toolbox',  'distrib_computing_toolbox');
if ~tbStatus
    issues{end+1} = 'Parallel Computing Toolbox is not available (required for parfor acceleration).';
end

%% -----------------------------------------------------------------------
%  2. External Libraries
%% -----------------------------------------------------------------------
fprintf('\n--- External Libraries ---\n\n');

% --- NIfTI and ANALYZE tools ---
ok = checkFunction('load_untouch_nii', 'NIfTI and ANALYZE tools');
if ~ok
    fprintf(['  -> Download from: https://mathworks.com/matlabcentral/fileexchange/8797\n' ...
             '     Unzip into a folder and add it to the MATLAB path:\n' ...
             '       addpath(genpath(''/path/to/NIfTI_tools''))\n\n']);
    issues{end+1} = 'NIfTI and ANALYZE tools not found (load_untouch_nii).';
end

% --- anisodiff3D ---
ok = checkFunction('anisodiff3D', 'anisodiff3D');
if ~ok
    fprintf(['  -> Download from: https://mathworks.com/matlabcentral/fileexchange/14995\n' ...
             '     Unzip into a folder and add it to the MATLAB path:\n' ...
             '       addpath(genpath(''/path/to/anisodiff3D''))\n\n']);
    issues{end+1} = 'anisodiff3D not found.';
end

% --- vesselness3D (Jerman) ---
%     Requires both the .m function AND the compiled MEX file eig3volume.
%     The MEX file must be compiled from eig3volume.c included in the library.
ok_m   = exist('vesselness3D',  'file') > 0;
ok_mex = exist('eig3volume',    'file') == 3;   % 3 = MEX file

if ok_m && ok_mex
    fprintf('  [OK] vesselness3D  (vesselness3D.m + eig3volume MEX)\n\n');
elseif ok_m && ~ok_mex
    fprintf('  [!!] vesselness3D.m found but eig3volume MEX is NOT compiled.\n');
    fprintf(['  -> In MATLAB, navigate to the vesselness3D folder and run:\n' ...
             '       mex eig3volume.c\n' ...
             '     A C compiler must be configured (run ''mex -setup'' if needed).\n\n']);
    issues{end+1} = 'vesselness3D: eig3volume MEX file not compiled (run mex eig3volume.c in the library folder).';
elseif ~ok_m && ok_mex
    fprintf('  [!!] eig3volume MEX found but vesselness3D.m is missing.\n');
    fprintf(['  -> Download from: https://mathworks.com/matlabcentral/fileexchange/63171\n' ...
             '     Unzip and add to the MATLAB path:\n' ...
             '       addpath(genpath(''/path/to/vesselness3D''))\n' ...
             '     Then compile: mex eig3volume.c\n\n']);
    issues{end+1} = 'vesselness3D: vesselness3D.m not found.';
else
    fprintf('  [!!] vesselness3D  NOT FOUND\n');
    fprintf(['  -> Download from: https://mathworks.com/matlabcentral/fileexchange/63171\n' ...
             '     Unzip and add to the MATLAB path:\n' ...
             '       addpath(genpath(''/path/to/vesselness3D''))\n' ...
             '     Then compile the required MEX file:\n' ...
             '       cd /path/to/vesselness3D\n' ...
             '       mex eig3volume.c\n' ...
             '     A C compiler must be configured (run ''mex -setup'' if needed).\n\n']);
    issues{end+1} = 'vesselness3D not found (vesselness3D.m + eig3volume MEX both missing).';
end

% --- Image Graphs ---
ok = checkFunction('binaryImageGraph', 'Image Graphs');
if ~ok
    fprintf(['  -> Download from: https://mathworks.com/matlabcentral/fileexchange/53614\n' ...
             '     Unzip into a folder and add it to the MATLAB path:\n' ...
             '       addpath(genpath(''/path/to/ImageGraphs''))\n' ...
             '     Note: binaryImageGraph3Weighted.m and binaryImageGraphWeighted.m\n' ...
             '     in LungScape/utils/ are extensions of this library and require it.\n\n']);
    issues{end+1} = 'Image Graphs not found (binaryImageGraph).';
end

%% -----------------------------------------------------------------------
%  3. Final Summary
%% -----------------------------------------------------------------------
fprintf('\n=================================================================\n');
if isempty(issues)
    fprintf('  All dependencies found. LungScape is ready to run.\n');
    fprintf('  -> Open LungScape.m, set DATA_DIRECTORY, and run LungScape.\n');
else
    fprintf('  %d issue(s) found. Please resolve the following:\n\n', numel(issues));
    for k = 1:numel(issues)
        fprintf('  %d. %s\n', k, issues{k});
    end
    fprintf('\n  After installing missing libraries, re-run setup_matlab.\n');
end
fprintf('=================================================================\n');


%% -----------------------------------------------------------------------
%  Helper functions
%% -----------------------------------------------------------------------

function ok = checkToolbox(name, licenseKey)
%CHECKTOOLBOX Check if a MATLAB toolbox is licensed and available.
    ok = license('test', licenseKey);
    if ok
        fprintf('  [OK] %s\n\n', name);
    else
        fprintf('  [!!] %s  NOT AVAILABLE\n', name);
        fprintf('  -> A valid MATLAB license including this toolbox is required.\n\n');
    end
end

function ok = checkFunction(funcName, libName)
%CHECKFUNCTION Check if a function (from an external library) is on the path.
    ok = exist(funcName, 'file') > 0;
    if ok
        fprintf('  [OK] %s  (%s)\n\n', libName, funcName);
    else
        fprintf('  [!!] %s  NOT FOUND  (%s)\n', libName, funcName);
    end
end
