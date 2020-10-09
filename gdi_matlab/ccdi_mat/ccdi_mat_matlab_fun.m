function [ est_CDI, est_CDI_list ] = ccdi_mat_matlab_fun( X, Y, Z, B )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes Classifier DI conditioned on Z
%
% B is the number of bootstrap iterations to use.
%
% Calls the python code for CCMI and uses file saving and loading to pass X
% Y, & Z, which should have dimensions as columns and rows as observations
%
% Copyright (C) 2020 Joseph Young - see GPLv2_note.txt for full notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Save file for python to access, with labindex included so that when
    % running this function parallel the same file will not be accessed by
    % multiple processes
    t = getCurrentTask(); 
    if isempty(t)
        t.ID = 1;
    end
    save([pwd '/' sprintf('ccdi_mat/data_for_ccdi_%i.mat', t.ID)],...
         'X','Y','Z','-v7')
    
    % Load python path - script should have variable named 'python_path'
    % which contains a string for the python path
    python_path_script
    
    % Call Python script to run on saved X,Y,&Z for MI to be estimated by
    % CCMI
    system(sprintf('export PATH=%s ;%s ccdi_mat/ccdi_mat_py_fun.py %i %i',...
      system_path,python_path,t.ID,B));
        
    % Load CDI estimate computed by python code for CCMI, which will have
    % name est_CDI. est_CDI_list will also be loaded, containing the
    % results of each iteration (est_CDI is the average across these
    % iterations which are bootstrapped).
    load([pwd '/' sprintf('ccdi_mat/ccdi_est_CDI_%i.mat', t.ID)])
    
    % Remove intermediate files
    delete([pwd '/' sprintf('ccdi_mat/data_for_ccdi_%i.mat', t.ID)])
    delete([pwd '/' sprintf('ccdi_mat/ccdi_est_CDI_%i.mat', t.ID)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
