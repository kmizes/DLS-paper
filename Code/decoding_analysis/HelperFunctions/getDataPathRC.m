function [data_path] = getDataPathRC(load_path)

holylfs_samba = '/n/holylfs02/LABS/olveczky_lab/Ashesh/';
holylfs_samba_SW = '/n/holylfs02/LABS/olveczky_lab/Steffen/';


if contains(load_path, 'Dhanashri')
    data_path = [holylfs_samba 'Dhanashri/'];
    
elseif contains(load_path, 'Hamir')
    data_path = [holylfs_samba 'Hamir/'];
    
elseif contains(load_path, 'Hindol')
    data_path = [holylfs_samba 'Hindol/'];
    
elseif contains(load_path, 'Kamod')
    data_path = [holylfs_samba 'Kamod/'];

elseif contains(load_path, 'JaunpuriL')
    data_path = [holylfs_samba 'JaunpuriL/'];
    
elseif contains(load_path, 'Jaunpuri')
    data_path = [holylfs_samba 'Jaunpuri/'];

elseif contains(load_path, 'Gara')
    data_path = [holylfs_samba 'Gara/'];
    
elseif contains(load_path, 'Gandhar')
    data_path = [holylfs_samba 'Gandhar/'];
    
elseif contains(load_path, 'GaudMalhar')
    data_path = [holylfs_samba 'GaudMalhar/'];
    
elseif contains(load_path, 'Gunakari')
    data_path = [holylfs_samba 'Gunakari/'];

elseif contains(load_path, 'Gorakh')
    data_path = [holylfs_samba 'Gorakh/'];
        
elseif contains(load_path, 'SW158')
    data_path = [holylfs_samba 'SW158/'];
    
elseif contains(load_path, 'SW163')
    data_path = [holylfs_samba 'SW163/'];
    
elseif contains(load_path, 'SW116')
    data_path = [holylfs_samba 'SW116/'];
%     data_path = [holylfs_samba_SW 'SW116/'];

elseif contains(load_path, 'SW160')
    data_path = [holylfs_samba 'SW160/'];
%     data_path = [holylfs_samba_SW 'SW160/'];
    
elseif contains(load_path, 'SW166')
    data_path = [holylfs_samba 'SW166/'];
%     data_path = [holylfs_samba_SW 'SW166/'];

elseif contains(load_path, 'SW233')
    data_path = [holylfs_samba 'SW233/'];

end

