function [map, sig, cohname] = getCohmap(cohtype, mapfname, sigfname)
%PARSECOHTYPE Summary of this function goes here
%   Detailed explanation goes here
    load(mapfname);
    
    if (exist(sigfname, 'file') == 2)
        load(sigfname);
        switch cohtype
            case 'CMC'
                map = CMC_Map;
                sig = CMC_Sig;
                cohname = 'CMC';
            case 'gPDC1'
                map = PDC1_Map;
                sig = PDC1_Sig;
                cohname = 'gPDC  (EEG -> EMG)';
            case 'gPDC2'
                map = PDC2_Map;
                sig = PDC2_Sig;
                cohname = 'gPDC  (EMG -> EEG)';
            case 'PDC1'
                map = gPDC1_Map;
                sig = gPDC1_Sig;
                cohname = 'PDC  (EEG -> EMG)';
            case 'PDC2'
                map = gPDC2_Map;
                sig = gPDC2_Sig;
                cohname = 'PDC  (EMG -> EEG)';
            case 'iCOH1' 
                map = iCOH1_Map;
                sig = iCOH1_Sig;
                cohname = 'iCOH  (EEG -> EMG)';
            case 'iCOH2' 
                map = iCOH2_Map;
                sig = iCOH2_Sig;
                cohname = 'iCOH  (EMG -> EEG)';
        end    
    else
        switch cohtype
            case 'CMC'
                map = CMC_Map;
                cohname = 'CMC';
            case 'gPDC1'
                map = PDC1_Map;
                cohname = 'gPDC  (EEG -> EMG)';
            case 'gPDC2'
                map = PDC2_Map;
                cohname = 'gPDC  (EMG -> EEG)';
            case 'PDC1'
                map = gPDC1_Map;
                cohname = 'PDC  (EEG -> EMG)';
            case 'PDC2'
                map = gPDC2_Map;
                cohname = 'PDC  (EMG -> EEG)';
            case 'iCOH1' 
                map = iCOH1_Map;
                cohname = 'iCOH  (EEG -> EMG)';
            case 'iCOH2' 
                map = iCOH2_Map;
                cohname = 'iCOH  (EMG -> EEG)';
        end   
        [x,y] = size(map);
        sig = zeros(x,y);
    % Get rid of first components which are sometimes zero
    switch cohtype
        case 'CMC'
            map = map;
            sig = sig;
        otherwise
            map = map(2:end, :);
            sig = sig(2:end, :);
    end
end

