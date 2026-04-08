function DGRF2020 = IGRF2020MagneticField
    
    % --- Data IGRF2020
    n_values = [1, 1, 2, 2, 2, 3, 3, 3]; 
    m_values = [0, 1, 0, 1, 2, 0, 1, 3]; 
    g_nm = [-29404.8, -1450.9, -2499.6, 2982.0, 1677.0, 1363.2, -2381.2, 525.7]; % [nT]
    h_nm = [0, 4652.5, 0, -2991.6, -734.6, 0, -82.1, -543.4]; % [nT]
    
    g_nm = g_nm * 10^-9; 
    h_nm = h_nm * 10^-9; 
    
    DGRF2020 = struct();
    
    for i = 1:length(n_values)
        DGRF2020(i).n = n_values(i);              
        DGRF2020(i).m = m_values(i);              
        DGRF2020(i).g_nm = g_nm(i);               
        DGRF2020(i).h_nm = h_nm(i);               
        DGRF2020(i).description = sprintf('Harmonic data for n=%d, m=%d', ...
                                           n_values(i), m_values(i)); 
    end
end 