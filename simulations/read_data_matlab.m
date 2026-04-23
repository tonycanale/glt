out_dir = 'path/to/your/csv/files';  % adjust this

m_dim   = [20, 50, 100];
T_dim   = [20, 40; 50, 100; 100, 200];
H_dim   = [3, 8, 10;
            5, 10, 15;
            5, 10, 15;
            8, 20, 40];
n_reps  = 20; 

scenario_names = {'scenario1', 'scenario2', 'scenario3', 'scenario4'};

data = struct();

for s = 1:4
    scen = scenario_names{s};
    for dimens = 1:3
        m = m_dim(dimens);
        H = H_dim(s, dimens);
        for samplesize = 1:2
            T = T_dim(dimens, samplesize);
            for rep = 1:n_reps
                key = sprintf('%s_m%d_H%d_T%d_rep%02d', scen, m, H, T, rep);

                f_Y      = fullfile(out_dir, sprintf('%s_Y_%s.csv',      scen, sprintf('m%d_H%d_T%d_rep%02d', m, H, T, rep)));
                f_Delta  = fullfile(out_dir, sprintf('%s_Delta_%s.csv',  scen, sprintf('m%d_H%d_T%d_rep%02d', m, H, T, rep)));
                f_Lambda = fullfile(out_dir, sprintf('%s_Lambda_%s.csv', scen, sprintf('m%d_H%d_T%d_rep%02d', m, H, T, rep)));

                data.(key).Y      = readmatrix(f_Y);
                data.(key).Delta  = readmatrix(f_Delta);
                data.(key).Lambda = readmatrix(f_Lambda);
            end
        end
    end
end