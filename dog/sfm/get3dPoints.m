function threeDPoints = get3dPoints(all_2d_points)


kp_names = {'REye', 'LEye', 'nosetip', 'RHead', 'MHead', 'LHead', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13'};
kp_perm = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

n_dogs = length(all_2d_points);

kps_all = [];
vis_all = [];
SCALE = 224;

lr_edges = [2 1; 6 4]; %left to right edges (along -X)
bf_edges = [5 3]; % back to front edges (along -Y)



box_trans = zeros(n_dogs, 2);
for b = 1:n_dogs

    kps_b = squeeze(all_2d_points(b, 1:2, :))/SCALE;
    
    % Add flipped data
    kps_b_flipped = kps_b(:, kp_perm);
    kps_b_flipped(1, :) = -kps_b_flipped(1, :);
    
    vis_b = squeeze(vertcat(all_2d_points(b, 3, :),all_2d_points(b, 3, :)));
    vis_b_flipped = vis_b(:, kp_perm);
    
    % Mean center here,,
    box_trans(b, :) = mean(kps_b(:, vis_b(1,:)>0), 2);
    kps_b = kps_b - box_trans(b, :)';
    kps_b_flipped = kps_b_flipped - mean(kps_b_flipped(:, vis_b_flipped(1,:)>0), 2);
    
    kps_all = vertcat(kps_all, kps_b, kps_b_flipped);
    
    vis_all = vertcat(vis_all, vis_b, vis_b_flipped);
    % keyboard
    % sfigure(2); clf;
    % scatter(kps_b(1,:), kps_b(2,:));
    % hold on;
    % scatter(kps_b_flipped(1,:), kps_b_flipped(2,:));
end

[~, S, ~] = sfmFactorization(kps_all, 30, 10);

    S = diag([1 1 -1])*S;
    
    R = alignSfmModel(S, lr_edges, bf_edges, []);
   
    S = R*S;

    if (S(3, 3) < S(3, 13))
        S = diag([1 1 -1])*S;
        R = alignSfmModel(S, lr_edges, bf_edges, []);
   
        S = R*S;
    end



    max_dist = max(pdist(S'));
    S_scale = 2. / max_dist;
    S = S*S_scale;
    % [M,T,~] = sfmFactorizationKnownShape(kps_all, S, 50);

    threeDPoints = S;
end