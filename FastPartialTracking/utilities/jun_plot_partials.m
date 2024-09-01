function jun_plot_partials(Partials, time, fs, num_tracks)
% Plots the partials
% 
%
% Julian Neri, 180914
% McGill University, Montreal, Canada

M = length(time);

track_F = -inf+zeros(M,num_tracks);
track_A = -inf+zeros(M,num_tracks);


for m = 1:M    
    active = Partials{m}(:,4);

    track_F(m,active) = Partials{m}(:,1);
    track_A(m,active) = Partials{m}(:,2);
end

plot(time/fs, fs/(2*pi)*track_F, 'linewidth',1.5);
axis tight; grid on; box on;
colororder(["#4c72b0";"#dd8452";"#55a868";"#c44e52";"#8172b3";"#937860";"#da8bc3";"#8c8c8c";"#ccb974";"#64b5cd"]);

end