count =[] ; 
for i =1:150
    count= [count sum(abs(result_hbic(i).tstat)>1.96)];
end

rate = count ./201;
figure;
bar(1:150, rate,  'FaceColor', [0.5, 0.5, 0.5]); % Line plot with black color
xlabel('Factor ID');
ylabel('significant Rate');
% title('significant rate  of 201 test');


% Set grid transparency to make it lighter
% set(gca, 'GridAlpha', 0.3); % Make grid lines lighter
% 
% % Optionally adjust the appearance
% set(gca, 'FontSize', 12);
