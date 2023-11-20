[v, u] = generateOrthogonalVectors(1);


n=100; k=2; nter=50; %nter is the maximum encoding time step t





figure; % Create a new figure
hold on; % Hold the plot
ii=1;
% Determine the color for the start markers
start_u_color = 'r'; % red for w
start_v_color = 'b'; % Blue for w'

% Plot starting points
start_u = plot(u(1), u(2), 'o', 'MarkerSize',8, 'MarkerFaceColor', start_u_color, 'MarkerEdgeColor', start_u_color, 'LineWidth', 2, 'DisplayName', '$w$ Start'); % Plot v1 start marker
start_v = plot(v(1), v(2), 'o', 'MarkerSize',8, 'MarkerFaceColor', start_v_color, 'MarkerEdgeColor', start_v_color, 'LineWidth', 2, 'DisplayName', '$w''$ Start'); % Plot v2 start marker

% Create variables to store old values
old_u = u;
old_v = v;


% entanglement take place here
% we apply encoding directly on u and v for each encoding time step
    
while ii<nter

    % Encoding for u
    [reduced_codeword1,reduced_matrix]=Encoding(u',n,k);
    fu1=reduced_codeword1';
    u=fu1;

    % Encoding for v
    reduced_codeword2=reduced_matrix*(v');
    fu2=reduced_codeword2';
    v=fu2;

    ii=ii+1;
 

    u = u/norm(u);  % Normalize the vector
    v = v/norm(v);  % Normalize the vector

    % Plot updated vectors
    % Draw lines from old points to new points
    line([old_u(1), u(1)], [old_u(2), u(2)], 'Color', 'r', 'LineWidth', 1);
    line([old_v(1), v(1)], [old_v(2), v(2)], 'Color', 'b', 'LineWidth', 1);

    % Store current points as old points for the next iteration
    old_u = u;
    old_v = v;

end



grid on; % Turn the grid on
xlabel('X'); % Label X-axis
ylabel('Y'); % Label Y-axis

% Dummy lines for legend
dummy_u = line([NaN], [NaN], 'Color', 'r', 'LineWidth', 1, 'DisplayName', '$w$ Path');
dummy_v = line([NaN], [NaN], 'Color', 'b', 'LineWidth', 1, 'DisplayName', '$w''$ Path');

% Plot end points
end_u = plot(u(1), u(2), 'ks', 'MarkerSize',10, 'MarkerFaceColor', 'none', 'LineWidth', 2, 'DisplayName', '$w$  End'); % Plot v1 end as black square
end_v = plot(v(1), v(2), 'kd', 'MarkerSize',10, 'MarkerFaceColor', 'none', 'LineWidth', 2, 'DisplayName', '$w''$ End'); % Plot v2 end as black diamond

legend([start_u, start_v, dummy_u, dummy_v, end_u, end_v],'Interpreter', 'latex'); % This will automatically grab the 'DisplayName' properties of your plot objects
hold off; % Release the plot








function [yfil,frmat]=Encoding(x,n,t)
k=length(x);

rmat=randn(n,k);
% rmat=orth(rmat);
y=rmat*x;
absy=abs(y);
[sorted_data, sortedindex ]= sort(absy, 'descend');
topindex=(sortedindex(1:t));
frmat=rmat(topindex,:);

yfil=y(topindex);

end






function [targetVec, orthogonalVecs] = generateOrthogonalVectors(m)
    targetVec = randn(1, 2);  % Generate target vector
    targetVec = targetVec / norm(targetVec);  % Normalize the vector

    orthogonalVecs = zeros(m, 2);  % Preallocate matrix for efficiency

    for i = 1:m
        % Generate a random angle in radians
        angle = rand() * 2 * pi;  % Angle in [0, 2*pi]
        % Create a rotation matrix
        R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
        % Rotate the target vector
        rotatedVec = (R * targetVec')';
        
        % Ensure the generated vector is orthogonal to the previous ones
        if i > 1
            for j = 1:i-1
                % Project the rotated vector onto the previously generated vectors
                projection = dot(rotatedVec, orthogonalVecs(j, :)) / dot(orthogonalVecs(j, :), orthogonalVecs(j, :));
                % Subtract the projection from the rotated vector
                rotatedVec = rotatedVec - projection * orthogonalVecs(j, :);
            end
        end
        
        % Normalize the orthogonal vector
        rotatedVec = rotatedVec / norm(rotatedVec);
        
        % Store the orthogonal vector in the matrix
        orthogonalVecs(i, :) = rotatedVec;
    end

end
