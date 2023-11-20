
clear all
n=10; k=2; nter=20; %nter is the encoding time step t

for jj=1:50  % generate 50 random pair
   [v, u] = generateOrthogonalVectors(1); % generate vector u and v othogonal to each other
   

    ii=1;
  
    % Create variables to store old values
    old_u = u;
    old_v = v;



    % entanglement take place here
    % we apply encoding directly on u and v for each encoding time step
    pointr=[];
    while ii < nter

        % Encoding for u
        [reduced_codeword,reduced_matrix]=Encoding(u',n,k);
        fu1=reduced_codeword';
        u=fu1;

        % Projection for v
        fu2=reduced_matrix*(v');
        fu2=fu2';
        v=fu2;

        ii=ii+1;


        u = u/norm(u);  % Normalize the vector
        v = v/norm(v);  % Normalize the vector
        
        
        dis=dot(u, v);% compute their dot product
        
        
        pointr=[pointr;dis];

      
        % Store current points as old points for the next iteration
        old_u = u;
        old_v = v;
    end

    pointrR{jj}=pointr;


end
% Create a figure
figure;

% Plot each set of values
for i = 1:numel(pointrR)
    plot(1:length(pointrR{i}), pointrR{i}, '-o', 'DisplayName', ['Time ' num2str(i)],'MarkerSize',15, 'MarkerIndices',1:3:nter, 'LineWidth', 4);
    hold on;
end

xlabel('Encoding Time Step $(t)$', 'Interpreter', 'latex');
ylabel('$\cos \theta_t\pi = (w_t \cdot w''_t)$', 'Interpreter', 'latex');


ylim([-1, 1]); xlim([1, nter]);
% Hold off to stop superimposing new plots
hold off;
% 






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




