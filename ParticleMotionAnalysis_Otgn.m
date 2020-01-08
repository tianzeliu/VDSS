%% Particle motion analysis with orthognal error estimation.
%
% History:
% The name was changed from polarization_otgn.m.
% 07/27/2015, Tianze Liu
%
% The code is copied to the current folder and becomes a universal code.
% 03/06/2017, Tianze Liu

function [P,SV,theta_p,theta_s] = ParticleMotionAnalysis_Otgn(R,Z,delta)
    
    % Reshape the inputs to column vectors
    if size(R,2) ~= 1
        R = R';
        Z = Z';
    end
    
    %R = R-mean(R);
    %Z = Z-mean(Z);
    % Define the searching range for P and S
    Theta = 22:delta:89;
    Phi = 112:delta:179;
    n_theta = length(Theta);
    n_phi = length(Phi);
    D = zeros(n_theta,n_phi);
    
    % Grid search
    for i = 1:n_theta
        for j = 1:n_phi
            RZ = [R,Z];
            E1 = [cosd(Theta(i)),cosd(90-Theta(i))];
            E2 = [cosd(Phi(j)),cosd(90-Phi(j))];
            PE1 = RZ*(E1'*E1);
            PE2 = RZ*(E2'*E2);
            DE1 = RZ - PE1;
            DE2 = RZ - PE2;
            D1 = sqrt(DE1(:,1).^2+DE1(:,2).^2);
            D2 = sqrt(DE2(:,1).^2+DE2(:,2).^2);
            Dmin = min(D1,D2);
            D(i,j) = sum(Dmin);
        end
    end
    
    % Find the best-fit directions
    [C1,I1] = min(D);
    [~,I2] = min(C1);
    alpha1 = Theta(I1(I2));
    alpha2 = Phi(I2);
    E1 = [cosd(alpha1),cosd(90-alpha1)];
    E2 = [cosd(alpha2),cosd(90-alpha2)];
    
    %if theta1<theta2
    %    Ep = E2;
    %    Es = E1;
    %    thetaP = theta2;
    %    thetaS = theta1;
    %else
    %    Ep = E1;
    %    Es = E2;
    %    thetaP = theta1;
    %    thetaS = theta2;
    %end
    
    % Project to the best-fit directions
    Ep = E1;
    Es = E2;
    theta_p = alpha1;
    theta_s = alpha2;
    Eps = [Ep',Es'];
    Trans1 = Eps\[1;0];
    Trans2 = Eps\[0;1];
    Trans = [Trans1,Trans2];
    RZ = [R';Z'];
    PSV = Trans*RZ;
    P = PSV(1,:);
    SV = PSV(2,:);

    % Transpose output to column vectors 
    P = P';
    SV = SV';
end