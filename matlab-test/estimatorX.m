close all
clear all
clc


currPhi = -1;

getSingleD    = @(phi, r) (sin(phi)*r(1) - cos(phi)*r(2))^2;
getSingleDPhi = @(phi, r) -2*(r(1)*cos(phi) + r(2)*sin(phi))*(r(2)*cos(phi) - r(1)*sin(phi));

M = 500;

numPoints = 20;


alpha = 1e-4;

residual = zeros(1, M);

for k = 1 : 1 : M
    mu = [10; 20] + [0.001; -0.002]*100*k;
    
    R = diag([20, 20])*1e-1;
    yMat = repmat(mu, 1, numPoints) + R*randn(numPoints, 2)';
    
    t = linspace(-50, 50, 2);
    lineX = cos(currPhi)*t;
    lineY = sin(currPhi)*t;
    
    plot(yMat(1, :), yMat(2, :), 'o');
    axis equal
    grid on
    xlim([-30, 30])
    ylim([-30, 30])
    hold on
    plot(lineX, lineY);
    hold off
    
    pause(0.01);    
    
    grad = 0;
    dist = 0;
    
    for i = 1 : 1 : numPoints
        grad = grad + getSingleDPhi(currPhi, yMat(:, i));
        dist = dist + getSingleD(currPhi, yMat(:, i));
    end
    grad = grad/numPoints; 
    dist = dist/numPoints;
        
    
    currPhi = currPhi - alpha*grad;
    disp(currPhi)
    residual(k) = dist;
end






figure
plot(residual)