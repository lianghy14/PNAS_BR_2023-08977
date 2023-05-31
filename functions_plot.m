function pfuns = functions_plot
pfuns.potential = @plot_potential;
end


function number = plot_potential(potential,x,y,boundary_H)
global h
gfuns = functions_given;
nx = length(x);
ny = length(y);
boundary_judge = zeros(ny,nx);
boundary_judge = gfuns.Boundary_value(x,y,boundary_judge,boundary_H,1);
[X,Y] = meshgrid(x,y);
potential(potential>=10^5) = nan;
potential =  potential.*(1-boundary_judge);
plotrange=[min(x),max(x),min(y),max(y),0,max(max(potential))];

figure;

subplot(2,1,1);
potential(potential<=0) = nan;
set(surf(X,Y,potential),'Facecolor','none'); axis(plotrange);
title(['dx = dy = ',num2str(h)])
xlabel('x(m)'); ylabel('y(m)');
subplot(2,1,2); 
imagesc(x,y, potential,'alphadata',~isnan(potential)); colorbar;
% title(['Method: ',method ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t)])
xlabel('x(m)'); ylabel('y(m)'); 

number = 0;
end