N = 10;
fine = 50;
theta = linspace(-90,90,N);
position = double.empty(0,2*fine+1);
potential = double.empty(2*fine+1,0);
voltage = double.empty(0,1);
%%
lx = 12;
ly = 20;
lz = 4;
j0 = -1e-9;
s = 0.0026;%sigma
b = 0.71;%beta
amr = 1/1.043;
%
s_tt = s*0.005;%sigma_she^tt
s_lt = s*0.0075;
s_tl = s*0.0075;
s_itt = s*0.005;%sigma_ishe^tt
s_ilt = s*0.0075;
s_itl = s*0.0075;
%
s_a = s_tt*b;
s_ia = s_itt*b;
%
lam_l = 2.8;
lam_t = 1.0;
%pdegplot(model,'faceLabels','on','FaceAlpha',0.5)
C = [amr*s;0;0;0;s;s_a;0;-s_a;s;%
     -b*amr*s;0;0;0;-b*s;-s_tt;0;s_tt;-b*s;
     0;0;s_lt;0;0;0;-s_tl;0;0;
     0;-s_lt;0;s_tl;0;0;0;0;0;
     b*amr*s;0;0;0;b*s;s_itt;0;-s_itt;b*s;%
     -amr*s;0;0;0;-s;-s_ia;0;s_ia;-s;
     0;0;0;0;0;0;0;0;0;
     0;0;0;0;0;0;0;0;0;
     0;0;-s_ilt;0;0;0;s_itl;0;0;%
     0;0;0;0;0;0;0;0;0;
     -amr*s;0;0;0;-s;0;0;0;-s;
     0;0;0;0;0;0;0;0;0;
     0;s_ilt;0;-s_itl;0;0;0;0;0;%
     0;0;0;0;0;0;0;0;0;
     0;0;0;0;0;0;0;0;0;
     -amr*s;0;0;0;-s;0;0;0;-s];
A = [0;-s*lam_l^-2;-s*lam_t^-2;-s*lam_t^-2];
F = [0;0;0;0];
G0 = [0;0;0;0];
%% 
% In case charge current flow into Face 6.

for i = 1:N
    gm = multicuboid(lx,ly,lz);
    rotate(gm,-theta(i),[0 0 lz/2],[0 1 lz/2]);
    model = createpde(4);
    model.Geometry = gm;
    pdegplot(model,'CellLabels', 'on', 'FaceAlpha', 0.5);
    set(gca, 'FontSize', 14);
    saveas(gcf,strcat('model',num2str(theta(i)),'.png'));
    applyBoundaryCondition(model,'neumann','Face',6,'g',[j0;0;0;0]);
    applyBoundaryCondition(model,'neumann','Face',4,'g',[-j0;0;0;0]);
    applyBoundaryCondition(model,'neumann','Face',[1,2,3,5],'g',G0);
    specifyCoefficients(model,'m',0,'d',0,'c',C,'a',A,'f',F);
    generateMesh(model,'Hmax',0.8);
    result = solvepde(model);
    u = result.NodalSolution;
    pdeplot3D(model,'ColorMapData',u(:,1));
    set(gca, 'FontSize', 14);
    saveas(gcf,strcat('u1',num2str(theta(i)),'.png'));
    pdeplot3D(model,'ColorMapData',u(:,2));
    set(gca, 'FontSize', 14);
    saveas(gcf,strcat('u2',num2str(theta(i)),'.png'));
    pdeplot3D(model,'ColorMapData',u(:,3));
    set(gca, 'FontSize', 14);
    saveas(gcf,strcat('u3',num2str(theta(i)),'.png'));
    pdeplot3D(model,'ColorMapData',u(:,4));
    set(gca, 'FontSize', 14);
    saveas(gcf,strcat('u4',num2str(theta(i)),'.png'));
    xq = linspace(-0.49*lx*cosd(theta(i)),0.49*lx*cosd(theta(i)),2*fine+1);
    yq = linspace(0,0,2*fine+1);
    zq = linspace(0.5*lz-0.49*lx*sind(theta(i)),0.5*lz+0.49*lx*sind(theta(i)),2*fine+1);
    uintrp = interpolateSolution(result,xq,yq,zq,[1,2,3,4]);
    position = cat(1, position, cat(2,-sqrt((zq(1:fine)-0.5*lz).^2 + xq(1:fine).^2), sqrt((zq(fine+1:2*fine+1)-0.5*lz).^2 + xq(fine+1:2*fine+1).^2)));
    potential = cat(2, potential, uintrp(:,1)+b*uintrp(:,2)-(uintrp(fine+1,1)+b*uintrp(fine+1,2)));
    voltage = cat(1, voltage, uintrp(2*fine+1,1)+b*uintrp(2*fine+1,2)-uintrp(1,1)-b*uintrp(1,2));
end
%%
for i = 1:N
    plot(position(i,:), potential(:,i),'LineWidth',2);
    hold on
end
hold off
box on
xlabel('Position');
ylabel('Voltage [V] (change from position zero)');
f = @(x) (strcat(num2str(x),'deg'));
legend(arrayfun(f, theta, 'UniformOutput', false), 'Location','northeastoutside');
ylim([-8e-9 8e-9]);
set(gca, 'FontSize', 16);
pbaspect([1 1 1]);
saveas(gcf,'ahe_profile.png');

scatter(theta,voltage,80,'filled');
box on
xlabel('Angle [deg]');
ylabel('Voltage difference [V]');
ylim([-2e-8 2e-8]);
xlim([-90 90]);
xticks(linspace(-90,90,7));
set(gca, 'FontSize', 18);
pbaspect([1 1 1]);
saveas(gcf,'ahe_angle.png');

writematrix(position);
writematrix(potential);
writematrix(theta);
writematrix(voltage);
%% 
%