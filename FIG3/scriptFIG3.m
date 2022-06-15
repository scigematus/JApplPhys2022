model = createpde(3);
gm = multicuboid(5,5,[5 5 5], 'Zoffset', [0 5 10]);
model.Geometry = gm;
pdegplot(model, 'FaceLabels',"on","FaceAlpha",0.5);
pdegplot(model, 'CellLabels',"on","FaceAlpha",0.5);
%%
sigmaFM = 0.010;
sigmaCU = 0.058;
sigmaPT = 0.093;

lambdaFM = 3.1;
lambdaCU = 500;
lambdaPT = 7.3;

betaFM = 0.37;

gammaFM = (28.0+17.0)/(28.0-17.0);

thetaPT = 0.1*2;

seeFM = 11e-6;
seeCU = 2.0*1.8e-6;
seePT = 2.0*(-5.3e-6);

pelFM = seeFM*300;
pelCU = seeCU*300;
pelPT = seePT*300;

kappaFM = 8.0e-8;
kappaCU = 40.0e-8;
kappaPT = 7.2e-8;

CFM = ccoeffunc(sigmaFM, betaFM, 0, gammaFM, 0, pelFM, seeFM, kappaFM);
CCU = ccoeffunc(sigmaCU, 0, 0, 0, 0, pelCU, seeCU, kappaCU);
CPT = ccoeffunc(sigmaPT, 0, 0, 0, thetaPT, pelPT, seePT, kappaPT);

AFM = [0;-sigmaFM/(lambdaFM^2);0];
ACU = [0;-sigmaCU/(lambdaCU^2);0];
APT = [0;-sigmaPT/(lambdaPT^2);0];
%%
specifyCoefficients(model, 'Cell', 1, 'm', 0, 'd', 0, 'c', CFM, 'a', AFM, 'f', @fcoefFM);
specifyCoefficients(model, 'Cell', 2, 'm', 0, 'd', 0, 'c', CCU, 'a', ACU, 'f', @fcoefCU);
specifyCoefficients(model, 'Cell', 3, 'm', 0, 'd', 0, 'c', CPT, 'a', APT, 'f', @fcoefPT);

H = [0 0 0;0 0 0;0 0 1];
Q = [0 0 0;0 0 0;0 0 0];
G = [0;0;0];

T0 = 300;
Td = 301;
applyBoundaryCondition(model, 'mixed', 'Face', 1, 'h', H, 'r', [0;0;T0],'q', Q, 'g', G);
applyBoundaryCondition(model, 'mixed', 'Face', 12, 'h', H, 'r', [0;0;Td],'q', Q, 'g', G);
generateMesh(model, 'Hmax', 0.4);
pdeplot3D(model);
result = solvepde(model);
u = result.NodalSolution;
%%
pdegplot(model,'CellLabels', 'on', 'FaceAlpha', 0.5);
set(gca, 'FontSize', 15);
saveas(gcf,'thermal_model.png');

pdeplot3D(model, "ColorMapData",u(:,1));
set(gca, 'FontSize', 14);
saveas(gcf,'thermal_u1.png');
pdeplot3D(model, "ColorMapData",u(:,2));
set(gca, 'FontSize', 14);
saveas(gcf,'thermal_u2.png');
pdeplot3D(model, "ColorMapData",u(:,3));
set(gca, 'FontSize', 14);
saveas(gcf,'thermal_u3.png');
%%
xq1 = linspace(0,0,501);
yq1 = linspace(-2.5,2.5,501);
zq1 = linspace(11,11,501);
uintrp = interpolateSolution(result,xq1,yq1,zq1,[1,2,3]);
plot(yq1,uintrp(:,1),'LineWidth',2);
xlabel('Position [nm]');
ylabel('{\it u}_1 (Pt) [V]');
set(gca, 'FontSize', 18);
pbaspect([1 1 1])
saveas(gcf,'thermal_pt.png');
%%
xq2 = linspace(0,0,1501);
yq2 = linspace(0,0,1501);
zq2 = linspace(0,15,1501);
uintrp = interpolateSolution(result,xq2,yq2,zq2,[1,2,3]);

[gradx, grady, gradz] = evaluateGradient(result, xq2, yq2, zq2, [1,2,3]);
grad = [gradx(:,1) grady(:,1) gradz(:,1) gradx(:,2) grady(:,2) gradz(:,2) gradx(:,3) grady(:,3) gradz(:,3)].';
matFM = ccoefmat(sigmaFM, betaFM, 0, gammaFM, 0, pelFM, seeFM, kappaFM);
matCU = ccoefmat(sigmaCU, 0, 0, 0, 0, pelCU, seeCU, kappaCU);
matPT = ccoefmat(sigmaPT, 0, 0, 0, thetaPT, pelPT, seePT, kappaPT);
gradFM = matFM*grad;
gradCU = matCU*grad;
gradPT = matPT*grad;

scatter(zq2(1:500), gradFM(6, 1:500).')
hold on
scatter(zq2(500:1000), gradCU(6, 500:1000).')
scatter(zq2(1000:1500), gradPT(6, 1000:1500).')
hold off
box on
xlabel('Position [nm]');
ylabel('z component of {\it j}_s [A/nm^2]');
legend({'Fe(bottom)','Cu','Pt(top)'})
ylim([-0.5e-8 2.e-8]);
set(gca, 'FontSize', 18);
pbaspect([1 1 1]);
saveas(gcf,'thermal_js.png');


plot(zq2, uintrp(:,1)+uintrp(:,2));
hold on
plot(zq2, uintrp(:,1)-uintrp(:,2));
hold off
plot(zq2, uintrp(:,3));
%%
function cmat = ccoefmat(sigma, beta, alpha, gamma, theta, pel, see, kappa)
    sk1 = theta*(alpha+beta)*0.5;
    sk2 = theta*(1+alpha*beta)*0.5;
    c13 = pel*(1+gamma*beta)*0.5;
    c23 = pel*(gamma+beta)*0.5;
    c31 = -see*(1+gamma*beta)*0.5;
    c32 = see*(gamma+beta)*0.5;
    c33 = -kappa/sigma;
    cmat = sigma*[
        1 0 0 beta 0 0 c31 0 0;
        0 1 -sk1 0 beta -sk2 0 c31 0;
        0 sk1 1 0 sk2 beta 0 0 c31;
        -beta 0 0 -1 0 0 c32 0 0;
        0 -beta sk2 0 -1 sk1 0 c32 0;
        0 -sk2 -beta 0 -sk1 -1 0 0 c32;
        c13 0 0 c23 0 0 c33 0 0;
        0 c13 0 0 c23 0 0 c33 0;
        0 0 c13 0 0 c23 0 0 c33];
end

function c_out = ccoeffunc(sigma, beta, alpha, gamma, theta, pel, see, kappa)
sk1 = theta*(alpha+beta)*0.5;
sk2 = theta*(1+alpha*beta)*0.5;
c13 = pel*(1+gamma*beta)*0.5;
c23 = pel*(gamma+beta)*0.5;
c31 = -see*(1+gamma*beta)*0.5;
c32 = see*(gamma+beta)*0.5;
c33 = -kappa/sigma;
c_out = sigma*...
    [1;0;0;0;1;sk1;0;-sk1;1;
    -beta;0;0;0;-beta;-sk2;0;sk2;-beta;
    c13;0;0;0;c13;0;0;0;c13;
    ...
    beta;0;0;0;beta;sk2;0;-sk2;beta;
    -1;0;0;0;-1;-sk1;0;sk1;-1;
    c23;0;0;0;c23;0;0;0;c23;
    ...
    c31;0;0;0;c31;0;0;0;c31;
    c32;0;0;0;c32;0;0;0;c32;
    c33;0;0;0;c33;0;0;0;c33
    ];
end

function f = fcoefFM(location,state)
sigma = 0.010;
beta = 0.37;
N = 3;
nr = length(location.x);
f = zeros(N,nr);
f(3,:) = -sigma*(state.ux(1,:).*state.ux(1,:) + state.uy(1,:).*state.uy(1,:) + state.uz(1,:).*state.uz(1,:) +...
    state.ux(2,:).*state.ux(2,:) + state.uy(2,:).*state.uy(2,:) + state.uz(2,:).*state.uz(2,:))...
    -2.0*beta*sigma*(state.ux(1,:).*state.ux(2,:) + state.uy(1,:).*state.uy(2,:) + state.uz(1,:).*state.uz(2,:));
end

function f = fcoefCU(location,state)
sigma = 0.058;
beta = 0.0;
N = 3;
nr = length(location.x);
f = zeros(N,nr);
f(3,:) = -sigma*(state.ux(1,:).*state.ux(1,:) + state.uy(1,:).*state.uy(1,:) + state.uz(1,:).*state.uz(1,:) +...
    state.ux(2,:).*state.ux(2,:) + state.uy(2,:).*state.uy(2,:) + state.uz(2,:).*state.uz(2,:))...
    -2.0*beta*sigma*(state.ux(1,:).*state.ux(2,:) + state.uy(1,:).*state.uy(2,:) + state.uz(1,:).*state.uz(2,:));
end

function f = fcoefPT(location,state)
sigma = 0.093;
beta = 0.0;
N = 3;
nr = length(location.x);
f = zeros(N,nr);
f(3,:) = -sigma*(state.ux(1,:).*state.ux(1,:) + state.uy(1,:).*state.uy(1,:) + state.uz(1,:).*state.uz(1,:) +...
    state.ux(2,:).*state.ux(2,:) + state.uy(2,:).*state.uy(2,:) + state.uz(2,:).*state.uz(2,:))...
    -2.0*beta*sigma*(state.ux(1,:).*state.ux(2,:) + state.uy(1,:).*state.uy(2,:) + state.uz(1,:).*state.uz(2,:));
end