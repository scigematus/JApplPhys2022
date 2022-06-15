tNM = 5;
tFM = 5;
width = 5;
model = createpde(2);
gm = multicuboid(width,width,[tNM tFM tNM],"Zoffset",[0 tNM tNM+tFM]);
model.Geometry = gm;
pdegplot(model,'Facelabels', 'on', 'FaceAlpha', 0.5);
pdegplot(model,'Edgelabels', 'on', 'FaceAlpha', 0.5);
pdegplot(model,'CellLabels', 'on', 'FaceAlpha', 0.5);
saveas(gcf,'nmfmnm_model.png');
%%
beta = 0.37;
sigmaNM = 9.3e-2;
sigmaFM = 1.0e-2;
lambdaNM = 7.3;
lambdaFM = 3.1;
CNM = ccoeffunc(sigmaNM, 0);
CFM = ccoeffunc(sigmaFM, beta);
ANM = [0;-sigmaNM/(lambdaNM^2)];
AFM = [0;-sigmaFM/(lambdaFM^2)];
F = [0;0];
%%
specifyCoefficients(model,"Cell",1,"m",0,"d",0,"c",CNM,"a",ANM,"f",F);
specifyCoefficients(model,"Cell",2,"m",0,"d",0,"c",CFM,"a",AFM,"f",F);
specifyCoefficients(model,"Cell",3,"m",0,"d",0,"c",CNM,"a",ANM,"f",F);
applyBoundaryCondition(model,'dirichlet','Face',1,'h',[1 0;0 0],'r',[0;0]);
applyBoundaryCondition(model,'dirichlet','Face',12,'h',[1 0;0 0],'r',[1;0]);
generateMesh(model,"Hmax",0.3);
pdeplot3D(model);
result = solvepde(model);
u = result.NodalSolution;
%%
pdegplot(model,'CellLabels', 'on', 'FaceAlpha', 0.5);
set(gca, 'FontSize', 15);
saveas(gcf,'nmfmnm_model.png');

pdeplot3D(model, "ColorMapData",u(:,1));
set(gca, 'FontSize', 14);
saveas(gcf,'nmfmnm_u13d.png');
pdeplot3D(model, "ColorMapData",u(:,2));
set(gca, 'FontSize', 14);
saveas(gcf,'nmfmnm_u23d.png');
%%
fine = 501;
xq = linspace(0,0,fine);
yq = linspace(0,0,fine);
zq = linspace(0,tNM+tFM+tNM,fine);
uintrp = interpolateSolution(result,xq,yq,zq,[1,2]);

plot(zq,uintrp(:,1),'LineWidth',2);
xlabel('Position [nm]');
ylabel('{\it u}_1 [V]');
set(gca, 'FontSize', 20);
pbaspect([1.2 1 1])
saveas(gcf,'nmfmnm_u1.png');

plot(zq,uintrp(:,2),'LineWidth',2);
xlabel('Position [nm]');
ylabel('{\it u}_2 [V]');
set(gca, 'FontSize', 20);
pbaspect([1.2 1 1])
saveas(gcf,'nmfmnm_u2.png');
%%
[gradx,grady,gradz] = evaluateGradient(result,xq,yq,zq,[1,2]);
grad = [gradx(:,1) grady(:,1) gradz(:,1) gradx(:,2) grady(:,2) gradz(:,2)].';
matNM = ccoefmat(sigmaNM, 0);
matFM = ccoefmat(sigmaFM, beta);
gradcnm = matNM*grad;
gradcfm = matFM*grad;

faceA = round(fine*tNM/(tNM+tFM+tNM));
faceB = round(fine*(tNM+tFM)/(tNM+tFM+tNM));
scatter(zq(1:faceA), gradcnm(6, 1:faceA).')
hold on
scatter(zq(faceA:faceB), gradcfm(6, faceA:faceB).')
scatter(zq(faceB:fine), gradcnm(6, faceB:fine).')
hold off
box on
xlabel('Position [nm]');
ylabel('z component of {\it j}_s [A/nm^2]');
legend({'NM(bottom)','FM','NM(top)'})
ylim([-0.7e-3 0.3e-3])
set(gca, 'FontSize', 20);
pbaspect([1.2 1 1])
saveas(gcf,'nmfmnm_js.png');
%%
function ccoef = ccoeffunc(sigma, beta)
    ccoef = sigma*[1;0;0;0;1;0;0;0;1;
    -beta;0;0;0;-beta;0;0;0;-beta;
    beta;0;0;0;beta;0;0;0;beta;
    -1;0;0;0;-1;0;0;0;-1];
end

function cmat = ccoefmat(sigma, beta)
    cmat = sigma*[
        1 0 0 beta 0 0;
        0 1 0 0 beta 0;
        0 0 1 0 0 beta;
        -beta 0 0 -1 0 0;
        0 -beta 0 0 -1 0;
        0 0 -beta 0 0 -1];
end