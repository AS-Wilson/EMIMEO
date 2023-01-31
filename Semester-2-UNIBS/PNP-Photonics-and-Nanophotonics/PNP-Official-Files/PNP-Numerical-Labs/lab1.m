

%Scenario 1
N_angles=90;
angles=linspace(0,89,N_angles);
for ii=1:N_angles
    T_R_A=tmm(600e-9,angles(ii),1,0,1.33,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end
figure; plot(angles, T,'LineWidth',3);
xlabel('Angle of Incidence (deg)', 'FontName', 'Times New Roman', 'FontSize', 20);
ylabel('T, R, A')

brewster(1)=atan(1/1.33)*180/pi;
critical(1)=asin(1/1.33)*180/pi;

% Scenario 2
N_angles=90;
angles=linspace(0,89,N_angles);
for ii=1:N_angles
    T_R_A=tmm(600e-9,angles(ii),1,0,1.33,1,1);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end
% figure; plot(angles, T,'LineWidth',3);
% xlabel('Angle of Incidence (deg)', 'FontName', 'Times New Roman', 'FontSize', 20);
% ylabel('T, R, A')

brewster(2)=atan(1/1.33)*180/pi;
critical(2)=asin(1/1.33)*180/pi;


%Scenario 3
N_angles=90;
angles=linspace(0,89,N_angles);
for ii=1:N_angles
    T_R_A=tmm(600e-9,angles(ii),1,0,1,1.33,1);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end
figure; plot(angles, T,'LineWidth',3);
xlabel('Angle of Incidence (deg)', 'FontName', 'Times New Roman', 'FontSize', 20);
ylabel('T, R, A')

brewster(3)=atan(1.33/1)*180/pi;
critical(3)=asin(1.33/1)*180/pi;