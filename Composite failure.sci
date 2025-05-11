clc
clear

    Nx = input("Insira o valor de Nx (N/mm): ")*10^3;
    Ny = input("Insira o valor de Ny (N/mm): ")*10^3;
    Nxy = input("Insira o valor de Nxy (N/mm): ")*10^3;
    Mx = input("Insira o valor de Mx (Momento em Nm/mm): ")*10^3;
    My = input("Insira o valor de My (Momento em Nm/mm): ")*10^3;
    Mxy = input("Insira o valor de Mxy (Momento de cisalhamento em Nm/mm): ")*10^3;
    numCamadas = input("Insira o número de camadas: ");
    espessura = input("Insira a espessura de cada camada (mm): ");
    angulos = zeros(1, numCamadas);
    for i = 1:numCamadas
        angulos(i) = input("Insira o ângulo da camada " + string(i) + " (graus): ");
    end
    E1 = input("Insira o valor de E1 (GPa): ") * 10^9;
    E2 = input("Insira o valor de E2 (GPa): ") * 10^9;
    v12 = input("Insira o valor de v12: ");
    G12 = input("Insira o valor de G12 (GPa): ") * 10^9;
    XT = input("Insira o valor de XT (MPa): ") * 10^6;
    XC = input("Insira o valor de XC (MPa): ") * 10^6;
    YT = input("Insira o valor de YT (MPa): ") * 10^6;
    YC = input("Insira o valor de YC (MPa): ") * 10^6;
    S12 = input("Insira o valor de S12 (MPa): ") * 10^6;

    v21 = (v12 * E2) / E1;
    Q11 = E1 / (1 - v12 * v21);
    Q22 = E2 / (1 - v12 * v21);
    Q12 = (v12 * E2) / (1 - v12 * v21);
    Q66 = G12;
    Qlocal = [Q11, Q12, 0; Q12, Q22, 0; 0, 0, Q66];


for k = 1:numCamadas
    theta = angulos(k) * %pi / 180;
    m = cos(theta);
    n = sin(theta);
    Q11 = Qlocal(1,1);
    Q12 = Qlocal(1,2);
    Q22 = Qlocal(2,2);
    Q66 = Qlocal(3,3);
    Q11k = Q11*m^4 + 2*(Q12 + 2*Q66)*m^2*n^2 + Q22*n^4;
    Q22k = Q11*n^4 + 2*(Q12 + 2*Q66)*m^2*n^2 + Q22*m^4;
    Q12k = (Q11 + Q22 - 4*Q66)*m^2*n^2 + Q12*(m^4 + n^4);
    Q16k = (Q11 - Q12) * n * m^3 + (Q12-Q22)*n^3 * m - 2*m*n*(m^2 + n^2) * Q66;
    Q26k = (Q11 - Q12) * n^3 * m + (Q12-Q22)*n * m^3 + 2*m*n*(m^2 + n^2) * Q66;
    Q66k = (Q11 + Q22 -2*Q12 - 2*Q66)*n^2 * m^2 + Q66 * (n^4 + m^4);
    Qk = [Q11k, Q12k, Q16k; Q12k, Q22k, Q26k; Q16k, Q26k, Q66k];
    Qks(:,:,k) = Qk;
end

zn=-(espessura*numCamadas)/2

    A = zeros(3, 3);
    B = zeros(3, 3);
    D = zeros(3, 3);
    for k = 1:numCamadas
        h_k = zn + espessura;
        h_k_minus_1 = zn;
        A = A +  Qks(:,:,k) * (h_k - h_k_minus_1);
        B = B +  Qks(:,:,k) * 0.5 * (h_k^2 - h_k_minus_1^2);
        D = D + Qks(:,:,k) * (1/3) * (h_k^3 - h_k_minus_1^3);
        zn = h_k
    end

R = [1 0 0; 0 1 0; 0 0 2];

r = inv(R);

    ABD_inv = inv([A, B; B, D]);
    N = [Nx; Ny; Nxy];
    M = [Mx; My; Mxy];
    deformacoes_momento = ABD_inv * [N; M];
    epsilon0 = deformacoes_momento(1:3);
    kappa = deformacoes_momento(4:6);
    
zb = zeros(1, numCamadas);
zt = zeros(1, numCamadas);
zm = zeros(1, numCamadas);

zb(1) = -espessura * numCamadas / 2; // Posição inicial (borda inferior da primeira lâmina)
for k = 1:numCamadas
    if k > 1 then
        zb(k) = zt(k-1); // Borda inferior da lâmina atual é a borda superior da lâmina anterior
    end
    zt(k) = zb(k) + espessura; // Borda superior da lâmina atual
    zm(k) = (zb(k) + zt(k)) / 2; // Posição média da lâmina atual
end

for k = 1:numCamadas
    z = zm(k); 
    sigma_global = Qks(:,:,k) * (epsilon0 - z * kappa);
end

tensoes_locais = cell(1, numCamadas);

for k = 1:numCamadas
    m = cos(theta);
    n = sin(theta);
    T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, (m^2 - n^2)];
    T_matrices(:,:,k) = T;
    sigma_local = T * sigma_global

    sigma1 = sigma_local(1);
    sigma2 = sigma_local(2);
    tau12 = sigma_local(3);
    
tensoes_locais{k} = sigma_local;
    
    disp("Tensões locais na lâmina " + string(k) + ":");
    disp("Sigma X: " + string(sigma_local(1)));
    disp("Sigma Y: " + string(sigma_local(2)));
    disp("Tau XY: " + string(sigma_local(3)));
    disp(" ");
    
end

for i = 1:numCamadas
    defgzm(:,:,i) = epsilon0 + (zm(i) * kappa);
end

for i = 1:numCamadas
    defgzt(:,:,i) = epsilon0 + (zt(i) * kappa);
end    

for i = 1:numCamadas
    defgzb(:,:,i) = epsilon0 + (zb(i) * kappa);
end

for i = 1:numCamadas
    defloczt(:,:,i) = R * T_matrices(:,:,i) * r * defgzt(:,:,i);
    defloczm(:,:,i) = R * T_matrices(:,:,i) * r * defgzm(:,:,i);
    defloczb(:,:,i) = R * T_matrices(:,:,i) * r * defgzb(:,:,i);
end 

    Q_global(:,:,i) = inv(T_matrices(:,:,i)) * Qlocal * R * T_matrices(:,:,i) * inv(R);

for i = 1:numCamadas
    epsilon_medio = deformacoes_momento(1:3); 
    tglobalzt(:,:,i) = epsilon_medio .* defgzt(:,:,i); 
    tglobalzm(:,:,i) = epsilon_medio .* defgzm(:,:,i); 
    tglobalzb(:,:,i) = epsilon_medio .* defgzb(:,:,i); 
end

for i = 1:numCamadas
    tlocalzt(:,:,i) = R * T_matrices(:,:,i) * r * tglobalzt(:,:,i);
    tlocalzm(:,:,i) = R * T_matrices(:,:,i) * r * tglobalzm(:,:,i);
    tlocalzb(:,:,i) = R * T_matrices(:,:,i) * r * tglobalzb(:,:,i);
end
    
z = linspace(0, numCamadas * espessura, numCamadas + 1);
    
tensoes_sigma1 = zeros(1, length(z));
tensoes_sigma2 = zeros(1, length(z));
tensoes_tau12 = zeros(1, length(z));

for k = 1:numCamadas
    margemSegurancaTensao = 1;
    if sigma1 >= XT then
        margemSegurancaTensao = (XT / sigma1) - 1;
    elseif sigma1 <= -XC then
        margemSegurancaTensao = (-XC / sigma1) - 1;
    elseif sigma2 >= YT then
        margemSegurancaTensao = (YT / sigma2) - 1;
    elseif sigma2 <= -YC then
        margemSegurancaTensao = (-YC / sigma2) - 1;
    elseif abs(tau12) >= S12 then
        margemSegurancaTensao = (S12 / abs(tau12)) - 1;
    end
    
    
    disp("Camada " + string(k));
    disp("Sigma1: " + string(sigma1) + ", Sigma2: " + string(sigma2) + ", Tau12: " + string(tau12));
    disp("Margem de Segurança (Máxima Tensão): " + string(margemSegurancaTensao));
    
    if margemSegurancaTensao < 0 then
        disp("Resultado: A camada " + string(k) + " falhou.");
    else
        disp("Resultado: A camada " + string(k) + " está segura.");
    end
end
    
for k = 1:length(z)
    z_k = z(k); 
    sigma_global = Qk * (epsilon0 + z_k * kappa); 
    sigma_local = T_matrices(:,:,i) * sigma_global; 

    tensoes_sigma1(k) = sigma_local(1); 
    tensoes_sigma2(k) = sigma_local(2); 
    tensoes_tau12(k) = sigma_local(3); 
end

scf(1); 
plot(z, tensoes_sigma1, 'r', 'LineWidth', 2); 
plot(z, tensoes_sigma2, 'g', 'LineWidth', 2); 
plot(z, tensoes_tau12, 'b', 'LineWidth', 2);  
title('Tensões ao Longo da Espessura');
xlabel('Posição ao longo da espessura (mm)');
ylabel('Tensão (Pa)');
legend(['\\sigma_1', '\\sigma_2', '\\tau_{12}']);
xgrid();
    
deformacoes_globais = [0; 0; 0];
tensoes_globais = [0; 0; 0];
h_total = numCamadas * espessura;
    
    F1 = 1/XT - 1/XC;
    F2 = 1/YT - 1/YC;
    F11 = 1/(XT*XC);
    F22 = 1/(YT*YC);
    F66 = 1/S12^2;
    F12 = -0.5 * sqrt(F11 * F22);

printf ("\n Critério Tsai-Wu\n");
for i = 1:numCamadas
    Xc = abs(XC);
    Yc = abs(YC);
    A1 = (tlocalzm(1,1,i)^2) / (XT * Xc) + (tlocalzm(2,1,i)^2 / (YT * Yc)) + ((tlocalzm(3,1,i)^2) / (S12^2)) - ((tlocalzm(1,1,i)*tlocalzm(2,1,i)) / sqrt(XT*Xc*YT*Yc));
    B1 = (((1 / XT) - (1 / Xc)) * tlocalzm(1,1,i)) + (((1 / YT) - (1 / Yc)) * tlocalzm(2,1,i));
    R3 = (-B1 / (2 * A1)) + sqrt((B1^2) / (2 * A1)^2 + (1 / A1));
    IFtsaiwu = 1 / R3;
    
    resultado = zeros(1, numCamadas);
    
for i = 1:numCamadas
    F1 = (1 / XT) + (1 / XC);
    F11 = -1 / (XT * XC);
    F2 = (1 / YT) + (1 / YC);
    F22 = -1 / (YT * YC);
    F66 = (1 / S12)^2;
    F12 = -0.5 * sqrt(F11 * F22);
    
    criterio = F1 * sigma1 + F2 * sigma2 + F11 * sigma1^2 + F22 * sigma2^2 + 2 * F12 * sigma1 * sigma2 + F66 * tau12^2;
 Fsigma(i) = 1 / criterio;

    if Fsigma(i) < 1 then
        resultado(i) = 0;
    else
        resultado(i) = 1;
    end
end


//scf(2);
//bar(resultado, 'g'); 
//title('Análise de Falha para Cada Lâmina');
//xlabel('Número da Lâmina');
//ylabel('Estado (1 = Segura, 0 = Falha)');
end

//TSAI HILL
function val=cth(x,y,j,zz)
    if x > 0 & y > 0 then
        val = (x/XT(j))^2 + (y/YT(j))^2 + x*y/XT(j)^2 - 1 + zz(j)^2/S12(j)^2
    elseif x > 0& y < 0
        val = x^2/XT(j)^2 + y^2/YC(j)^2 - x*y/XT(j)^2 - 1 + zz(j)^2/S12(j)^2
    elseif x < 0 & y < 0 // 3o q
        val = x^2/XC(j)^2 + y^2/YC(j)^2 + x*y/XC(j)^2 - 1 + zz(j)^2/S12(j)^2
    else
        val = x^2/XC(j)^2 + y^2/YT(j)^2 - x*y/XC(j)^2 - 1 + zz(j)^2/S12(j)^2
    end
endfunction;

//TSAI WU
function val=ctw(x,y,j,zz)
    val = F1(j)*x + F2(j)*y + F11(j)*x^2 + F22(j)*y^2 + 2*F12(j)*x*y - 1  + F66(j)*zz^2
endfunction

scf(3);
for j = 1:1

    x1 = [XC(j), 0, XT(j)];
    y1 = [YT(j), YT(j), YT(j)];
    y2 = [YC(j), YC(j), YC(j)];
    x3 = [XC(j), XC(j), XC(j)];
    y3 = [YT(j), 0, YC(j)];
    x4 = [XT(j), XT(j), XT(j)];
    plot(x1, y1, 'g', x1, y2, 'g', x3, y3, 'g', x4, y3, 'g');

    data = linspace(-2 * XT(j), 2 * XT(j), 500);
    zz = 0;

    contour2d(data, data, cth, [0, 0], style = color("red"));
    contour2d(data, data, ctw, [0, 0], style = color("blue"));

    xlabel("$\sigma_1$", "fontsize", 4, "color", "black");
    ylabel("$\sigma_2$", "fontsize", 4, "color", "black");

    title("Material " + string(j), "color", "black", "fontsize", 4);

    for i = 1:numCamadas
        sigma_global = Qks(:,:,i) * (epsilon0 + zm(i) * kappa);
        sigma1 = sigma_global(1)*10^3;
        sigma2 = sigma_global(2)*10^3;

        plot(sigma1, sigma2, ".r", "MarkerSize", 8);
    end

    set(gca(), "grid", [1 1]);
    fig = gca();
    fig.y_location = "origin";
    fig.x_location = "origin";
end

    disp("Matriz de Rigidez Local (Qlocal):");
    disp(Qlocal);

    disp("Submatriz A:");
    disp(A);
    disp("Submatriz B:");
    disp(B);
    disp("Submatriz D:");
    disp(D);
    disp("Deformações no plano médio (epsilon0):");
    disp(epsilon0);
    disp("Curvaturas (kappa):");
    disp(kappa);
    disp('R:', R);
    disp('r:', r);
    disp(deformacoes_momento);

    theta_range = 0:0.1:2*%pi;
    Aglobal = zeros(3, 3, length(theta_range));
    Bglobal = zeros(3, 3, length(theta_range));
    Dglobal = zeros(3, 3, length(theta_range));
    
    for a = 1:length(theta_range)
        T_rot = [cos(theta_range(a))^2, sin(theta_range(a))^2, 2*cos(theta_range(a))*sin(theta_range(a));                 sin(theta_range(a))^2, cos(theta_range(a))^2, -2*cos(theta_range(a))*sin(theta_range(a));                 -cos(theta_range(a))*sin(theta_range(a)), cos(theta_range(a))*sin(theta_range(a)), cos(theta_range(a))^2 - sin(theta_range(a))^2];
        Aglobal(:,:,a) = inv(T_rot) * A * T_rot;
        Bglobal(:,:,a) = inv(T_rot) * B * T_rot;
        Dglobal(:,:,a) = inv(T_rot) * D * T_rot;
    end
    
    scf(4);
    
    subplot(331);
    polarplot(theta_range, squeeze(Aglobal(1, 1, :)));
    title("A11 global");
    
    subplot(332);
    polarplot(theta_range, squeeze(Aglobal(1, 2, :)));
    title("A12 global");
    
    subplot(333);
    polarplot(theta_range, squeeze(Aglobal(1, 3, :)));
    title("A16 global");
    
    subplot(334);
    polarplot(theta_range, squeeze(Aglobal(2, 1, :)));
    title("A21 global");
    
    subplot(335);
    polarplot(theta_range, squeeze(Aglobal(2, 2, :)));
    title("A22 global");
    
    subplot(336);
    polarplot(theta_range, squeeze(Aglobal(2, 3, :)));
    title("A26 global");
    
    subplot(337);
    polarplot(theta_range, squeeze(Aglobal(3, 1, :)));
    title("A61 global");
    
    subplot(338);
    polarplot(theta_range, squeeze(Aglobal(3, 2, :)));
    title("A62 global");
    
    subplot(339);
    polarplot(theta_range, squeeze(Aglobal(3, 3, :)));
    title("A66 global");

    scf(5);
    
    subplot(331);
    polarplot(theta_range, squeeze(Bglobal(1, 1, :)));
    title("B11 global");
    
    subplot(332);
    polarplot(theta_range, squeeze(Bglobal(1, 2, :)));
    title("B12 global");
    
    subplot(333);
    polarplot(theta_range, squeeze(Bglobal(1, 3, :)));
    title("B16 global");
    
    subplot(334);
    polarplot(theta_range, squeeze(Bglobal(2, 1, :)));
    title("B21 global");
    
    subplot(335);
    polarplot(theta_range, squeeze(Bglobal(2, 2, :)));
    title("B22 global");
    
    subplot(336);
    polarplot(theta_range, squeeze(Bglobal(2, 3, :)));
    title("B26 global");
    
    subplot(337);
    polarplot(theta_range, squeeze(Bglobal(3, 1, :)));
    title("B61 global");
    
    subplot(338);
    polarplot(theta_range, squeeze(Bglobal(3, 2, :)));
    title("B62 global");
    
    subplot(339);
    polarplot(theta_range, squeeze(Bglobal(3, 3, :)));
    title("B66 global");

    scf(6);
    
    subplot(331);
    polarplot(theta_range, squeeze(Dglobal(1, 1, :)));
    title("D11 global");
    
    subplot(332);
    polarplot(theta_range, squeeze(Dglobal(1, 2, :)));
    title("D12 global");
    
    subplot(333);
    polarplot(theta_range, squeeze(Dglobal(1, 3, :)));
    title("D16 global");
    
    subplot(334);
    polarplot(theta_range, squeeze(Dglobal(2, 1, :)));
    title("D21 global");
    
    subplot(335);
    polarplot(theta_range, squeeze(Dglobal(2, 2, :)));
    title("D22 global");
    
    subplot(336);
    polarplot(theta_range, squeeze(Dglobal(2, 3, :)));
    title("D26 global");
    
    subplot(337);
    polarplot(theta_range, squeeze(Dglobal(3, 1, :)));
    title("D61 global");
    
    subplot(338);
    polarplot(theta_range, squeeze(Dglobal(3, 2, :)));
    title("D62 global");
    
    subplot(339);
    polarplot(theta_range, squeeze(Dglobal(3, 3, :)));
    title("D66 global");
