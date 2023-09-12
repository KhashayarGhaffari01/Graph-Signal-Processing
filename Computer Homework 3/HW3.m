%file location should be corrected.
v = csvread("Graph-Signal-Processing\Computer-Homework-3\distance.csv");
realDistance = zeros(31);
for i=1:31
    for j=1:31
        if i == j
            realDistance(i,j) = 0;
        else
        realDistance(i,j) = getDistance(v(i,1), v(i,2), v(j,1), v(j,2));
        end
    end
end
W = (realDistance*(-0.95/max(realDistance, [], 'all')) + 1)-eye(31);
G = gsp_graph(W,gsp_sensor(31).coords);
G = gsp_compute_fourier_basis(G);
%gsp_plot_graph(G);
temp_signal = v(:,3);
[Graphs, m, y, x] = MyPyramidAnalysis(G.L, temp_signal, 3, 0.005, 1);
x_baz = MyPyramidSynthesis(Graphs, y, m, x{3}, 3, G, 0.005);


function G = lap2Graph(gLap)
G = gsp_graph(-gLap + diag(diag(gLap)));
end


function n = GraphSize(G)
v = size(G.W);
n = v(1);
end


function V_1 = MyVertexSelection(G)
n = GraphSize(G);
u_max = G.U(:, n);
V_1 = u_max;
for i=1:n
    if V_1(i) < 0
        V_1(i) = -i;
    else
        V_1(i) = i;
    end
end
V_1(V_1<0) = [];
end


function G_reduced = MySKReduction(G, V_1)
V_1_complement = Vcomplement(V_1,GraphSize(G));
kronL = G.L(V_1, V_1) - (G.L(V_1, V_1_complement)*(G.L(V_1_complement, V_1_complement)^(-1))*(G.L(V_1_complement, V_1)));
kronW = -kronL + diag(diag(kronL));
SKReducedW = zeros(length(V_1));

d_R_G = zeros(length(V_1));
S = (1:length(V_1))';
for i=1:length(V_1)
    for j=1:length(V_1)
       d_R_G(i,j) = (kroneckerDelta(S,sym(i))-kroneckerDelta(S,sym(j)))'*(pinv(kronL))*(kroneckerDelta(S,sym(i))-kroneckerDelta(S,sym(j))); 
    end
end

P_distribution = zeros(length(V_1));
for i=1:length(V_1)
    for j=1:length(V_1)
      P_distribution(i,j) = (d_R_G(i,j)*kronW(i,j))/(sum(sum(d_R_G.*kronW)));
    end
end

Q = ceil(log10(GraphSize(G))*4*GraphSize(G));
for q=1:Q
    P_distribution1 = reshape(P_distribution,[1, length(V_1)^2]);
    y = randsample(length(V_1)^2,1,true,P_distribution1);
    j = rem(y, length(V_1));
    i = ((y-j)/length(V_1))+1;
    if j == 0
        j = length(V_1);
    end
    if i ~= j
    SKReducedW(i,j) = SKReducedW(i,j) + kronW(i,j)/(Q*P_distribution(i,j));
    end
end
G_reduced = gsp_graph(SKReducedW).L;
end


function xhat = GFT(G,x)
xhat = (G.U)^(-1)*x;
end


function x = RGFT(G,xhat)
x = G.U * xhat;
end


function V_1_Complement = Vcomplement(V,n)
S = (1:n)';
V_1_Complement = setdiff(S, V);
end


function y = MyHfilter(G, x)
G = gsp_compute_fourier_basis(G);
hhat = ((diag(G.e))*2 + eye(GraphSize(G)))^(-1);
y = RGFT(G, hhat*GFT(G, x));
end


function xDsample = MyDS(V_1, x)
xDsample = zeros(length(V_1), 1);
for i=1:length(V_1)
    xDsample(i) = x(V_1(i));
end
end


function f_interp = MyInterpolate(G, V_1, f, epsilon, flag)
if flag == 0
    phi_V_1 = zeros(GraphSize(G), length(V_1));
    for j=1:length(V_1)
        phi_j = zeros(GraphSize(G),1);
        for l=1:GraphSize(G)
            phi_j = phi_j + ((1/(G.e(l)+epsilon))*conj(G.U(j,l))*G.U(:,l));
        end
        phi_V_1(:, j) = phi_j;
    end
    epsilonL = G.L + eye(GraphSize(G))*epsilon;
    V_1_complement = Vcomplement(V_1,GraphSize(G));
    alphaStar = (epsilonL(V_1,V_1)-(epsilonL(V_1, V_1_complement))*((epsilonL(V_1, V_1_complement))^(-1))*(epsilonL(V_1, V_1_complement)))*f;
    f_interp = phi_V_1*alphaStar;
end
if flag == 1
    f_interp = G.U * pinv(G.U(V_1, :)) * f;
end
end


function [G_j1Lap, m_j, y_j, x_j1] = MyAnalysis(G_jLap, x_j, epsilon, flag)
G_j = lap2Graph(G_jLap);
G_j = gsp_compute_fourier_basis(G_j);
m_j = MyVertexSelection(G_j);
G_j1Lap = (MySKReduction(G_j, m_j));
x_j1 = MyDS(m_j, MyHfilter(G_j, x_j));
y_j = x_j - MyInterpolate(G_j, m_j, x_j1, epsilon, flag);
end


function [GraphsLapCells, m, y, x] = MyPyramidAnalysis(G_0, x_0, N, epsilon, flag)
GraphsLapCells = cell(N,1);
y = cell(N,1);
m = cell(N,1);
x = cell(N,1);
[GraphsLapCells{1}, m{1}, y{1}, x{1}] = MyAnalysis(G_0, x_0, epsilon, flag);
for i=2:N
    [GraphsLapCells{i}, y{i}, m{i}, x{i}] = MyAnalysis(GraphsLapCells{i-1}, x{i-1}, epsilon, flag);
end
end


function x_j = MySynthesis(GL_j, x_j1, y_j1, m_j1, epsilon)
G_j = gsp_graph(-GL_j+diag(diag(GL_j)));
G_j = gsp_compute_fourier_basis(G_j);
phi_m_j1 = zeros(GraphSize(G_j), length(m_j1));
    for j=1:length(m_j1)
        phi_j = zeros(GraphSize(G_j),1);
        for l=1:GraphSize(G_j)
            phi_j = phi_j + ((1/(G_j.e(l)+epsilon))*conj(G_j.U(j,l))*G_j.U(:,l));
        end
        phi_m_j1(:, j) = phi_j;
    end
   
x_j = [phi_m_j1 * ((phi_m_j1(m_j1, :))^(-1)), eye(GraphSize(G_j))] *[x_j1; y_j1];
end


function x = MyPyramidSynthesis(GraphsLapCells, y, m, x_N, N, G, epsilon)
x = x_N;
for i=N:2:-1
    x = MySynthesis(GraphsLapCells{i-1}, x, y{i}, m{i}, epsilon);
end
x = MySynthesis(G.L, x, y{1}, m{1}, epsilon);
end


function D = getDistance(lat1, lon1, lat2, lon2)
rad = pi / 180;
radius = 6371; %earth radius in kilometers
D = abs(acos(sin(lat2 * rad) * sin(lat1 * rad) + cos(lat2 * rad) * cos(lat1 * rad) * cos(lon2 * rad - lon1 * rad)) * radius); %result in Kilometers
end