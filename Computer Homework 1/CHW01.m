%Problem 2
G1 =  gsp_graph([0 1.6 2.4; 1.6 0 0.8; 2.4 0.8 0], [1 1.7; 0 0; 2 0]);
figure(1);
gsp_plot_graph(G1);

G2 = gsp_graph([0 0.7 1.1 2.3; 0.7 0 0 0; 1.1 0 0 0; 2.3 0 0 0], [0 0; -1 -1.8; 1 -1.8; 0 2]);
figure(2);
gsp_plot_graph(G2);


%Problem 3
Gs = gsp_graph_product(G1, G2);
Gs.coords = [26 -97; -26 -97; -71 -71; -97 -26; -97 26; -71 71; -26 97; 26 97; 71 71; 97 26; 97 -26; 71 -71];
figure(3);
gsp_plot_graph(Gs);

param.rule = 'kronecker';
Gt = gsp_graph_product(G1, G2, param);
Gt.coords = [26 -97; -26 -97; -71 -71; -97 -26; -97 26; -71 71; -26 97; 26 97; 71 71; 97 26; 97 -26; 71 -71];
figure(4);
gsp_plot_graph(Gt);


%Problem 4
MyG = Gs;
MySignal = 20*rand(12, 1)-10;
figure(5);
gsp_plot_signal(MyG, MySignal);


%Problem 5
MyG = gsp_compute_fourier_basis(MyG);
disp(MyG.U);
disp(MyG.e);

figure(6);
gsp_plot_signal_spectral(MyG,MySignal);


%Problem 6
figure(7);
gsp_plot_signal(MyG, MyG.U(:,2));
figure(8);
gsp_plot_signal(MyG, MyG.U(:,3));
figure(9);
gsp_plot_signal(MyG, MyG.U(:,11));
figure(10);
gsp_plot_signal(MyG, MyG.U(:,12));


%Problem 7
GL = gsp_logo;
signal = zeros(1130,1);
for i = 1:1130
    if GL.coords(i , 1) > 400
        signal(i) = -1;
    elseif GL.coords(i , 1) < 200
            signal(i) = -0.5;
        else 
            signal(i) = 1;
    end
end

figure(11);
gsp_plot_signal(GL, signal);


%Problem 8
GL = gsp_compute_fourier_basis(GL);
disp(GL.U);
disp(GL.e);


%Problem 9 & 10
twoDimCoord = zeros(1130, 2);
for i = 1:1130
    twoDimCoord(i, :) = [GL.U(i, 2) GL.U(i,3)];
end
figure(12);
gsp_plot_graph(gsp_graph(zeros(1130), twoDimCoord));


%Problem 13
threeDimCoord = zeros(1130, 3);
for i = 1:1130
    threeDimCoord(i, :) = [GL.U(i, 2) GL.U(i,3) GL.U(i, 4)];
end
figure(13);
gsp_plot_graph(gsp_graph(zeros(1130), threeDimCoord));