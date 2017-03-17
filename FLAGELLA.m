clear all
clc

%Flagella length and partitioning
L = 1;
s = 0.1; Z%forcing point spacing
N = L/s + 1; % # of forcing points
x = 0:s:L;
%Time scale of motion and partitioning
T = 10;
tic = 0.01;
ticks = T/tic + 1;
t = 0:tic:T;
t = 0;
%Prescribed wave function for Flagella's motion
Yc = zeros(N,1);
Yc(1) = 10*s;
k = 10;
w = 0.001;
B = w/k;
Yc = Yc(1)*cos(k*x - w*t);

%Filling in coefficient matrices a,b,c
a = zeros(N,N);
b = zeros(N,N);
c = zeros(N,N);
for i = 1:N
    for j = 1:N
        y = Yc;
        r = sqrt((x(j) - x(i))^2 + (y(j) - y(i))^2);
        a(i,j) = (x(j) - x(i))*(y(j) - y(i))*(mob_par_approx(r) - mob_perp_approx(r))/...
            r^2;
        b(i,j) = (mob_par_approx(r)*((y(j) - y(i))^2) - mob_perp_approx(r)*(((x(j) - x(i))^2)))/...
            r^2;
        c(i,j) = (mob_par_approx(r)*(((x(j) - x(i))^2) + mob_perp_approx(r)*((y(j) - y(i))^2)))/...
            r^2;
    end
end

%Filling in MATRIX with coefficients of unknown variables fx1,fx2,...,fxN 
%and fy1,fy2,...,fyN, and VyL and VxL (the flagella velocity components in the lab frame) in (2N+2)x(2N+2) linear
%system
M = zeros(2*N + 2, 2*N + 2);
%Upper Left corner of MATRIX
for i = 1:N
    for j = 1:N
        M(i,j) = a(i,j);
    end
end

%Upper Right corner of MATRIX
for i = 1:N
    for j = (N+1):2*N
        M(i,j) = b(i,j-N);
    end
end

%Bottom Left corner of MATRIX
for i = (N+1):2*N
    for j = 1:N
        M(i,j) = c(i-N,j);
    end
end

%Bottom Right corner of MATRIX
for i = (N+1):2*N
    for j = (N+1):2*N
        M(i,j) = a(i-N,j-N);
    end
end

%Filling Last 2 columns except for last 2 rows
%Upper Half of Last 2 columns (j = 2N + 1, 2N + 2)
for i = 1:N 
    M(i,2*N + 1) = 0;
    M(i,2*N + 2) = -1;
end

%Bottom Half of Last 2 columns (j = 2N + 1, 2N + 2)
for i = (N+1):2*N
    M(i,2*N + 1) = -1;
    M(i,2*N + 2) = 0;
end

%Filling Last 2 Row except for last 2 columns
%Left Half of Last 2 Rows (i = 2N + 1, 2N + 2)
for j = 1:N
    M(2*N + 1, j) = 1;
    M(2*N + 2, j) = 0;
end

%Right Half of Last 2 Rows (i = 2N + 1, 2N + 2)
for j = (N+1):2*N
    M(2*N + 1, j) = 0;
    M(2*N + 2, j) = 1;
end

%Filling in the very tip of the bottom right corner
%[(i,j) = (2N + 1, 2N + 1),(2N + 1, 2N +2), (2N + 2, 2N + 1), (2N + 2, 2N + 2)]
%for i = (2*N + 1):(2*N + 2)
%    for j = (2*N + 1):(2*N + 2)
%        M(i,j) = 0;
%    end
%end

%Filling column vector V, with known constant terms, from eq. V = MS
%(2N)Boundary conditons and (2) zero net force eqs. give us these.
V = zeros(2*N + 2, 1);
for i = 1:N
    V(i) = w*Yc(1)*sin(k*x(i) - w*t);
end
for i = (N+1):2*N
    V(i) = B;
end

%LINEAR EQUATION V = MS solved for S, using MATLAB's MATRIX LINEAR EQUATION SOLVER
S = zeros(2*N + 2, 1);
S = linsolve(M,V);
