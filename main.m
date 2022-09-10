clear all;
digits(4);

% Quadrotor Constants - Section IV 
% https://ieeexplore.ieee.org/document/5717652
% Mass
quad.m = 4.34;
% Acceleration due to gravity
quad.g = 9.98;
% distance between Centre of mass to center of rotor 
quad.d = 0.315;
% Drag Constant
quad.ctf = 8.004e-4;
% Ineteria Tensor
quad.J = diag([0.0820,0.0845,0.1377]);
% Controller parameters
quad.kx = 16 * quad.m;
quad.kv = 5.6 * quad.m;
quad.kr = 8.81;
quad.kw = 2.54;
% Axes
quad.e1 = [1;0;0];
quad.e2 = [0;1;0];
quad.e3 = [0;0;1];


% Initial Conditions
% Position
x0 = [0;0;0];
% Velocity
v0 = [0;0;0];
% Rotaion matrix
R0 = eye(3);
% Angular Velocity
W0 = [0;0;0];
t0 = 0;


% 4x4 matrix from Eq 1
% [ 1   1   1      1]  
% | 0  -d   0      d|  
% | d   0   -d     0| 
% [-ctf ctf -ctf ctf] 
quad.T = [1, 1, 1, 1; 0, -quad.d, 0, quad.d; quad.d, 0, -quad.d, 0;
    -quad.ctf, quad.ctf, -quad.ctf, quad.ctf];


syms a
% Section IV 
% Desired trajectory
quad.xd = [0.4*a; 0.4*sin(pi*a); 0.6*cos(pi*a)];
% Desired Heading
% heading direction of the quadrotor UAV in the plane normal to b3d
quad.b1d = [cos(pi*a);sin(pi*a);0];

% desired position - numerical derivatives
quad.x1dot = diff(quad.xd,1);
quad.x2dot = diff(quad.xd,2);
quad.x3dot = diff(quad.xd,3);
quad.x4dot = diff(quad.xd,4);

% desired body 1 axis - numerical derivatives
quad.b1d_dot = diff(quad.b1d,1);
quad.b1d_2dot = diff(quad.b1d,2);

y0 = [x0',v0',W0',R0(1,:),R0(2,:),R0(3,:)];
dyn = @(t,y)quadDynamics(t,y,quad);

[t,y] = ode45(dyn,[0 4],y0);

% Plot data here
%% 
% https://github.com/sir-avinash/geometry-toolbox
function hat_x = hat(x)
% Mapping vector to Skew Symmetric Matrix 
hat_x = [0, -x(3), x(2);
         x(3), 0, -x(1);
         -x(2), x(1), 0];
end
%% 
% https://github.com/sir-avinash/geometry-toolbox
function vee_x = vee(x)
% Mapping Skew Symmetric Matrix to vector
vee_x = [x(3,2); x(1,3); x(2,1)];
end

%% 
function dy = quadDynamics(t,y,quad)
x = y(1:3);
v = y(4:6);
W = y(7:9);
R = [y(10:12)';y(13:15)';y(16:18)'];
controllerOuput = controller(t,y,quad); 
f = controllerOuput(1);
M = controllerOuput(2:4);
dx = v;
dv = quad.g * quad.e3 - (1/quad.m) * f * R * quad.e3;
dW = (quad.J)\(M - cross(W, quad.J * W));
dR = R * hat(W);
dy = [dx',dv',dW',dR(1,:),dR(2,:),dR(3,:)]';
end
%% 
function controllerOut = controller(t,y,quad)
digits(4)

x = y(1:3);
v = y(4:6);
W = y(7:9);
R = [y(10:12)';y(13:15)';y(16:18)'];

xd = vpa(subs(quad.xd,t));
vd = vpa(subs(quad.x1dot,t));
x2dot = vpa(subs(quad.x2dot,t));
x3dot = vpa(subs(quad.x3dot,t));
x4dot = vpa(subs(quad.x4dot,t));

b1d = vpa(subs(quad.b1d,t));
b1d_dot = vpa(subs(quad.b1d_dot,t));
b1d_2dot = vpa(subs(quad.b1d_2dot,t));

% Eq 6
ex = x - xd;
% Eq 7
ev = v - vd;

% Eq 15 
A = vpa(-quad.kx * ex - quad.kv * ev - quad.m * quad.g * quad.e3 + quad.m * x2dot);
% Eq 15 Total Thrust
f = dot(-A, R * quad.e3);

% Eq 12
% desired direction of the third body-fixed frame b3d
b3d = -A/norm(A);

% error in accleration
% Eq 3
dv =  quad.g * quad.e3 - (1/quad.m) * f * R * quad.e3;
% actual - desired
ea = dv - x2dot;

% Computing b2d and b1d
% b2d = (b3d x b1d )/||b3d x b1d||
C = cross(b3d, b1d);
b1d = -(1/norm(C)) * cross(b3d,C);
b2d = C/norm(C);

% Desired Attitude
Rd = [b1d b2d b3d];

% https://arxiv.org/pdf/1003.2005v3.pdf
% Time derivative of body axis - First
dA = -quad.kx*ev - quad.kv*ea + quad.m * x3dot;

dR = R * hat(W);
df = -dot(dA, R * quad.e3) - dot(A, dR * quad.e3);
ddv = (1/quad.m) * (df * R * quad.e3 + f * dR * quad.e3);
% Jerk error
ej = ddv - x3dot;

b3d_dot = -dA/norm(A) + (dot(A,dA)/norm(A)^3)*A;
dC = cross(b3d_dot,b1d) + cross(b3d,b1d_dot);
b2d_dot = dC/norm(C) - (dot(C,dC)/norm(C)^3)*C; 
b1d_dot = cross(b2d_dot,b3d) + cross(b2d,b3d_dot);

% Time derivative of body axis - Second
ddA = - quad.kx * ea - quad.kv * ej + quad.m * x4dot;
b3d_2dot = - ddA/norm(A) + (2/norm(A)^3)*dot(A,dA)*dA + ((norm(dA)^2 + dot(A,ddA))/norm(A)^3)*A - (3/norm(A)^5)*(dot(A,dA)^2)*A;
ddC = cross(b3d_2dot,b1d) + cross(b3d_dot,b1d_dot) + cross(b3d_dot,b1d_dot) + cross(b3d,b1d_2dot);
b2d_2dot = ddC/norm(C) - (2/norm(C)^3)*dot(C,dC)*dC - ((norm(dC)^2 + dot(C,ddC))/norm(C)^3)*C + (3/norm(C)^5)*(dot(C,dC)^2)*C;  
b1d_2dot = cross(b2d_2dot,b3d) + cross(b2d_dot,b3d_dot) + cross(b2d_dot,b3d_dot) + cross(b2d,b3d_2dot);

dRd = [b1d_dot b2d_dot b3d_dot];
ddRd = [b1d_2dot b2d_2dot b3d_2dot];
Wd = vee(Rd.'* dRd);
dWd = vee(Rd.'*ddRd - hat(Wd)*hat(Wd));   

% Eq 10
% Attitude tracking error er
er = 0.5 * vee(Rd.'* R - R.'*Rd);

% Eq 11
% tracking error for angular velocity
ew = W - R.'* Rd * Wd;

% Eq 16
M = - quad.kr * er - quad.kw * ew + cross(W, quad.J * W) - quad.J * (hat(W) * R.' * Rd * Wd - R.' * Rd * dWd);

controllerOut = double([f;M]);
end
