% Wave1d
%
% Using initial and boundary conditions provided by the user we are 
% approximating the wave equation in one dimension using an 
% iterative method accounting for insulated boundary conditions
% 
%
% Parameters
% ==========
%    kappa The initial solution to the diffusion equation 
%    x_rng The range of physical spatial boundary
%    t_rng The range of times we are solving upon
%
%    u_init The initial solution to the diffusion equation
%    u_bndry The final solution to the diffusion equation
%    du_init The inital values of he derivative   
%
%    nx Number of evenly spaced spatial points we are calculating the approx. on
%    nt Number of evenly spaced temporal points we are calculating the approx. on
%
% Return Values
% =============
%    x_out The vector x is each value of x at which a solution was found
%    t_out The vector t is the corresponding time at which solution was found
%    U_out The matrix U is the value of the solution at the corresponding position and time

function [x_out, t_out, U_out] = wave1d( c, x_rng, nx, t_rng, nt, u_init, du_init,  u_bndry )
%Argument Checking
%=================
if ~isscalar( c ) 
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument kappa is not a scalar' ) );
end
if ~all( size( x_rng ) == [1, 2] ) 
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_rng is not a 2-dimensional row vector' ) );
end
if ~all( size( t_rng ) == [1, 2] ) 
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument t_rng is not a 2-dimensional row vector' ) );
end
if ~isscalar(nx) || (nx <= 0) || (nx ~= round(nx))
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument nx is not a positive integer' ) );
end
if ~isscalar(nt) || (nt <= 0) || (nt ~= round(nt))
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument nt is not a positive integer' ) );
end
if ~isa( u_init, 'function_handle' )
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_init is not a function handle' ) );
end
if ~isa( u_bndry, 'function_handle' )
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a function handle' ) );
end
if ~isa( du_init, 'function_handle' )
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument du_init is not a function handle' ) );
end

%error checking
%==============
%
%Need to check if method will converge by checking the value of the
%Coeffecient
%Tell user an acceptable value for dt so that the ratio will be less than 1
%If the value is too high, then the solution may be unstable
delta_t = (t_rng(2)-t_rng(1))/(nt-1);
h = (x_rng(2) - x_rng(1))/(nx-1);
r = (c*delta_t/h)^2;
n_t = ceil((c*(t_rng(2)-t_rng(1))/(h))+1);
if r > 1 
    throw(MException( 'MATLAB:questionable_argument', ...
        'the ratio (c*dt/h)^2 = %f >= 1, consider using nt =  %i', r, n_t )) ;
end

%Initialization
%==============
%
% Initialize matricies that are being returned
x_out = linspace(x_rng(1), x_rng(2),  nx)';
t_out = linspace(t_rng(1), t_rng(2), nt);
U_out = zeros(nx, nt);
for i = 2:nt
    a = u_bndry(t_out(i));
    U_out(1,i) = a(1);
    U_out(nx,i) =a(2);
end
for i = 1:nx
    U_out(i,1) = u_init(x_out(i));
end

%solving for t2
%==============
%
% Solving the second column of the U_out matrix
for i = 1:nx-2
    U_out(i+1,2) = U_out(i+1,1) + delta_t*du_init(x_out(i));
end
if isnan(U_out(1,2)) ==1
    U_out(1,2) = U_out(2,2);
end
if isnan(U_out(nx,2)) == 1
    U_out(nx,2) = U_out(nx-1,2);
end
%Solving
%===========
%
%iterate through every time column solving them based off of the previous
%ones
%if boundary is insulated, adjust it accordingly
for k = 2:nt-1
    for i = 2:nx-1
        U_out(i, k+1) = 2*U_out(i,k) -U_out(i,k-1) + r*(U_out(i-1,k) -2*U_out(i,k) + U_out(i+1,k));
    end
    if isnan(U_out(1,k)) ==1
        U_out(1,k) = U_out(2,k);
    end
    if isnan(U_out(nx,k)) == 1
        U_out(nx,k) = U_out(nx-1,k);
    end
end
end
