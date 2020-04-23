function [Y] = RK4(iteration,Y,A,F,dt)
    % Runge-Kutta of order 4
    k_1 = iteration(Y,A,F);
    k_2 = iteration(Y+0.5*dt*k_1,A,F);
    k_3 = iteration(Y+0.5*dt*k_2,A,F);
    k_4 = iteration(Y+k_3*dt,A,F);

    % output
    Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
end