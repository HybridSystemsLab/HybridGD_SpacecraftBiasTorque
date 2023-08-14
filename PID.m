classdef PID < HybridSystem
    
    properties(SetAccess=immutable)
        theta
        gammac
        gammad
        lambdac
        lambdad
        Kp
        Kd
        Ki
        Js
        Jw
        delta
        Mt
        Omega_max
        taus
        u
        intInhibit
    end

    properties 
        %
    end
    
    methods 
        function this = PID(parameters)
            state_dim = 5; % (x,xdot,thetahat,psi,eta,u)
            this = this@HybridSystem(state_dim);
            
            this.theta = parameters.theta;
            this.gammac = parameters.gammac;
            this.gammad = parameters.gammad;
            this.lambdac = parameters.lambdac;
            this.lambdad = parameters.lambdad;
            this.Kp = parameters.Kp;
            this.Kd = parameters.Kd;
            this.Ki = parameters.Ki;
            this.Js = parameters.Js;
            this.Jw = parameters.Jw;
            this.delta = parameters.delta;
            this.Mt = parameters.Mt;
            this.Omega_max = parameters.Omega_max;
            this.taus = parameters.taus;
            this.u = parameters.u;
            this.intInhibit = 500;
        end

        function xidot = flowMap(this, xi, t, j)            
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            xerrInt = xi(5);
            
            z = xi(1:4);

            %%            
            alpha = -this.Kp*(this.u - x) + this.Kd*xdot - this.Ki*xerrInt;

            fc = [xdot
                  -this.Js\alpha
                  this.Jw\alpha
                  1];
            phic = [0
                    this.Js\1
                    0
                    0];

            zdot = fc + phic*this.theta;
            
            if (tau <= this.intInhibit) && (j >= 1)
                xerr = 0;
            else
                xerr = this.u - x;
            end

            xidot = [zdot; xerr];
            
        end

        function xiplus = jumpMap(this, xi, t, j)
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            xerrInt = xi(5);
            
            z = xi(1:4);
            
            %%
            fd = [x
                  xdot + this.Js\this.delta*this.Mt
                  omega
                  0];
            phid = [0
                    this.Js\this.delta
                    0
                    0];

            zplus = fd + phid*this.theta;

            xerrIntplus = xerrInt;

            %%
            xiplus = [zplus; xerrIntplus];
        end
        
        function inC = flowSetIndicator(this, xi, t, j)
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            xerrInt = xi(5);
            
            z = xi(1:4);
            
            if omega <= this.Omega_max || tau <= this.taus
                inC = 1;
            else
                inC = 0;
            end
        end

        function inD = jumpSetIndicator(this, xi, t, j)
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            xerrInt = xi(5);
            
            z = xi(1:4);
            
            if omega >= this.Omega_max && tau >= this.taus
                inD = 1;
            else
                inD = 0;
            end
        end
    end
end
