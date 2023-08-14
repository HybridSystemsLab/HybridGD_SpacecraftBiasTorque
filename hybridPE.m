classdef hybridPE < HybridSystem
    
    properties(SetAccess=immutable)
        theta
        gammac
        gammad
        lambdac
        lambdad
        Kp
        Kd
        Js
        Jw
        delta
        Mt
        Omega_max
        taus
        u
    end

    properties 
        %
    end
    
    methods 
        function this = hybridPE(parameters)
            state_dim = 13; % (x,xdot,thetahat,psi,eta,u)
            this = this@HybridSystem(state_dim);
            
            this.theta = parameters.theta;
            this.gammac = parameters.gammac;
            this.gammad = parameters.gammad;
            this.lambdac = parameters.lambdac;
            this.lambdad = parameters.lambdad;
            this.Kp = parameters.Kp;
            this.Kd = parameters.Kd;
            this.Js = parameters.Js;
            this.Jw = parameters.Jw;
            this.delta = parameters.delta;
            this.Mt = parameters.Mt;
            this.Omega_max = parameters.Omega_max;
            this.taus = parameters.taus;
            this.u = parameters.u;
        end

        function xidot = flowMap(this, xi, t, j)            
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            thetahat = xi(5);
            psi = xi(6:9);
            eta = xi(10:13);
            
            z = xi(1:4);

            %%            
            alpha = -this.Kp*(this.u - x) + this.Kd*xdot + thetahat;

            fc = [xdot
                  -this.Js\alpha
                  this.Jw\alpha
                  1];
            phic = [0
                    this.Js\1
                    0
                    0];

            zdot = fc + phic*this.theta;

            y_c = z + eta;
            thetahatdot = this.gammac*psi.'*(y_c - psi*thetahat);
            psidot = -this.lambdac*psi + phic;
            etadot = -this.lambdac*(z + eta) - fc;

            xidot = [zdot; thetahatdot; psidot; etadot];
            
        end

        function xiplus = jumpMap(this, xi, t, j)
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            thetahat = xi(5);
            psi = xi(6:9);
            eta = xi(10:13);
            
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

            %%
%             y_d = zplus - fd;
%             thetahatplus = thetahat + (phid.'/(this.gammad + norm(phid)^2))*(y_d - phid*thetahat);
%             psiplus = (1 - this.lambdad)*psi;
%             etaplus = (1 - this.lambdad)*(z + eta) - zplus;
            
            psiplus = (1-this.lambdad)*psi + phid;
            etaplus = (1-this.lambdad)*(z + eta) - fd;
            yplus = zplus + etaplus;
            thetahatplus = thetahat + (psiplus.'/(this.gammad + norm(psiplus)^2))*(yplus - psiplus*thetahat);
            % thetahatplus = thetahat;

            %%
            xiplus = [zplus; thetahatplus; psiplus; etaplus];
        end
        
        function inC = flowSetIndicator(this, xi, t, j)
            x = xi(1);
            xdot = xi(2);
            omega = xi(3);
            tau = xi(4);
            thetahat = xi(5);
            psi = xi(6:9);
            eta = xi(10:13);
            
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
            thetahat = xi(5);
            psi = xi(6:9);
            eta = xi(10:13);
            
            z = xi(1:4);
            
            if omega >= this.Omega_max && tau >= this.taus
                inD = 1;
            else
                inD = 0;
            end
        end
    end
end
