function [z_WLC, g_WLC_kbT] = WLC_z_G_twist(f,L,dLk_DNA)                                   

% OUTCOME: extension and internal free energy for stretching and twisting the WLC by
% magnetic tweezers

% INPUT PARAMETERS: f - force ramp (pN), L - contour length (nm), dLk_DNA -number of applied
% turns by the magnet

% Artur Kaczmarczyk, akaczmarczyk88@gmail.com
% June 2018 


                                                                    
kbT = 4.114;                                                                 % Boltzmann constant in room temperature (pNnm);
P = 50;                                                                      % persistence length (nm);
S = 1000;                                                                    % stretch modulus (pN)
C = 100;                                                                     % torsional stiffness (nm)
Ct = C * (1 - C ./ (4 .* P) .* sqrt(kbT ./ (P .* f)));                       % force-dependent effective torsional stiffness (nm)
Cp = 24;                                                                     % torsional stiffness of plectonemic DNA (nm)

tau = Ct .* (2 .* pi .* dLk_DNA ./ L) .* kbT;                                % DNA restoring torque (pN*nm)
tau_melt = -11;                                                              % melting torque (pN*nm)

z_WLC = L.*((1 - 0.5 .* sqrt(kbT./(P .* f))) - (C.^2 ./ 16) .* (2 .* pi .* dLk_DNA ./ L).^2 .* ( kbT ./ (P.*f)) .^ (3/2));  % twisted DNA length in nanometers

g = (f - sqrt(f.*kbT./P));                                                   % "stretching part" for energy term (same as 1 -0.5sqrt(kbT/P*f) used later) 
%z_WLC = (L*(1-0.5*sqrt(kbT./(f.*P)))+f./S);                                 % no twist component                                 
%g_WLC =  z_WLC.*f - L.*f.* (1-sqrt(kbT./(f.*P))+f./(2.*S));                 % no twist component


%% border conditions for plectonemes - different energy calculation for twisted and plectonemic DNA (free energy of plectonemic DNA increases linearly with twist; the gain in free energy per turn is: )

% positive twist
if dLk_DNA >= 0                                                                 
    
for i = 1:length(f)
                                         
    tau_buck(i) = sqrt (2 .* kbT .* Cp .* g(i) ./ (1 - Cp./Ct(i)));        % buckling torque is force-dependent
        
    if tau > tau_buck(i)                                                   % torque is exerted due to the applied twist, if it exceeds the buckling torque for given force, plectonemes are formed
            
        
            z_slope(i) = (2 .* pi .* (1 - 0.5 .* sqrt(kbT./(P .* f(i))) - C.^2./(16 .* Ct(i).^2) .* (kbT ./(P.*f(i))).^(3/2).*(tau_buck(i)./kbT).^2))/((tau_buck(i)./kbT).*(1./Cp - 1./Ct(i)));  % slope that described DNA buckling (nm/turn)
            
            z_WLC(i) = z_WLC(i) - abs(dLk_DNA) .* z_slope(i);
            
            g_WLC(i) =   z_WLC(i).*f(i) - L.*f(i).* (1-sqrt(kbT./(f(i).*P))+f(i)./(2.*S)) + dLk_DNA .* 2 .* pi .* tau_buck(i);
            
            if z_WLC(i) < 0;

                z_WLC(i) = 0;
                g_WLC(i) =  g_WLC(i);
                
            else
                
                z_WLC(i) = z_WLC(i);
                g_WLC(i) =  g_WLC(i);
            
            end
    
       
                 
    else
        
            z_slope(i) = 0;                                                 % if the buckling torque has not been reached, the extension follows the WLC model
            
            z_WLC(i) = z_WLC(i);
           
            g_WLC(i) =   z_WLC(i).*f(i) - L.*f(i).* (1-sqrt(kbT./(f(i).*P))+f(i)./(2.*S)) + L .*( (0.5 * kbT .* Ct(i)) .* (2 .* pi .* dLk_DNA ./ L).^2); % energy of twisted DNA
            
            
    end

    end

else
    
    %% negative twist
    
    for i = 1:length(f)
        
        if tau > tau_melt      
        
            z_slope(i) = (2 .* pi .* (1 - 0.5 .* sqrt(kbT./(P .* f(i))) - C.^2./(16 .* Ct(i).^2) .* (kbT ./(P.*f(i))).^(3/2).*(tau_buck(i)./kbT).^2))/((tau_buck(i)./kbT).*(1./Cp - 1./Ct(i)));  % slope that described DNA buckling (nm/turn)
            
            z_WLC(i) = z_WLC(i) - abs(dLk_DNA) .* z_slope(i);
            
            g_WLC(i) =   z_WLC(i).*f(i) - L.*f(i).* (1-sqrt(kbT./(f(i).*P))+f(i)./(2.*S)) + dLk_DNA .* 2 .* pi .* tau_buck(i);
            
            if z_WLC(i) < 0;

                z_WLC(i) = 0;
                g_WLC(i) =  g_WLC(i);
                
            else
                
                z_WLC(i) = z_WLC(i);
                g_WLC(i) =  g_WLC(i);
            
            end
            
            
            
        else
         
            z_WLC(i) = z_WLC(i);
           
            g_WLC(i) =   z_WLC(i).*f(i) - L.*f(i).* (1-sqrt(kbT./(f(i).*P))+f(i)./(2.*S)) + dLk_DNA .* 2 .* pi .* tau_melt; % energy of twisted DNA
    
            if z_WLC(i) < 0;

                z_WLC(i) = 0;
                g_WLC(i) =  g_WLC(i);
                
            else
                
                z_WLC(i) = z_WLC(i);
                g_WLC(i) =  g_WLC(i);
            
            end
            
            
            
            
        end
    end
    
end

    
g_WLC_kbT = g_WLC./kbT;                                                     % internal free energy (kbT)


%% plotting for checking purposes; switch off for regular use    

                                                      
figure (1)
plot(z_WLC,f);
hold on;
legend ([num2str(dLk_DNA) ' turns'],'Location','southeast')
xlabel ('extension (um)')
ylabel ('force (pN)')

figure(2)
plot(g_WLC_kbT,f);
legend ([num2str(dLk_DNA) ' turns'],'Location','northwest')
xlabel ('energy (kT)')
ylabel ('force (pN)')
hold on;

%plot(tau_buck,f);

end
