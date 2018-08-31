function [z_final, f] = fit_FE_fiber_twist(f,dLk,CL,Nnuc,Ntet,NRL,kf,k22,degeneracy,dG1,dG2)  


%   PLOTTING CHROMATIN FORCE - EXTENSION (torsionally constrained)

%   based on: Meng, Andresen, van Noort - NAR, 2015
%   modification with respect to the paper: 
%   - transition from singly wrapped to extended state is described only with WLC function
%   - added twist dependence (torsional energy factor and change in linking number induced by fiber unfolding)

%   input parameters: 
%   force ramp (pN), number of applied turns, contour length (nm), number of
%   nucleosomes, number of tetrasomes, Nucleosome Repeat Length (bp), fiber stiffness (pN/nm), 
%   twist modulus (pN * nm^2), degeneracy, free energy of unstacking (dG1), free energy of the intermediate transition (dG2)

%   Artur Kaczmarczyk, akaczmarczyk88@gmail.com, July 2018

%% TO DO: FIX THE LEVEL OF PLATEAU AFTER APPLYING TWIST

%%

%f = (0.03:0.01:7);
kbT = 4.114;                                     % Boltzmann constant in room temperature (pNnm);
%Nnuc = 15;                                      % number of assembled nucleosomes
%Ntet = 0;                                       % number of tetrasomes
%NRL = 167;                                      % Nucleosome Repeat Length (bp)
%CL = 4535;                                      % length of the DNA template used for reconstitution(bp)
Lextended = 5;                                   % additional gain in extension after unstacking and 1st nucleosomal turn unwrapping (nm/nuc)
%kf = 1;                                         % stiffness of a folded chromatin fiber per nucleosome (pN/nm)
Lwrap = 89;                                      % length of DNA unwrapped from the first nucleosomal turn (bp)
%degeneracy = 0;                                

%dG1 = 22;                                       % free energy of unstacking (dG1)
%dG2 = 11;                                       % free energy of the intermediate transition (dG2)
dG3 = 100;                                       % high force transition (not in equlibrium, value just for plotting purposes)


L_tether = CL - Nnuc * NRL;                      % length of bare DNA handles (bp)


P = 50;                                                                      % persistence length of DNA (nm);
S = 1000;                                                                    % stretch modulus of DNA (pN)
C = 100;                                                                     % torsional stiffness of DNA (nm)
ftwisted = 0.5;                                                              % force at which the twist was applied (pN)
Ct = C * (1 - C ./ (4 .* P) .* sqrt(kbT ./ (P .* ftwisted)));                % force-dependent effective torsional stiffness of DNA (nm)


z0_f = 1.7;                                                                  % length of the fiber in the absence of force (nm/nuc)
%k11 = 1.2;                                                                  % stretch spring constant per nucleosome (pN/nm);
k11N = kf ./ ((Nnuc-Ntet).*z0_f);                                            % 
%k22 = 10;                                                                   % twist modulus per nucleosome(pN*nm);
k22N = k22 ./ ((Nnuc-Ntet).*z0_f);                                           % 
k12 = 0.02;                                                                  % twist-stretch coupling factor (pN/rad);
k12N = k12 ./ ((Nnuc-Ntet).*z0_f);                                           % 
Lkf0 = -0.15;



%% to get  dLk_f

%{
a1 = 2 .* pi .* Nnuc .* k11N .* k22N .* Lkf0;
a2 = k12N .* (ftwisted - 2 .* pi .* Nnuc .* k12N .* Lkf0);
a3 = 2 .* pi .* ((k12N).^2 - k11N .* k22N) .* dLk;
a4 = 2.* pi .* (k12N .^2 - k11N .* k22N - (k11N .* kbT ./ L_tether).*Ct);

Lk_DNA_fake = -((a1+a2-a3)./a4)
%}
%Lk_DNA_fake = (Nnuc .* Lkf0 + dLk) ./ ((Ct .* kbT)./(k22N.*L_tether)+1)

%Lk_f = dLk - Lk_DNA_fake

%Lk_DNA_zero = (Nnuc .* Lkf0 + 0) ./ ((Ct .* kbT)./(k22N.*L_tether)+1)

%% to get DNA slope (offseted to zero)

%dLk_DNA = Lk_DNA_fake - Lk_DNA_zero


Lk_DNA = ((dLk)) ./ ((Ct .* kbT)./(k22N .* L_tether.*0.34)+1)
dLk_DNA = dLk./ ((Ct .* kbT)./(k22N .* L_tether.*0.34)+1)
dLkf = (Nnuc - Ntet).* Lkf0 + (dLk - Lk_DNA)

%dLk_fiber_fake = dLk - dLk_DNA;
%dLk_DNA = 0.15 .* dLk
%dLk_f = 0.85 .* dLk

%dLk_f = dLk - Lk_DNA_fake
Lkf  = dLkf ./ (Nnuc-Ntet);

gain = 0.3;

%% help plots with transition borders

[zet_singlywrapped] = WLC_z_G(f,(CL - Nnuc.*Lwrap).* 0.34);             
[zet_extended] = WLC_z_G(f,(CL - Nnuc.*Lwrap) .* 0.34  + Nnuc .* Lextended);
[zet_unwrapped] = WLC_z_G(f,CL .* 0.34);


%plot(f,zet_singlywrapped./1000,'--','linewidth',2,'color',[0.8 0.8 1]);

plot(f,zet_extended./1000,'--','linewidth',2,'color',[0.8 0.8 1]);
hold on;
plot(f,zet_unwrapped./1000,'--','linewidth',2,'color',[0.8 0.8 1]);

%% tether length at different states (per nucleosome)

[z_wrapped, G_wrapped] = WLC_z_G_twist(f,L_tether .* 0.34,0);                        % total DNA handles (nm)

[z_fiber, G_fiber] = fiber_z_G_twist(f,kf,k22,Lkf);                                  % fiber extension and internal free energy per nucleosome [nm]

[z_singlewrap, g_singlewrap] = WLC_z_G_twist(f,(NRL - Lwrap) .* 0.34,0);             % DNA extension and internal free energy of unwrapped base pairs; per nucleosome [nm]
G_singlewrap = g_singlewrap + dG1;                                                   % correction with free energy G1 that overcomes the energy barrier

[z_extended, g_extended] = WLC_z_G_twist(f,(NRL - Lwrap) .* 0.34 + Lextended,0);     % DNA extension and internal free energy of unwrapped base pairs; per nucleosome [nm]
G_extended = g_extended +  dG1 + dG2;                                                % correction with free energy G2 that overcomes the energy barrier

[z_unwrapped, g_unwrapped] = WLC_z_G_twist(f,NRL .* 0.34,0);                         % DNA extension and internal free energy of unwrapped base pairs; per nucleosome [nm]
G_unwrapped = g_unwrapped + dG1 + dG2 + dG3;                                         % correction with free energy G3 that overcomes the energy barrier


%% calculation of all possible states that nucleosome in an array can form, its extension and free energy (!!!work is substracted -z * F !!!)
%  rows: states, columns: elements of the force ramp

D = [];
state = [];
s = 1;

for i = 1:(Nnuc - Ntet +1)
   
    for j = 1: (Nnuc - Ntet - i +3)  
        
        for k = 1: (Nnuc - i - j +3)
            
               state(s,:) = [i-1,j-1,k-1,Nnuc - (i-1) - (j-1) - (k-1)];
          
               D(s) = factorial(Nnuc)./(factorial(i-1) .* factorial(Nnuc-(i-1))) .* (factorial(Nnuc-(i-1))./(factorial(j-1).*factorial(Nnuc-(i-1)-(j-1)))) .* (factorial(Nnuc-(i-1)-(j-1))./(factorial(k-1) .*factorial(Nnuc-(i-1)-(j-1)-(k-1)))); % this particular state exists in D combinations
               
               D(s) = 1+(D(s)-1) .* degeneracy;                             % a trick to switch off the degeneracy if not needed                      
            
               D_array(s,:) = repmat(D(s),1,length(f));                     % copying the degeneracy to all columns                 
               
               z_tot(s,:) = z_wrapped + state(s,1) .* z_fiber + state(s,2) .* z_singlewrap + state(s,3) .* z_extended + state(s,4) .* z_unwrapped;          % extension of the states  multiplied by the amount of elements in this particular state

               dLk_DNA_cor = (state(s,2)+state(s,3)).* gain;
               
               L_unwrapped(s) = L_tether .* 0.34  + state(s,2) .* ((NRL - Lwrap - Lextended./0.34) .* 0.34) + state(s,3) .* Lextended;
             
               [z_unwrapped_cor, g_unwrapped_cor] = WLC_z_G_twist(f,L_unwrapped(s),dLk_DNA + dLk_DNA_cor);
               %z_tot(s,:) = z_unwrapped_cor + state(s,1) .* z_fiber + state(s,4) .* z_unwrapped; 
               G_unwrapped_2(s,:) = g_unwrapped_cor + state(s,2) .* dG1 + state(s,3) .* (dG1 + dG2); % combined 3 states in one equations as it all goes to WLC function (with corrected dLk upon unwrapping)
               
               G_tot(s,:) = ((G_unwrapped_2(s,:) + state(s,1) .* G_fiber + state(s,4) .* G_unwrapped)) - f.*z_tot(s,:)./kbT;     % internal energy + energy barriers + work done by the tether (-z * F)
               %G_tot(s,:) = ((G_wrapped + state(s,1) .* G_fiber + state(s,2) .* G_singlewrap + state(s,3) .* G_extended + state(s,4) .* G_unwrapped)) - f.*z_tot(s,:)./kbT;     % internal energy + energy barriers + work done by the tether (-z * F)
               s = s+1;
               
               
        end
        
    end
    
end

%% Boltzmann weighing - the state with the lowest total free energy dominates in total extension

for i = 1:length(f)
    
    G_total(:,i) = G_tot(:,i) - min(G_tot(:,i));                            % offseting the absolute values of energies to avoid large negative numbers that go to exponent
                 
end
    
%plot(f,G_total)
Bolz_top = z_tot .* D_array .* exp (-G_total);                              
Bolz_bot = D_array .* exp (-G_total);

[row col] = size(z_tot);

for i = 1:col
    z_final(i) = (sum(Bolz_top(:,i))./sum(Bolz_bot(:,i))./1000);
   
end

%plot(f,z_final,'linewidth',3,'color','k')
%plot(f,z_final,'linewidth',3)
%plot(G_tot(s,:),f)
%xlim([0,1.5])
hold on;
