function [z_final] = DNA_fiber_twisted(F,L0,dLk,length_correction,k11,k22,k12,Lkf0,Nnuc,Ntet)  

% fitting rotation curves of chromatin fibers at contant force


% define: force at which the twist experiment was performed (F), 
        % length of the tether when the fiber is folded (nm) (L0),
        % twist ramp (for example from -50:0.1:50) (dLk), 
        % length correction (due to off-center bead attachment),
        % stretch stiffness (k11), twist modulus (k22) (pN *nm^2), twist-stretch coupling (k12)
        % number of nucleosomes (Nnuc) and tetrasomes in the assembled fiber

% Artur Kaczmarczyk, akaczmarczyk88@gmail.com, August 2018


kbT = 4.1;                                              % kBT - Boltzmann constant in room temperature (pN*nm);                                                          
C = 100;                                                % DNA twist persistence length (nm)
Cp = 12;                                                % twist persistence length of plectonemic DNA in the PRESENCE of nucleosomes (nm) (for bare DNA: 24 nm)
At = 50;                                                % DNA stretch persistence length (nm)
Ct = C * (1 - C ./ (4 .* At) .* sqrt(kbT ./ (At.*F)));  % twist persistence length as a function of force with constant C = 100 nm


S = 1000;
g = (F - sqrt(F*kbT/At));                               % free energy per nm of torsionally unconstrained DNA;
tau = Ct .* (2 .* pi .* dLk ./ L0) .* kbT;              % DNA restoring torque 
tau_buck = sqrt(2 .* kbT .* Cp .* g ./ (1 - Cp/Ct));    % DNA buckling torque
tau_melt = -11;                                         % DNA melting torque
Lk_buck = tau_buck .* L0 /(2 .*pi .* kbT .* Ct);        % maximal linking number that can be absorbed by DNA before buckling
Lk_melt = tau_melt .* L0 /(2 .*pi .* kbT .* Ct);        % maximal linking number that can be absorbed by DNA before melting


%% DNA extension as a function of twist

% extension of twisted DNA    (z/L)
z_t =  (1 - 0.5 .* sqrt(kbT/(At .* F))) - (C.^2 ./ 16) .* (2 .* pi .* dLk ./ L0).^2 .* ( kbT ./ (At*F)) .^ (3/2);   

% decrease in DNA extension due to plectoneme formation (nm/turn)
z_slope = (2 .* pi .* (1 - 0.5 .* sqrt(kbT/(At .* F)) - C.^2/(16 .* Ct.^2) * (kbT ./(At.*F)).^(3/2).*(tau_buck./kbT).^2))/((tau_buck./kbT).*(1/Cp - 1/Ct));  

% free energy of twisted DNA    (Gt / L0)
G_t = - g + (0.5 * kbT .* Ct) .* (2 .* pi .* dLk ./ L0).^2;                 

% free energy of plectonemic DNA gained per each turn;
dG_DNA = 2 .* pi .* tau_buck;                                               


%% conditions for DNA buckling at positive twist (omitted if it's a chromatin sample)

%{
for i = 1:length(dLk)
    if dLk(i) >= Lk_buck 
        z_t(i) = z_t(i-1);                                                  % once buckling torque is reached, the change in extension is attributed only to the plectoneme formation
        z(i) = (z_t(i).* L0  - (abs(dLk(i)-Lk_buck)) .* z_slope) ./ L0;     % keep in mind that the extension of twisted DNA is expressed in z/L
    
    else  
        z(i) = z_t(i);
    end
    
end

%% DNA buckling and melting at negative twist
 
k = 1;
start = dLk(k);

while start ~= 0                                                                
    k = k+1;
    start = dLk(k);
    k;                                                                      % the position of point in the turns array that starts to go up from negative to positive twist
end
                                                                            
for i = 1:(k-2)
    
  if -tau_buck > tau_melt                                                   % DNA buckles when melting torque is not reached
      
      if dLk(k-i) <= - Lk_buck 
            z(k-i-1) = z(k-i);   
            z(k-i) = (z(k-i).* L0  - (abs(dLk(k-i))-Lk_buck) .* z_slope) ./ L0; 
            
      else
            z(k-i) = z(k-i);
                
      end
      
  else
      
      if tau(k-i) <= tau_melt  
            z(k-i) = z(k-i+1);
      else
            z(k-i) = z(k-i);
      end
      
  end
  
end
 z(1) = z(2)                                                                % correction of the first data point as the end point in this particular loop 
                                                                 
%%% erasing negative extensions (function goes to negative infinity)

for i = 1:length(dLk)  
    
    if z(i) < length_correction .* 1000;                                     % adding some offset 
        z(i) = length_correction .* 1000;
    else
        z(i) = z(i);
    end
end

%}
%% chromatin fiber properties 

NRL = 167;                                                                            
%Nnuc = 26;                                                                 % number of nucleosomes
%N601 = Nnuc                                                                % for plotting purposes
z0_f = 1.7;                                                                 % length of the fiber in the absence of force (nm/nuc)

%k11 = 0.32;                                                                % stretch spring constant per nucleosome (pN/nm);
k11N = k11 ./ ((Nnuc).*z0_f);                                               % per whole fiber
%k22 = 2;                                                                   % twist modulus per nucleosome(pN*nm);
k22N = k22 ./ ((Nnuc).*z0_f);                                               % per whole fiber
%k12 = 0.1;                                                                 % twist-stretch coupling factor (pN/rad);
k12N = k12 ./ ((Nnuc).*z0_f);                                              % per whole fiber
%Lkf0 = -0.15;                                                              % linking number of a fiber, per nucleosome;


%% Combination of DNA and fiber responce to twisting - distribution of the applied twist
 
% L0 is the measured tether length at max. 2 pN (when fiber is not extended, serves to obtain the length of DNA handles) 

% L0 = L0 + (N601 - Nnuc) .* NRL .* 0.34;                                    % to correct for missing nucleosomes (for plotting purposes)

% distribution of dLk between DNA and fiber (extended formula that includes twist-stretch coupling from Meng's thesis)

%a1 = 2 .* pi .* Nnuc .* k11N .* k22N .* Lkf0;
%a2 = k12N .* (F - 2 .* pi .* Nnuc .* k12N .* Lkf0);
%a3 = 2 .* pi .* ((k12N).^2 - k11N .* k22N) .* dLk;
%a4 = 2.* pi .* (k12N .^2 - k11N .* k22N - (k11N .* kbT ./ L0tether).*Ct);
%Lk_DNA = -((a1+a2-a3)./a4);

% FORMULA DERIVED BY COMPARING TORQUE IN THE DNA AND IN THE FIBER 

Lk_DNA = dLk ./ ((Ct .* kbT) ./ (k22N .* L0) + 1);        

Lk_f = dLk - Lk_DNA;

%Lk_f = Nnuc .* Lkf0 + dLk - Lk_DNA;     % added offset: 'Nnuc.*Lkf0' to account for fiber chirality (not necessary at low forces when no unwrapping occurs)

%sigma = dLk/(L0/0.34/10.4);


%% Characterization of DNA chromatin fiber under twist


         %%z_t =  (1 - 0.5 .* sqrt(kbT/(At .* F))) - (C.^2 ./ 16) .* (2 .* pi .* dLk ./ L0).^2 .* ( kbT ./ (At*F)) .^ (3/2);      % if the total twist would go to DNA (in the absence of chromatin fiber)  
z_handle =  (L0) .* (1 - 0.5 .* sqrt(kbT/(At .* F))) - (C.^2 ./ 16) .* (2 .* pi .* Lk_DNA ./ L0).^2 .* ( kbT ./ (At*F)) .^ (3/2);   % extension of twisted DNA (size of the fiber is neglected here), small offset to plot the fit in the middle of the data

% change in extension of a fiber  due to twist-stretch coupling (includes all nucleosomes)
dZ_f = (F - k12N .* 2 .* pi .* Lk_f) ./ k11N;                          

% fiber dimension
Z_f = (Nnuc - Ntet)  .* z0_f + dZ_f;                                                 

% torque in a fiber
tau_f = k12N .* dZ_f  + k22N .* 2 .* pi .* Lk_f;                            

% torque in DNA handles
tau_handle = Ct .* (2 .* pi .* Lk_DNA ./ L0) .* kbT;                        

% free energy of a fiber
G_f = 0.5 .* k11 .* (dZ_f).^2 + 0.5 .* k22 .* (2 .* pi .* Lk_f ).^2 + k12 .* dZ_f .* (2 .* pi .* Lk_f);   

% free energy of DNA gained per each 1 turn after buckling;
dG_DNA = 2 .* pi .* tau_buck;                                               %


%% Constant torque when buckling or melting torque is reached

tau_total = [];

for i=1:length(tau_f)
    if tau_f(i) > tau_buck
        tau_total(i) = tau_buck;
    elseif  tau_f(i) < -11
        tau_total(i) = -11;
    else
        tau_total(i) = tau_f(i);
    end
    
end

%plot(dLk,tau_total)

%% Plotting DNA and fiber together

for i = 1:length(dLk)
    
    z_FIB(i) = z_handle(i) + Z_f(i);
    
    if Lk_DNA(i) >= Lk_buck && i>1
       
        z_FIB(i) = z_FIB(i-1);                                               % once buckling torque is reached, the change in extension is attributed only to the plectoneme formation
        z_FIB(i) = z_FIB(i) - ((dLk(i)-dLk(i-1)) .* z_slope);                % shortening of the tether lengh due to plectoneme formation (all the twist goes to DNA now)
    else  
        z_FIB(i) = z_FIB(i);
    end
    
end

%% DNA buckling and melting at negative twist

k = 1;
start = dLk(k);

while start ~= 0                                                                
    k = k+1;
    start = dLk(k);
    k;                                                                      % the position of point in the turns array that starts to go up from negative to positive twist
end
                                                                            

for i = 1:(k-2)
    
  if -tau_buck > tau_melt                                                   % DNA buckles when melting torque is not reached
      
      if dLk(k-i) <= - Lk_buck 
          
            z_FIB(k-i-1) = z_FIB(k-1);
            z_FIB(k-i-1) = z_FIB(k-i) - (abs(dLk(k-i))-abs(dLk(k-i+1))) .* z_slope; 
            
      else
            z_FIB(k-i) = z_FIB(k-i);
                
      end
      
  else
      
      if tau(k-i) <= tau_melt  
            z_FIB(k-i) = z_FIB(k-i+1);
      else
            z_FIB(k-i) = z_FIB(k-i);
      end
      
  end
  
end
 z_FIB(1) = z_FIB(2);                                                       % correction of the first data point which is the end point in this particular for-loop 

%% Erasing negative extensions

for i = 1:length(dLk)  
    
    if z_FIB(i) < (length_correction .* 1000) + 100;
        z_FIB(i) = (length_correction .* 1000) + 100;                       % manually corrected offset due to the off-center attachment (trace does not goes to 0)
    else
        z_FIB(i) = z_FIB(i);
    end
end


%% Boltzmann weighing for twisted fiber above 3 pN (unwrapping needs to be included)


if F > 3;                                                                   % apply Boltzman only when unwrapping occurs (above 3 pN at positive twist)
    
    % unstacking energy per nucleosome
    Gu = 4.1 .* 18;                                                         
   
    % varying the number of folded nucleosomes
    %Nnuc = Nnuc -2 % manual offset because the twisting at 4.5 pN is done already at the level when some nucs might have unwrapped
    for i = 1 : (Nnuc - Ntet +1)   
    
        N(i) = Nnuc - Ntet - (i-1);                                             % number of folded nucleosomes 
    
        z_unfolded(i) = (i-1) .* (((NRL - 147) + 58) .* 0.34);                  % gain in extension upon when one nucleosome unfolds
    
        Z_DNA(i,:) = z_handle + z_unfolded(i);                                  % tether length increases when fiber unfolds 
        % stretch spring constant (pN/nm) (k11 is per nucleosome)
        k11N(i) = k11 ./ N(i) .* z0_f;                                          % per fiber
        % twist modulus (pN*nm) (k22 is per nucleosome);
        k22N(i) = k22 ./ N(i) .* z0_f;                                          % per fiber
        % twist-stretch coupling factor (pN/rad);
        k12N(i) = k12 ./ N(i) .* z0_f;                                          % per fiber
      
        % DISTRIBUTION OF TWIST BETWEEN DNA AND FIBER
        
        dLk_DNA(i,:) =  (dLk) ./ ((Ct .* kbT)./(k22N(1) .* Z_DNA(i))+1); 
       
        dLk_f(i,:) = (Nnuc-Ntet).* Lkf0  +  dLk - dLk_DNA(i,:);                 % HERE IT IS IMPORTANT TO OFFSET LK_F DUE TO FIBER'S CHIRALITY 
   
        
        % free energy DNA with added F*z_DNA (because work will be substracted at the end in the Boltzmann formula)
        G_DNAtether(i,:) = F .* (Z_DNA(i,:)) -  Z_DNA(i,:) .* F .* (1-sqrt(kbT./(F.*At))+F./(2.*S)) + Z_DNA(i,:) .* ((0.5 * kbT .* Ct) .* (2 .* pi .* dLk_DNA(i,:) ./ Z_DNA(i,:)).^2);
                                         %G_bareDNA = - L0 .* F .* (1-sqrt(kbT./(F.*At))+F./(2.*S)) +     L0     .* ((0.5 * kbT .* Ct) .* (2 .* pi .* dLk ./ L0).^2);      
                              
        % fiber dimension
        dZ_fiber(i,:) = (F - (k12N(i)) .* 2 .* pi .* dLk_f(i)) ./ k11N(i);
        Z_fiber(i,:) = N(i) .* z0_f + dZ_fiber(i,:);                  
        
        % free energy of the fiber with different number of nucleosomes
        G_fiber(i,:) = ((F .^2 ./ (2 .* k11N(i))) + 0.5 .* k22N(i) .* (2 .* pi .* dLk_f(i,:)).^2);  
   
        % total tether length
        Z_total(i,:) = Z_DNA(i,:) + Z_fiber(i);
        
        % tau_fiber(i,:) = k12 ./ ((Nnuc - (i-1)).*z0_f) .* dZ_fiber(i,:)  + k22 ./ ((Nnuc - (i-1)).*z0_f) .* 2 .* pi .* dLk_f2(i,:); 
        
        
        
        % degeneracy
        D(i) = factorial(Nnuc) ./ (factorial(i-1) .* factorial(Nnuc-i+1));         
        D(i) = 1+(D(i)-1) .* 1; 
        D_array(i,:) = repmat(D(i),1,length(dLk)); 
    
        G_all (i,:) = (G_DNAtether(i,:) + G_fiber(i,:) +  ((i-1) .* Gu) - F .* Z_total(i))./kbT;
    
   
    end

    
    % offseting the absolute values of energies to avoid large negative numbers that go to exponent; they serve only for Boltmann weighing so it's fine ;]
    
    for i = 1:length(dLk)
        
        G_all_offset(:,i) = G_all(:,i) - min(G_all(:,i));     
        
    end

    
    
    % Boltzmann
    
    for i = 1: Nnuc - Ntet - 1
        
        Z(i,:) = D_array(i,:) .* exp(-(G_all_offset(i,:)));
        zZ(i,:) = Z_total(i,:) .* Z(i,:);
    
    end

    
    [row, col] = size(zZ);

    
    for i = 1:col
        
      z_final(i) = sum(zZ(:,i))./sum(Z(:,i));
      
    end

    
plot(dLk,z_final,'linewidth',3);
box on;
hold on;

else
    
        % analytical formulae for low force when there's no unfolding
        
        z_final = z_FIB;                                                    
        plot(dLk,z_final,'linewidth',3);
        box on;
        hold on;
        
end


