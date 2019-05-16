function dydt = TMZ_equations(t,y,p)

% PARAMETERS
kA         = p(1);  % hr-1       - Rate constant of absorption of TMZ in gut
fu         = p(2);  % unitless   - Fraction of free TMZ
PA_csf     = p(3);  % mL/hr      - Permeability area proudct of the blood brain barrier
Qh_brain   = p(4);  % mL/hr      - Healthy brain bulk flow rate
Qt_brain   = p(5);  % mL/hr      - Brain tumor bulk flow rate
kM_T_MTIC  = p(6);  % hr-1       - Rate constant of TMZ metabolism into MTIC
kT_h_vev   = p(7);  % hr-1       - Rate constant of TMZ transport from healthy vascular to extravascular 
kT_h_evv   = p(8);  % hr-1       - Rate constant of TMZ transport from healthy extravascular to vascular
kT_t_vev   = p(9);  % hr-1       - Rate constant of TMZ transport from tumor vascular to extravascular
kT_t_evv   = p(10); % hr-1       - Rate constant of TMZ transport from tumor extravascular to vascular
kCL_T      = p(11); % hr-1       - Rate constant of TMZ clearance
kM_MTIC_C  = p(12); % hr-1       - Rate constant of MITC metabolism into methylating cation
kCL_C      = p(13); % hr-1       - Rate constant of clearance of methylating cation
kDNA       = p(14); % hr-1       - Rate constant of DNA adduct formation
Vplasma    = p(15); % mL         - Volume of plasma compartment
Vcsf       = p(16); % mL         - Volume of CSF compartment
VhBrainV   = p(17); % mL         - Volume of healthy brain vascular compartment
VhBrainEV  = p(18); % mL         - Volume of healthy brain extravascular compartment
VtBrainV   = p(19); % mL         - Volume of tumor vascular compartment (no extravascular)
rho        = p(20); % hr-1       - Growth rate
K          = p(21); % cm-3       - Carrying capacity
a2         = p(22); % mL/ug/hr   - Kill strength in resistant tumor
pm         = p(23); % unitless   - Probability of a tumor cell being marked for destruction
pd         = p(24); % unitless   - Probability of a tumor cell dying
d          = p(25); % hr-1       - Time to death
konMGMT    = p(26); % M-1hr-1    - forward rate for MGMT-adduct binding
koffMGMT   = p(27); % hr-1       - reverse rate for MGMT-adduct binding
%VtBrainEV  = VtBrainEVS + VtBrainEVR + VtBrainEVD;
%f = (pm*pd/d)*(1-exp(-y(19)*(6.02*10^15)/VtBrainEV)); % cell kill

dydt = zeros(25,1);

% VALUES OF TEMOZOLOMIDE (prodrug)
% Everything in terms of mass
% 1  = TMZ in gut
% 2  = TMZ in plasma
% 3  = TMZ in CSF
% 4  = TMZ in vascular healthy brain tissue
% 5  = TMZ in extravascular healthy brain tissue
% 6  = TMZ in vascular brain tumor tissue
% 7  = TMZ in extravascular brain tumor tissue
% 8  = TMZ cleared from plasma (virtual compartment)

dydt(1)  = -kA*y(1);
dydt(2)  =  kA*y(1) - fu*PA_csf*y(2)/Vplasma + PA_csf*y(3)/Vcsf                  - fu*Qh_brain*y(2)/Vplasma + Qh_brain*y(4)/VhBrainV                                                                   - fu*Qt_brain*y(2)/Vplasma + Qt_brain*y(6)/VtBrainV                                                  - fu*kM_T_MTIC*y(2) - kCL_T*y(2)/Vplasma; 
dydt(3)  =            fu*PA_csf*y(2)/Vplasma - PA_csf*y(3)/Vcsf - kM_T_MTIC*y(3);
dydt(4)  =                                                                         fu*Qh_brain*y(2)/Vplasma - Qh_brain*y(4)/VhBrainV - kT_h_vev*y(4) + kT_h_evv*y(5) - kM_T_MTIC*y(4);
dydt(5)  =                                                                                                                             kT_h_vev*y(4) - kT_h_evv*y(5)                  - kM_T_MTIC*y(5);
dydt(6)  =                                                                                                                                                                                               fu*Qt_brain*y(2)/Vplasma - Qt_brain*y(6)/VtBrainV - kT_t_vev*y(6) + kT_t_evv*y(7) - kM_T_MTIC*y(6);
dydt(7)  =                                                                                                                                                                                                                                                   kT_t_vev*y(6) - kT_t_evv*y(7) - kM_T_MTIC*y(7);                                                                                                                              
dydt(8)  =                                                                                                                                                                                                                                                                                                                        kCL_T*y(2)/Vplasma;


% VALUES OF MTIC (Drug)
% Everything in terms of mass
% 9  = MITC in plasma
% 10 = MITC in CSF
% 11 = MITC in vascular healthy brain tissue
% 12 = MITC in extravascular healthy brain tissue
% 13 = MITC in vascular brain tumor tissue
% 14 = MITC in extravascular brian tumor tissue

dydt(9)  =                                                                                                                                                                                                                                                                                                  + fu*kM_T_MTIC*y(2)                                                                                         - Qh_brain*y(9)/Vplasma + Qh_brain*y(11)/VhBrainV - Qt_brain*y(9)/Vplasma + Qt_brain*y(13)/VtBrainV;
dydt(10) =                                                        kM_T_MTIC*y(3);   
dydt(11) =                                                                                                                                                            kM_T_MTIC*y(5)                                                                                                                                                                                                                                    + Qh_brain*y(9)/Vplasma - Qh_brain*y(11)/VhBrainV;
dydt(12) =                                                                                                                                                                              kM_T_MTIC*y(4)                                                                                                                                               - kM_MTIC_C*y(12);                                                  
dydt(13) =                                                                                                                                                                                                                                                                                                                                                                                kM_T_MTIC*y(6)                                                                  + Qt_brain*y(9)/Vplasma - Qt_brain*y(13)/VtBrainV;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
dydt(14) =                                                                                                                                                                                                                                                                                                                                                             - kM_MTIC_C*y(14)                 + kM_T_MTIC*y(7);


% VALUES OF METHYLATING CATION
% Everything in terms of mass
% 15 = Methylating Cation in extravascular healthy brain tissue
% 16 = Methylating Cation in extravascular brain tumor tissue
% 17 = Methylating Cation cleared (virtual compartment)

dydt(15) =                                                                                                                                                                                                                                                                                                                                             kM_MTIC_C*y(12)                                                                                                                                                       - kCL_C*y(15)               - kDNA*y(15);
dydt(16) =                                                                                                                                                                                                                                                                                                                                                              kM_MTIC_C*y(14)                                                                                                                                                    - kCL_C*y(16)              - kDNA*y(16);
dydt(17) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     kCL_C*y(15) + kCL_C*y(16);


% VALUES OF DNA ADDUCT
% Everything in terms of mass
% 18 = DNA Adduct in extravascular healthy brain tissue
% 19 = DNA Adduct in extravascular brain tumor tissue
% 23 = DNA adduct in dead tissue
% 24 = repaired DNA adduct 
% 25 = MGMT
dydt(18) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 kDNA*y(15);
dydt(19) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              kDNA*y(16) - (rho*(y(21)+y(20))*(y(20)+y(21))/K+ (pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*y(20)+(pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*a2*y(21))*y(19)/(y(20)+y(21)) - konMGMT*(y(25)/(y(20)+y(21)))*(y(19)/(y(20)+y(21)))+ koffMGMT*y(24)/(y(20)+y(21));
dydt(23) = (rho*(y(21)+y(20))*(y(20)+y(21))/K+ (pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*y(20)+(pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*a2*y(21))*y(19)/(y(20)+y(21));
dydt(24) = konMGMT*(y(25)/(y(20)+y(21)))*(y(19)/(y(20)+y(21))) - koffMGMT*y(24)/(y(20)+y(21));
dydt(25) = -konMGMT*(y(25)/(y(20)+y(21)))*(y(19)/(y(20)+y(21))) + koffMGMT*y(24)/(y(20)+y(21));
% VALUES OF TUMOR
% Everything below is in terms of volume
% 20 = Volume of tumor that is susceptible to TMZ
% 21 = Volume of tumor that is resistant to TMZ
% 22 = Volume of tumor that is dead (virtual compartment)
dydt(20) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           rho*y(20)*(1-(y(20)+y(21))/K)                                  - (pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*y(20);
dydt(21) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          rho*y(21)*(1-(y(20)+y(21))/K)                                                                 - (pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*a2*y(21);
dydt(22) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          rho*(y(21)+y(20))*(y(20)+y(21))/K         + (pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*y(20)  + (pm*pd/d)*(1-exp(-y(19)*(6.02*10^12)/(y(20)+y(21))))*a2*y(21);
