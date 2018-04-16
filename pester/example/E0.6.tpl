ptf $
12/02/18 User:_Gib
GUI2.1      GUI_VERSION_NAME                      GUI program version number.
DLL2.1      DLL_VERSION_NAME                      DLL version number.
60000       INITIAL_COUNT                         Initial number of tumour cells
0           USE_DIVIDE_TIME_DIST                  Use divide time distribution
24          DIVIDE_TIME_1_MEDIAN                  Median (h)
1.2         DIVIDE_TIME_1_SHAPE                   Shape parameter
24          DIVIDE_TIME_2_MEDIAN                  Division time median parameter
1.2         DIVIDE_TIME_2_SHAPE                   Division time shape parameter
0           V_DEPENDENT_GROWTH_RATE               V-dependent growth rate
1           RANDOMISE_INITIAL_V                   Randomise initial cell volumes
0.2         NDAYS                                 Number of days
600         DELTA_T                               Time step
1           NT_CONC                               Number of ODE solver sub-steps.
1           VCELL_PL                              Cell volume
0.33        WELL_AREA                             Well area
0.05        MEDIUM_VOLUME                         Medium volume
0           FULLY_MIXED                           Medium is fully mixed?
1.6         VDIVIDE0                              Nominal divide volume
0.3         DVDIVIDE                              Divide volume variation
0.1         MM_THRESHOLD                          Michaelis-Menten O2 threshold
0.15        ANOXIA_THRESHOLD                      Tag threshold
3           ANOXIA_TAG_TIME                       Tag time limit
3           ANOXIA_DEATH_TIME                     Death delay time
0.15        AGLUCOSIA_THRESHOLD                   Aglucosia threshold
3           AGLUCOSIA_TAG_TIME                    Aglucosia time limit
3           AGLUCOSIA_DEATH_TIME                  Aglucosia death delay time
0           TEST_CASE                             Test case #
1234        SEED1                                 First RNG seed
5678        SEED2                                 Second RNG seed
4           NCPU                                  Number of CPUs
2           NCELLTYPES                            Number of cell types
100         CELLPERCENT_1                         Percentage of cell type 1
0           CELLPERCENT_2                         Percentage of cell type 2
1           NT_ANIMATION                          Animation interval (timesteps)
0           SHOW_PROGENY                          Show descendants of cell #
1           USE_OXYGEN                            Use Oxygen?
1           OXYGEN_GROWTH                         Oxygen growth?
1           OXYGEN_DEATH                          Anoxia death?
2e-05       OXYGEN_DIFF_COEF                      Spheroid diffusion coeff
5e-05       OXYGEN_MEDIUM_DIFF                    Medium diffusion coeff
600         OXYGEN_CELL_DIFF_IN                   Cell influx parameter Kin
600         OXYGEN_CELL_DIFF_OUT                  Cell efflux parameter Kout
0.18        OXYGEN_BDRY_CONC                      Boundary concentration
0           OXYGEN_CONSTANT                       Constant concentration
6.25e-17    OXYGEN_CONSUMPTION                    Max consumption rate
1.33        OXYGEN_MM_KM                          Michaelis-Menten Km
1           OXYGEN_HILL_N                         Hill function N
1           USE_GLUCOSE                           Use Glucose?
1           GLUCOSE_GROWTH                        Glucose growth?
1           GLUCOSE_DEATH                         Aglucosia death?
3e-07       GLUCOSE_DIFF_COEF                     Spheroid diffusion coeff
8e-06       GLUCOSE_MEDIUM_DIFF                   Medium diffusion coeff
100         GLUCOSE_CELL_DIFF_IN                  Membrane diff constant
100         GLUCOSE_CELL_DIFF_OUT                 Membrane diff constant
0.66485     GLUCOSE_BDRY_CONC                     Boundary concentration
0           GLUCOSE_CONSTANT                      Constant concentration
$ g_vmax         $ GLUCOSE_CONSUMPTION                   Max consumption rate
$ g_km           $ GLUCOSE_MM_KM                         Michaelis-Menten Km
2           GLUCOSE_HILL_N                        Hill function N
1           USE_LACTATE                           Use Lactate?
3e-07       LACTATE_DIFF_COEF                     Spheroid diffusion coeff
6e-06       LACTATE_MEDIUM_DIFF                   Medium diffusion coeff
400         LACTATE_CELL_DIFF_IN                  Membrane diff constant
400         LACTATE_CELL_DIFF_OUT                 Membrane diff constant
0.17206     LACTATE_BDRY_CONC                     Boundary concentration
3.8e-17     LACTATE_CONSUMPTION                   Max consumption rate
20          LACTATE_MM_KM                         Michaelis-Menten Km
1           LACTATE_HILL_N                        Hill function N
0           USE_TRACER                            Use Tracer?
6e-07       TRACER_DIFF_COEF                      Spheroid diffusion coeff
6e-06       TRACER_MEDIUM_DIFF                    Medium diffusion coeff
20          TRACER_CELL_DIFF_IN                   Membrane diff constant
20          TRACER_CELL_DIFF_OUT                  Membrane diff constant
1           TRACER_BDRY_CONC                      Boundary concentration
1           TRACER_CONSTANT                       Constant concentration
0           TRACER_CONSUMPTION                    Consumption rate
0           TRACER_MM_KM                          Michaelis-Menten Km
0           TRACER_HILL_N                         Hill function N
0.0738      RADIATION_ALPHA_H_1                   Alpha (hypoxia)
0.00725     RADIATION_BETA_H_1                    Beta (hypoxia)
2.5         RADIATION_OER_ALPHA_1                 OER alpha
2.5         RADIATION_OER_BETA_1                  OER beta
0.0043      RADIATION_KM_1                        Km for radiosensitivity
1           RADIATION_DEATH_PROB_1                Death prob
0           RADIATION_GROWTH_DELAY_FACTOR_1       Growth delay factor
0           RADIATION_GROWTH_DELAY_N_1            Growth delay cycles
0.0473      RADIATION_ALPHA_H_2                   Alpha (hypoxia)
0.0017      RADIATION_BETA_H_2                    Beta (hypoxia)
2.5         RADIATION_OER_ALPHA_2                 OER alpha
3           RADIATION_OER_BETA_2                  OER beta
0.0043      RADIATION_KM_2                        Km for radiosensitivity
1           RADIATION_DEATH_PROB_2                Death prob
0           RADIATION_GROWTH_DELAY_FACTOR_2       Growth delay factor
0           RADIATION_GROWTH_DELAY_N_2            Growth delay cycles
0           RADIATION_GROWTH_DELAY_ALL            Delay growth of all cells
1           USE_CELL_CYCLE                        Use cell cycle with G1, S, G2, M phases
6           T_G1_1                                G1 phase base time (h)
6           T_G1_2                                G1 phase base time (h)
8           T_S_1                                 S phase base time (h)
8           T_S_2                                 S phase base time (h)
1           T_G2_1                                G2 phase base time (h)
1           T_G2_2                                G2 phase base time (h)
0.5         T_M_1                                 M phase base time (h)
0.5         T_M_2                                 M phase base time (h)
1.5         G1_MEAN_DELAY_1                       G1 mean delay (h)
1.5         G1_MEAN_DELAY_2                       G1 mean delay (h)
1           G2_MEAN_DELAY_1                       G2 mean delay (h)
1           G2_MEAN_DELAY_2                       G2 mean delay (h)
0.1         APOPTOSIS_RATE_1                      Apoptosis rate/hr
0.1         APOPTOSIS_RATE_2                      Apoptosis rate/hr
0.6         RMR_ETA_PL_1                          PL lesion L1 creation rate
0           RMR_ETA_L_1_1                         Lesion L2b creation rate
0.14        RMR_ETA_L_2_1                         Lesion L2c creation rate
0.1         RMR_KREP_BASE_1                       Base lesion L1 repair rate
0.8         RMR_KREP_MAX_1                        Max lesion L1 repair rate
0.0278      RMR_KMIS_1_1                          Lesion misrepair rate to L2b
0.0278      RMR_KMIS_2_1                          Lesion misrepair rate to L2c
0.13        RMR_KCP_1                             Checkpoint time limit factor
1           USE_METABOLISM                        Use glucose metabolism
2           N_GA_1                                ATP moles produced per glucose mole
14          N_PA_1                                ATP moles produced per pyruvate mole
0.5         N_GI_1                                Intermediate moles produced per glucose mole
0.4         N_PI_1                                Intermediate moles produced per pyruvate mole
3           N_PO_1                                Oxygen moles consumed per pyruvate mole
140         K_H1_1                                K_H1
0.001       K_H2_1                                K_H2
0.2         K_HB_1                                K_HB
4.63e-05    K_PDK_1                               K_PDK
0.3         PDKMIN_1                              PDKmin
0.05        C_O2_NORM_1                           Nominal normal IC O2 concentration
2.5         C_G_NORM_1                            Nominal normal IC glucose concentration
1           C_L_NORM_1                            Nominal normal IC lactate concentration
0.2         ATP_S_1                               ATP production threshold for survival (fraction of peak)
0.7         ATP_G_1                               ATP production threshold for growth (fraction of peak)
1.3         ATP_RAMP_1                            Ramp factor for reducing r_G, r_P based on ATP
1.0         K_PL_1                                Pyruvate -> lactate rate constant
0.01        K_LP_1                                Lactate -> pyruvate rate constant
50          PYRUVATE_MM_KM_1                      Pyruvate Michaelis-Menten Km (uM)
2           N_GA_2                                ATP moles produced per glucose mole
14          N_PA_2                                ATP moles produced per pyruvate mole
0.6         N_GI_2                                Intermediate moles produced per glucose mole
0.6         N_PI_2                                Intermediate moles produced per pyruvate mole
3           N_PO_2                                Oxygen moles consumed per pyruvate mole
140         K_H1_2                                K_H1
0.001       K_H2_2                                K_H2
0.2         K_HB_2                                K_HB
4.63e-05    K_PDK_2                               K_PDK
0.3         PDKMIN_2                              PDKmin
0.05        C_O2_NORM_2                           Nominal normal IC O2 concentration
2.5         C_G_NORM_2                            Nominal normal IC glucose concentration
0           C_L_NORM_2                            Nominal normal IC lactate concentration
0.3         ATP_S_2                               ATP production threshold for survival (fraction of peak)
0.55        ATP_G_2                               ATP production threshold for growth (fraction of peak)
1.5         ATP_RAMP_2                            Ramp factor for reducing r_G, r_P based on ATP
0.001       K_PL_2                                Pyruvate -> lactate rate constant
0.001       K_LP_2                                Lactate -> pyruvate rate constant
20          PYRUVATE_MM_KM_2                      Pyruvate Michaelis-Menten Km
0.1         HYPOXIA_1                             Hypoxia threshold 1
1           HYPOXIA_2                             Hypoxia threshold 2
4           HYPOXIA_3                             Hypoxia threshold 3
4           HYPOXIA_THRESHOLD                     Hypoxia threshold
0.25        GROWTH_FRACTION_1                     Growth fraction threshold 1
0.1         GROWTH_FRACTION_2                     Growth fraction threshold 2
0.01        GROWTH_FRACTION_3                     Growth fraction threshold 3
1e-06       DRUG_THRESHOLD                        Drug Threshold
200         SPCRAD                                Spectral radius
0               NDRUGS_USED
0           SAVE_FACS_DATA                        Save FACS data
facs_data   SAVE_FACS_DATA_FILE_NAME              Base file name for saving FACS data
0           SAVE_FACS_DATA_TSTART                 Start time
0           SAVE_FACS_DATA_INTERVAL               Interval
1           SAVE_FACS_DATA_NUMBER                 Number
0           DUMMY_HYPOXIA_THRESHOLD               Hypoxia threshold
0           DUMMY_GROWTH_FRACTION                 Growth fraction
1           nlive                                 
0           nviable                               
0           nanoxiadead                           
0           naglucosiadead                        
0           ndrugAdead                            
0           ndrugBdead                            
0           nradiationdead                        
0           nanoxiatagged                         
0           naglucosiatagged                      
0           ndrugAtagged                          
0           ndrugBtagged                          
0           nradiationtagged                      
0           hypoxicfraction                       
0           clonohypoxicfraction                  
0           growthfraction                        
1           platingefficiency                     
1           ECoxygen                              
1           ECglucose                             
1           EClactate                             
0           ECdrugA                               
0           ECdrugAmet1                           
0           ECdrugAmet2                           
0           ECdrugB                               
0           ECdrugBmet1                           
0           ECdrugBmet2                           
1           ICoxygen                              
1           ICglucose                             
1           IClactate                             
1           ICpyruvate                            
0           ICdrugA                               
0           ICdrugAmet1                           
0           ICdrugAmet2                           
0           ICdrugB                               
0           ICdrugBmet1                           
0           ICdrugBmet2                           
1           Medglucose                            
0           Medglucose                            
1           Medlactate                            
0           MeddrugA                              
0           MeddrugAmet1                          
0           MeddrugAmet2                          
0           MeddrugB                              
0           MeddrugBmet1                          
0           MeddrugBmet2                          
0           doublingtime                          
1           Grate                                 
1           Prate                                 
1           Arate                                 
1           Irate                                 
0           f_G                                   
0           f_P                                   
0           HIF-1                                 
0           PDK1                                  
0           dividerate                            
1           MULTI                                 
0           Oxygen                                
0           Glucose                               
0           Drug_A                                
0           Drug_A_metab1                         
0           Drug_A_metab2                         
0           Drug_B                                
0           Drug_B_metab1                         
0           Drug_B_metab2                         
PROTOCOL
0
