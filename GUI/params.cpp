#include <qstring.h>
#include "params.h"
#include "QDebug"

Params::Params()
{
    static infoStruct label_info[] = {
        {"PARENT_0", "Diffusion coefficient within the blob."},
        {"PARENT_1", "Diffusion coefficient in the medium."},
        {"PARENT_2", "Cell influx parameter Kin.  The rate of mass transport into the cell is Kin.Cex - Kout.Cin (currently no dependence on cell surface area)."},
        {"PARENT_3", "Cell efflux parameter Kout.  The rate of mass transport into the cell is Kin.Cex - Kout.Cin (currently no dependence on cell surface area)."},
        {"PARENT_4", "Half-life of the compound, used to calculate the decay rate.  This is the same in the cell and in the medium."},
        {"PARENT_CT1_0", "Kmet0 is the maximum rate of metabolism.  The actual rate is the product of drug concentration Cdrug, Kmet0 and a sigmoid function of O2 concentration C_O2, with parameters C2 and KO2:\n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_1", "C2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate.  When C_O2 = 0, the function = 1, when C_O2 >> KO2, the function = 1 - C2: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_2",  "KO2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_3", "Vmax and Km are parameters that determine the dependence of the maximum rate of metabolism on drug concentration: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug))"},
         {"PARENT_CT1_4", "Vmax and Km are parameters that determine the dependence of the maximum rate of metabolism on drug concentration: \n\
metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug))"},
        {"PARENT_CT1_5", "Klesion is currently unused."},
        {"PARENT_CT1_6", "The O2 concentration in the kill experiment."},
        {"PARENT_CT1_7", "The drug concentration in the kill experiment."},
        {"PARENT_CT1_8", "The duration the kill experiment."},
        {"PARENT_CT1_9", "The kill fraction achieved in the kill experiment (1 - SF)."},
        {"PARENT_CT1_10", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)"},
        {"PARENT_CT1_11", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)"},
        {"PARENT_CT1_12", "Sensitisation of the cells to radiation is determined by three parameters.  The usual radiation kill parameters OER_alpha and OER_beta are multiplied by the sensitisation enhancement ratio SER: \n\
         SER = (C_O2 + SER_KO2*(Cdrug*SER_max + SER_Km)/(Cdrug + SER_Km))/(C_O2 + SER_KO2)"},
        {"PARENT_CT1_13",  "n_O2 is one of the three parameters of the basic sigmoid function of C_O2 that determines metabolism rate: \n\
 metabolism rate = dMdt = Cdrug.(1 - C2 + C2.KO2^n_O2/(KO2^n_O2 + C_O2^n_O2)).Kmet0   If Vmax > 0, Kmet0 is replaced by (Kmet0 + Vmax*Cdrug/(Km + Cdrug)) \n\
 This the rate of transformation of parent drug to metabolite 1, or of metabolite 1 to metabolite 2, or removal of metabolite 2"},
        {"PARENT_CT1_14", "The death probability of a drug-tagged cell at time of division."},
        {"PARENT_CT1_15", "The kill probability rate parameter.\n\
 E.g. for kill model 1, kill probability rate r = Kd*dM/dt = Kd*kmet*Cdrug, and for duration t, SF = exp(-rt)"},
        {"PARENT_CT1_16", "This box is ticked if the drug is cytotoxic and kill parameters are provided."},
        {"PARENT_CT1_17", "Using dMdt = Cdrug*(1 - C2 + C2*KO2/(KO2 + C_O2))*Kmet0, the kill probability Pkill in time dt for each model is: \n\
1. Kd*dMdt*dt  2. Kd*Cdrug*dMdt*dt  3. Kd*dMdt^2*dt  4. Kd*Cdrug*dt  5. Kd*Cdrug^2*dt"},
        {"PARENT_CT1_18", "This box is ticked if the drug sensitises the cells to radiation."},
    };

    PARAM_SET params[] = {

{"GUI_VERSION_NAME", 0, 0, 0,
 "GUI0.00",
 "GUI program version number."},

{"DLL_VERSION_NAME", 0, 0, 0,
 "DLL0.00",
 "DLL version number."},

{"INITIAL_COUNT", 1000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"USE_DIVIDE_TIME_DIST", 1, 0, 1,
"Use divide time distribution",
"The divide time will be a random variate from a log-normal distribution."},

{"DIVIDE_TIME_1_MEDIAN", 18, 0, 0,
"Median (h)",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_1_SHAPE", 1.1, 0, 0,
"Shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE_TIME_2_MEDIAN", 18, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_2_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"V_DEPENDENT_GROWTH_RATE", 0, 0, 1,
"V-dependent growth rate",
"The growth rate of a cell is proportional to the volume."},

{"RANDOMISE_INITIAL_V", 1, 0, 1,
"Randomise initial cell volumes",
"The volumes of the initial cell population are randomised."},

{"NDAYS", 20.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 600, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[mins]"},

{"NT_CONC", 1, 0, 0,
"Number of ODE solver sub-steps.",
"The number of subdivisions of the major time step, for the ODE diffusion-reaction solver."},

{"VCELL_PL", 1.0, 0, 0,
"Cell volume",
"Nominal cell volume."},

{"WELL_AREA", 0.33, 0, 0,
"Well area",
"Cross-sectional area of the well."},

{"MEDIUM_VOLUME", 0.2, 0, 0,
"Medium volume",
"Volume of the medium in which the spheroid is growing."},

{"FULLY_MIXED", 0, 0, 1,
"Medium is fully mixed?",
"The medium is fully mixed"},

//{"UNSTIRRED_LAYER", 0.01, 0, 0,
//"Unstirred layer width",
//"Thickness of the unstirred layer around the spheroid (cm)."},

{"VDIVIDE0", 1.6, 0, 0,
"Nominal divide volume",
"Nominal multiple of normal cell volume at which division occurs."},

{"DVDIVIDE", 0.3, 0, 0,
"Divide volume variation",
"Variation (+/-) about nominal divide volume multiple."},

{"MM_THRESHOLD", 0.1, 0, 0,
"Michaelis-Menten O2 threshold",
"O2 concentration at which the 'soft-landing' adjustment to the Michaelis-Menten function kicks in.\n\
[uM]"},

{"ANOXIA_THRESHOLD", 0.15, 0, 0,
"Tag threshold",
"A cell begins to experience starvation of oxygen (anoxia) or glucose (aglucosia) leading to cell death at the oxygen/glucose concentration given by this threshold value."},

{"ANOXIA_TAG_TIME", 3.0, 0, 0,
"Tag time limit",
"Length of time under anoxia (O2 < anoxia threshold) or aglucosia (glucose < aglucosia threshold) after which a cell is tagged to die of anoxia or aglucosia."},

{"ANOXIA_DEATH_TIME", 3.0, 0, 0,
"Death delay time",
"Time taken for a cell to die after it is tagged to die of anoxia or aglucosia."},

{"AGLUCOSIA_THRESHOLD", 0.15, 0, 0,
"Aglucosia threshold",
"A cell begins to experience aglucosia leading to cell death at the glucose concentration given by this threshold value."},

{"AGLUCOSIA_TAG_TIME", 3.0, 0, 0,
"Aglucosia time limit",
"Length of time under aglucosia (glucose < aglucosia threshold) after which a cell is tagged to die of aglucosia.\n\
[h]"},

{"AGLUCOSIA_DEATH_TIME", 3.0, 0, 0,
"Aglucosia death delay time",
"Time taken for a cell to die after it is tagged to die of aglucosia.\n\
[h]"},

{"TEST_CASE", 0, 0, 0,
"Test case #",
"Number of the test case to run.  The default value of 0 is for a normal run"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NCELLTYPES", 2, 0, 0,
"Number of cell types",
"Maximum number of cell types in the spheroid.  The initial percentage of each type must be specified"},

{"CELLPERCENT_1", 100, 0, 100,
"Percentage of cell type 1",
"Percentage of cell type 1"},

{"CELLPERCENT_2", 0, 0, 100,
"Percentage of cell type 2",
"Percentage of cell type 2"},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"SHOW_PROGENY", 0, 0, 0,
 "Show descendants of cell #",
 "All the descendants of cell with the specified ID are highlighted.  (0 = no selection)"},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

{"OXYGEN_GROWTH", 1, 0, 1,
"Oxygen growth?",
"The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose"},

{"OXYGEN_DEATH", 1, 0, 1,
"Anoxia death?",
"Oxygen controls death by anoxia"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "Spheroid diffusion coeff",
 "Constituent diffusion coefficient in the spheroid"},

{"OXYGEN_MEDIUM_DIFF", 5.0e-5, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"OXYGEN_CELL_DIFF_IN", 600, 0, 0,
 "Cell influx parameter Kin",
 "Cell membrane diffusion constant Kin"},

{"OXYGEN_CELL_DIFF_OUT", 600, 0, 0,
 "Cell efflux parameter Kout",
 "Cell membrane diffusion constant Kout"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"OXYGEN_CONSUMPTION", 6.25e-17, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

{"OXYGEN_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"OXYGEN_HILL_N", 1, 1, 2,
 "Hill function N",
 "Oxygen uptake rate Hill function N"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

{"GLUCOSE_GROWTH", 1, 0, 1,
"Glucose growth?",
"The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose"},

{"GLUCOSE_DEATH", 1, 0, 1,
"Aglucosia death?",
"Glucose controls death by aglucosia"},

{"GLUCOSE_DIFF_COEF", 3.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_MEDIUM_DIFF", 8.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"GLUCOSE_CELL_DIFF_IN", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kin"},

{"GLUCOSE_CELL_DIFF_OUT", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kout"},

{"GLUCOSE_BDRY_CONC", 5.5, 0, 0,
 "Boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"GLUCOSE_CONSUMPTION", 7.1e-17, 0, 0,
 "Max consumption rate",
 "GLUCOSE consumption rate"},

{"GLUCOSE_MM_KM", 220, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"GLUCOSE_HILL_N", 1, 1, 2,
 "Hill function N",
 "Glucose uptake rate Hill function N"},

{"USE_LACTATE", 1, 0, 1,
"Use Lactate?",
"Lactate is simulated"},

//    {"LACTATE_GROWTH", 1, 0, 1,
//    "Lactate growth?",
//    "The rate of growth of a cell is the maximum rate multiplied by the fractional rates of metabolism of both O2 and glucose, or something"},

    {"LACTATE_DIFF_COEF", 3.0e-7, 0, 0,
     "Spheroid diffusion coeff",
     "LACTATE diffusion coefficient"},

    {"LACTATE_MEDIUM_DIFF", 6.0e-6, 0, 0,
     "Medium diffusion coeff",
     "Constituent diffusion coefficient in the medium"},

    {"LACTATE_CELL_DIFF_IN", 400, 0, 0,
     "Membrane diff constant",
     "Cell membrane diffusion coefficient Kin"},

    {"LACTATE_CELL_DIFF_OUT", 400, 0, 0,
     "Membrane diff constant",
     "Cell membrane diffusion coefficient Kout"},

    {"LACTATE_BDRY_CONC", 3, 0, 0,
     "Boundary concentration",
     "LACTATE boundary concentration"},

    {"LACTATE_CONSUMPTION", 3.8e-17, 0, 0,
     "Max consumption rate",
     "LACTATE consumption rate"},

    {"LACTATE_MM_KM", 20, 0, 0,
     "Michaelis-Menten Km",
     "Michaelis-Menten Km (uM)"},

    {"LACTATE_HILL_N", 1, 1, 2,
     "Hill function N",
     "Lactate uptake rate Hill function N"},

{"USE_TRACER", 0, 0, 1,
"Use Tracer?",
"Tracer is simulated"},

{"TRACER_DIFF_COEF", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "TRACER diffusion coefficient"},

{"TRACER_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"TRACER_CELL_DIFF_IN", 20, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kin"},

{"TRACER_CELL_DIFF_OUT", 20, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kout"},

{"TRACER_BDRY_CONC", 1.0, 0, 0,
 "Boundary concentration",
 "TRACER boundary concentration"},

{"TRACER_CONSTANT", 1, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"TRACER_CONSUMPTION", 0, 0, 0,
 "Consumption rate",
 "TRACER consumption rate"},

{"TRACER_MM_KM", 0, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"TRACER_HILL_N", 0, 0, 2,
 "Hill function N",
 "Tracer uptake rate Hill function N"},

//==========================
// Radiotherapy parameters
//==========================

{"RADIATION_ALPHA_H_1", 0.0738, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_1", 0.00725, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_1", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_1", 2.5, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_1", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_DEATH_PROB_1", 1.0, 0, 0,
"Death prob",
"Probability of death at mitosis for a cell tagged for damage by radiation"},

{"RADIATION_GROWTH_DELAY_FACTOR_1", 0.0, 0, 0,
"Growth delay factor",
"For a damaged cell, cell cycle is delayed for a number of hours given by this factor x radiation dose"},

{"RADIATION_GROWTH_DELAY_N_1", 0, 0, 0,
"Growth delay cycles",
"For a damaged cell, cell cycle delay persists for a number of cell cycles"},

{"RADIATION_ALPHA_H_2", 0.0473, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_2", 0.0017, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_2", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_2", 3.0, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_2", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_DEATH_PROB_2", 1.0, 0, 0,
"Death prob",
"Probability of death at mitosis for a cell tagged for damage by radiation"},

{"RADIATION_GROWTH_DELAY_FACTOR_2", 0.0, 0, 0,
"Growth delay factor",
"For a damaged cell, cell cycle is delayed for a number of hours given by this factor x radiation dose"},

{"RADIATION_GROWTH_DELAY_N_2", 0, 0, 0,
"Growth delay cycles",
"For a damaged cell, cell cycle delay persists for a number of cell cycles"},

{"RADIATION_GROWTH_DELAY_ALL", 0, 0, 0,
"Delay growth of all cells",
"Cell cycle delay is also applied to cells that are not fated to die"},

    {"USE_CELL_CYCLE", 1,0,1,
     "Use cell cycle with G1, S, G2, M phases",
     "Cell cycle parameters determine the time spent in each phase.\n\
     In the case of G1 and G2, an exponentially distributed random delay is added"},

     {"T_G1_1", 6, 0, 0,
     "G1 phase base time (h)",
     "Deterministic component of time spent in phase G1"},

     {"T_G1_2", 6, 0, 0,
     "G1 phase base time (h)",
     "Deterministic component of time spent in phase G1"},

     {"T_S_1", 8, 0, 0,
     "S phase base time (h)",
     "Time spent in phase S"},

     {"T_S_2", 8, 0, 0,
     "S phase base time (h)",
     "Time spent in phase S"},

     {"T_G2_1", 1, 0, 0,
     "G2 phase base time (h)",
     "Deterministic component of time spent in phase G2"},

     {"T_G2_2", 1, 0, 0,
     "G2 phase base time (h)",
     "Deterministic component of time spent in phase G2"},

     {"T_M_1", 0.5, 0, 0,
     "M phase base time (h)",
     "Time spent in phase M"},

     {"T_M_2", 0.5, 0, 0,
     "M phase base time (h)",
     "Time spent in phase M"},

     {"G1_MEAN_DELAY_1", 1.5, 0, 0,
     "G1 mean delay (h)",
     "Mean of the random component of time spent in phase G1 (exponentially distributed)"},

     {"G1_MEAN_DELAY_2", 1.5, 0, 0,
     "G1 mean delay (h)",
     "Mean of the random component of time spent in phase G1 (exponentially distributed)"},

     {"G2_MEAN_DELAY_1", 1, 0, 0,
     "G2 mean delay (h)",
     "Mean of the random component of time spent in phase G2 (exponentially distributed)"},

     {"G2_MEAN_DELAY_2", 1, 0, 0,
     "G2 mean delay (h)",
     "Mean of the random component of time spent in phase G2 (exponentially distributed)"},

     {"APOPTOSIS_RATE_1", 0.1, 0, 0,
     "Apoptosis rate/hr",
     "The rate of passage from state = DYING to state = DEAD, probability/hour"},

     {"APOPTOSIS_RATE_2", 0.1, 0, 0,
     "Apoptosis rate/hr",
     "The rate of passage from state = DYING to state = DEAD, probability/hour"},


     {"RMR_ETA_PL_1", 0.6, 0, 0,
     "PL lesion L1 creation rate",
     "Coefficient of rate of creation of L1 potentially lethal lesions: eta_PL"},

     {"RMR_ETA_L_1_1", 0.0, 0, 0,
     "Lesion L2b creation rate",
     "Coefficient of rate of creation of L2b lesions: eta_L(1)"},

     {"RMR_ETA_L_2_1", 0.14, 0, 0,
     "Lesion L2c creation rate",
     "Coefficient of rate of creation of L2c lesions: eta_L(2)"},

     {"RMR_KREP_BASE_1", 0.1, 0, 0,
     "Base lesion L1 repair rate",
     "Base coefficient of rate of repair of L1 potentially lethal lesions: Krepair_base (Curtis's epsilon_PL)\n\
     The true repair rate varies linearly from Krepair_base to Krepair_max over the S phase"},

     {"RMR_KREP_MAX_1", 0.8, 0, 0,
     "Max lesion L1 repair rate",
     "Maximum coefficient of rate of repair of L1 potentially lethal lesions: Krepair_max (Curtis's epsilon_PL)\n\
     The true repair rate varies linearly from Krepair_base to Krepair_max over the S phase"},

     {"RMR_KMIS_1_1", 0.0278, 0, 0,
     "Lesion misrepair rate to L2b",
     "Coefficient of rate of misrepair of L1 potentially lethal lesions to L2b: (Curtis's epsilon_2PL)"},

     {"RMR_KMIS_2_1", 0.0278, 0, 0,
     "Lesion misrepair rate to L2c",
     "Coefficient of rate of misrepair of L1 potentially lethal lesions to L2c: (Curtis's epsilon_2PL)"},

     {"RMR_KCP_1", 0.13, 0, 0,
     "Checkpoint time limit factor",
     "Factor for computing maximum time spent in the checkpoint, as a function of # of L1 lesions"},

// Metabolism parameters

      {"USE_METABOLISM", 1,0,1,
       "Use glucose metabolism",
       "If metabolism is turned on, cell cycle is used, and lactate is simulated.  If metabolism is turned off, lactate is not simulated"},

      {"N_GA_1", 2, 0, 0,
      "ATP moles produced per glucose mole",
      "Number of ATP moles produced by the glycolysis of one glucose mole"},

      {"N_PA_1", 14, 0, 0,
      "ATP moles produced per pyruvate mole",
      "Number of ATP moles produced by the oxidation of one pyruvate mole"},

      {"N_GI_1", 0.4, 0, 0,
      "Intermediate moles produced per glucose mole",
      "Number of moles of anabolic intermediates produced the glycolysis of one glucose mole"},

      {"N_PI_1", 0.4, 0, 0,
      "Intermediate moles produced per pyruvate mole",
      "Number of moles of anabolic intermediates produced the oxidation of one pyruvate mole"},

      {"N_PO_1", 3, 0, 0,
      "Oxygen moles consumed per pyruvate mole",
      "Number of moles of oxygen consumed the oxidation of one pyruvate mole"},

      {"K_H1_1", 3, 0, 0,
      "Exponent N_H in Hss equation",
      "With x = C_O/C_O_max, steady-state HIF-1 level is Hss = (1-x)^N_H, where C_O = oxygen concentration"},

      {"K_H2_1", 0.00005, 0, 0,
      "HIF-1 rate constant K_H",
      "The rate of change of HIF-1 level H is given by: dH/dt = K_H*(Hss - H), where C_O = oxygen concentration, Hss is steady-state level for C_O"},

      {"K_HB_1", 0.2, 0, 0,
      "K_HB",
      "Glycolysis rate = K_HA*(1 + K_HB*H)*C_G^N/(C_G^N + Km^N)\n\
        where: H = HIF-1 level, C_G = glucose concentration, K_HA is maximum glucose consumption rate when H=0, Km and N are the glucose consumption Hill function parameters"},

      {"K_PDK_1", 4.63e-5, 0, 0,
      "K_PDK",
      "Representing PDK1 factor level by y in the range (PDKmin,1), the rate of change of y is: dy/dt = -K_PDK*(y - 1 + cH) where H = HIF-1 level, and c = 1 - PDKmin"},

       {"PDKMIN_1", 0.3, 0, 0,
       "PDKmin",
       "Minimum value of the PDK1 factor.\n\
       Representing PDK1 factor level by y in the range (PDKmin,1), the rate of change of y is: dy/dt = -K_PDK*(y - 1 + cH) where H = HIF-1 level, and c = 1 - PDKmin"},

       {"C_O2_NORM_1", 0.05, 0, 0,
       "Nominal normal IC O2 concentration",
       "Nominal normal IC O2 concentration, used to set normal metabolic rates for unconstrained growth"},

       {"C_G_NORM_1", 2.5, 0, 0,
       "Nominal normal IC glucose concentration",
       "Nominal normal IC glucose concentration, used to set normal metabolic rates for unconstrained growth"},

       {"C_L_NORM_1", 1.0, 0, 0,
       "Nominal normal IC lactate concentration",
       "Nominal normal IC lactate concentration, used to set normal metabolic rates for unconstrained growth"},

      {"ATP_S_1", 0.3, 0, 0,
      "ATP production threshold for survival (fraction of peak)",
      "Cell death occurs when the ATP production rate falls below the fraction ATP_S of the maximum (no nutrient constraints) production rate"},

      {"ATP_G_1", 0.75, 0, 0,
       "ATP production threshold for growth (fraction of peak)",
       "Cell growth stops when the ATP production rate falls below the fraction ATP_G of the maximum (no nutrient constraints) production rate.\n\
        Intermediates production from glycolysis and pyruvate ramps down when ATP rate is below ATPramp*ATP_G, to 0 when ATP rate = ATP_G."},

      {"ATP_RAMP_1", 1.3, 0, 0,
        "Ramp factor for reducing r_G, r_P based on ATP",
       "Fractional rates of production of intermediates from glycolysis (r_G) and pyruvate (r_P) are reduced as ATP production rate falls.\n\
        If the normal ATP rate is r_An, and the current rate is r_A, then the fractional production rate f_A = r_A/r_An.\n\
        A multiplying factor on both r_G and r_P, f, is calculated as a ramp from 0 - 1 as f_A varies from ATP_G - ATPramp*ATP_G\n\
        i.e. the reduction starts when f_A = ATPramp*ATP_G (f = 1), and when f_A = ATP_G, f = 0.  f_G = f*f_Gn, f_P = f*f_Pn"},

      {"K_PL_1", 0.001, 0, 0,
      "Pyruvate -> lactate rate constant",
      "The forward rate constant K_PL of the pyruvate-lactate reaction, i.e. the rate constant for conversion of pyruvate to lactate"},

      {"K_LP_1", 0.001, 0, 0,
      "Lactate -> pyruvate rate constant",
      "The reverse rate constant K_LP of the pyruvate-lactate reaction, i.e. the rate constant for conversion of lactate to pyruvate"},

       {"PYRUVATE_MM_KM_1", 20, 0, 0,
       "Pyruvate Michaelis-Menten Km (uM)",
       "Pyruvate Michaelis-Menten Km (uM). fMM = C_P/(Km + C_P) is a multiple of pyruvate oxidation rate, ensuring that the rate -> 0 as C_P -> 0"},

//       {"N_GA_2", 2, 0, 0,
//       "ATP moles produced per glucose mole",
//       "Number of ATP moles produced by the glycolysis of one glucose mole"},

//       {"N_PA_2", 14, 0, 0,
//       "ATP moles produced per pyruvate mole",
//       "Number of ATP moles produced by the oxidation of one pyruvate mole"},

//       {"N_GI_2", 0.6, 0, 0,
//       "Intermediate moles produced per glucose mole",
//       "Number of moles of anabolic intermediates produced the glycolysis of one glucose mole"},

//       {"N_PI_2", 0.6, 0, 0,
//       "Intermediate moles produced per pyruvate mole",
//       "Number of moles of anabolic intermediates produced the oxidation of one pyruvate mole"},

//       {"N_PO_2", 3, 0, 0,
//       "Oxygen moles consumed per pyruvate mole",
//       "Number of moles of oxygen consumed the oxidation of one pyruvate mole"},

//       {"K_H1_2", 140, 0, 0,
//       "K_H1",
//       "The rate of change of HIF-1 level H is given by: dH/dt = K_H2*(1 - H*exp(K_H1*C_O)), where C_O = oxygen concentration"},

//       {"K_H2_2", 0.001, 0, 0,
//       "K_H2",
//        "The rate of change of HIF-1 level H is given by: dH/dt = K_H2*(1 - H*exp(K_H1*C_O)), where C_O = oxygen concentration"},

//       {"K_HB_2", 0.2, 0, 0,
//       "K_HB",
//        "Glycolysis rate = K_HA*(1 + K_HB*H)*C_G^N/(C_G^N + Km^N)\n\
//         where: H = HIF-1 level, C_G = glucose concentration, K_HA is maximum glucose consumption rate when H=0, Km and N are the glucose consumption Hill function parameters"},

//       {"K_PDK_2", 4.63e-5, 0, 0,
//       "K_PDK",
//       "Representing PDK1 factor level by y in the range (0,1), the rate of change of y is: dy/dt = -K_PDK*(y - 1 + H) where H = HIF-1 level"},

//        {"PDKMIN_2", 0.3, 0, 0,
//        "PDKmin",
//        "Minimum value of the PDK1 factor"},

//        {"C_O2_NORM_2", 0.05, 0, 0,
//        "Nominal normal IC O2 concentration",
//        "Nominal normal IC O2 concentration, used to set normal metabolic rates for unconstrained growth"},

//        {"C_G_NORM_2", 2.5, 0, 0,
//        "Nominal normal IC glucose concentration",
//        "Nominal normal IC glucose concentration, used to set normal metabolic rates for unconstrained growth"},

//        {"C_L_NORM_2", 0.0, 0, 0,
//        "Nominal normal IC lactate concentration",
//        "Nominal normal IC lactate concentration, used to set normal metabolic rates for unconstrained growth"},

////        {"CO_H_2", 0.005, 0, 0,
////         "Oxygen H threshold for Ofactor",
////         "When C_O2 < CO_H both f_G and f_P are reduced by Ofactor, which ranges from 0 to 1 when CO_L < C_O2 < CO_H.  CO_L = 0.8*CO_H"},

////        {"CG_H_2", 0.05, 0, 0,
////         "Glucose H threshold for Gfactor",
////         "When C_G < CG_H f_G is reduced by Gfactor, which ranges from 0 to 1 when CG_L < C_G < CG_H.  CG_L = 0.8*CG_H"},

//       {"ATP_S_2", 0.3, 0, 0,
//       "ATP production threshold for survival (fraction of peak)",
//        "Cell death occurs when the ATP production rate falls below the fraction ATP_S of the maximum (no nutrient constraints) production rate"},

//        {"ATP_G_2", 0.55, 0, 0,
//        "ATP production threshold for growth (fraction of peak)",
//        "Cell growth stops when the ATP production rate falls below the fraction ATP_G of the maximum (no nutrient constraints) production rate"},

//        {"ATP_RAMP_2", 1.5, 0, 0,
//          "Ramp factor for reducing r_G, r_P based on ATP",
//         "Fractional rates of production of intermediates from glycolysis (r_G) and pyruvate (r_P) are reduced as ATP production rate falls.\n\
//          If the normal ATP rate is r_An, and the current rate is r_A, then the fractional production rate f_A = r_A/r_An.\n\
//          A multiplying factor on both r_G and r_P, f, is calculated as a ramp from 0 - 1 as f_A varies from ATP_G - ATPramp*ATP_G\n\
//          i.e. the reduction starts when f_A = ATPramp*ATP_G (f = 1), and when f_A = ATP_G, f = 0.  f_G = f*f_Gn, f_P = f*f_Pn"},

//       {"K_PL_2", 0.001, 0, 0,
//       "Pyruvate -> lactate rate constant",
//       "The forward rate constant of the pyruvate-lactate reaction, i.e. the rate constant for conversion of pyruvate to lactate"},

//       {"K_LP_2", 0.001, 0, 0,
//       "Lactate -> pyruvate rate constant",
//        "The reverse rate constant of the pyruvate-lactate reaction, i.e. the rate constant for conversion of lactate to pyruvate"},

//        {"PYRUVATE_MM_KM_2", 20, 0, 0,
//        "Pyruvate Michaelis-Menten Km",
//        "Pyruvate Michaelis-Menten Km (uM). fMM = C_P/(Km + C_P) is a multiple of pyruvate oxidation rate, ensuring that the rate -> 0 as C_P -> 0"},



{"HYPOXIA_1", 0.1, 0, 0,
"Hypoxia threshold 1",
"Hypoxia threshold 1"},

{"HYPOXIA_2", 1.0, 0, 0,
"Hypoxia threshold 2",
"Hypoxia threshold 2"},

{"HYPOXIA_3", 4.0, 0, 0,
"Hypoxia threshold 3",
"Hypoxia threshold 3"},

{"HYPOXIA_THRESHOLD", 4.0, 0, 0,
"Hypoxia threshold",
"Hypoxia threshold"},

{"GROWTH_FRACTION_1", 0.25, 0, 0,
"Growth fraction threshold 1",
"Growth fraction threshold 1"},

{"GROWTH_FRACTION_2", 0.1, 0, 0,
"Growth fraction threshold 2",
"Growth fraction threshold 2"},

{"GROWTH_FRACTION_3", 0.01, 0, 0,
"Growth fraction threshold 3",
"Growth fraction threshold 3"},

{"DRUG_THRESHOLD", 1.0e-6, 0, 0,
 "Drug Threshold",
 "Threshold drug concentration - when all intracellular and extracellular concentrations fall below this level, the drug concentrations everywhere are set to zero"},

{"DRUG_LABEL_THRESHOLD", 0, 0, 0,
"Label Threshold",
"Threshold label-drug concentration - when a labelling drug (e.g. EDU) is used, this is the threshold for a cell to be considered as labelled"},

{"SPCRAD", 200.0, 0, 0,
"Spectral radius",
"Spectral radius value used by RKC solver"},

 {"SAVE_FACS_DATA",0,0,1,
  "Save FACS data",
  "Save data for FACS at a specified interval"},

 {"SAVE_FACS_DATA_FILE_NAME",0,0,0,
  "facs_data",
  "Base file name for saving FACS data"},

 {"SAVE_FACS_DATA_TSTART",0,0,0,
  "Start time",
  "Start time for saving FACS data"},

 {"SAVE_FACS_DATA_INTERVAL",0,0,0,
  "Interval",
  "Time interval for saving FACS data"},

 {"SAVE_FACS_DATA_NUMBER",1,0,0,
  "Number",
  "Number of times to save FACS data"},


//{"USE_EXTRA", 0, 0, 1,
//"Use extra conc",
//"Use extracellular O2 and glucose concentrations to determine cell death"},

//{"USE_RELAX", 1, 0, 1,
//"Use O2 relaxation solver",
//"Use over- and under-relaxation to solve reaction-diffusion for oxygen"},

//{"USE_PAR_RELAX", 1, 0, 1,
//"Use parallel O2 relaxation solver",
//"Use over- and under-relaxation to solve reaction-diffusion for oxygen, with parallelized over-relaxation"},

//{"FD_SOLVER_1", 1, 0, 1,
//"Use FD solver?",
//"Use the FD solver in the far field"},

//{"USE_DROP", 0, 0, 1,
//"Account for drop deformation",
//"Account for drop deformation when it is released to sit at the bottom of the well"},

//{"NDROP", 1000, 0, 0,
//"Dropping cell count",
//"Number of cells in the spheroid when it is dropped."},

//{"DROP_ALPHA", 0.4, 0, 0,
//"Contact_diameter/diameter",
//"Drop parameter alpha = initial (surface contact diameter)/(blob diameter).  Must be < 1."},

//{"DROP_BETA", 0.6, 0, 0,
//"Height/diameter",
//"Drop parameter beta = initial (blob height)/(blob diameter).  Must be < 1."},

//    {"SAVE_PROFILE_DATA",0,0,1,
//     "Save profile data",
//     "Save data for profile plots at a specified interval"},

//    {"SAVE_PROFILE_DATA_FILE_NAME",0,0,0,
//     "profile_data",
//     "Base file name for saving profile data"},

//    {"SAVE_PROFILE_DATA_INTERVAL",0,0,0,
//     "Interval",
//     "Time interval for saving profile data"},

//    {"SAVE_PROFILE_DATA_NUMBER",1,0,0,
//     "Number",
//     "Number of times to save profile data"},

//    {"SAVE_SLICE_DATA",0,0,1,
//     "Save slice data",
//     "Save data for z-slices at a specified interval"},

//    {"SAVE_SLICE_DATA_FILE_NAME",0,0,0,
//     "slice_data",
//     "Base file name for saving slice data"},

//    {"SAVE_SLICE_DATA_INTERVAL",0,0,0,
//     "Interval",
//     "Time interval for saving slice data"},

//    {"SAVE_SLICE_DATA_NUMBER",1,0,0,
//     "Number",
//     "Number of times to save slice data"},

// This is the end of the parameters that are actually read by the DLL
// Entries after this point are QMyLabel dummies, to enable display of explanatory info  - no input data is transmitted,
// followed by the list of time-series and profile plots selected for this run.

{"DUMMY_HYPOXIA_THRESHOLD", 0, 0, 0,
"Hypoxia threshold",
"Select the intracellular O2 level below which the cell is counted as hypoxic"},

{"DUMMY_GROWTH_FRACTION", 0, 0, 0,
"Growth fraction",
"Select the threshold fraction of average growth rate (i.e. with no nutrient limits) used to count slow-growing cells"},

// Time-series plots
    {"nlive",                     1, 0,1,"","Number of live cells"},
    {"nviable",                   0, 0,1,"","Number of viable cells"},
    {"nanoxiadead",               0, 0,1,"","Total number of cells that have been killed by anoxia"},
    {"naglucosiadead",            0, 0,1,"","Total number of cells that have been killed by aglucosia"},
    {"ndrugAdead",                0, 0,1,"","Total number of cells that have been killed by drugA"},
    {"ndrugBdead",                0, 0,1,"","Total number of cells that have been killed by drugB"},
    {"nradiationdead",            0, 0,1,"","Total number of cells that have been killed by radiation"},
    {"nanoxiatagged",             0, 0,1,"","Current number of cells tagged to die by anoxia"},
    {"naglucosiatagged",          0, 0,1,"","Current number of cells tagged to die by aglucosia"},
    {"ndrugAtagged",              0, 0,1,"","Current number of cells tagged to die by drugA"},
    {"ndrugBtagged",              0, 0,1,"","Current number of cells tagged to die by drugB"},
    {"nradiationtagged",          0, 0,1,"","Current number of cells tagged to die by radiation"},
    {"hypoxicfraction",           0, 0,1,"","Fraction of cells with oxygen level below the specified threshold for hypoxia"},
    {"clonohypoxicfraction",      0, 0,1,"","Fraction of clonogenic cells with oxygen level below the specified threshold for hypoxia"},
    {"growthfraction",            0, 0,1,"","Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits"},
    {"platingefficiency",         1, 0,1,"","Percentage of live cells that are viable"},
    {"ECoxygen",                  1, 0,1,"","EC concentration of oxygen in the medium (bottom)"},
    {"ECglucose",                 1, 0,1,"","EC concentration of glucose in the medium (bottom)"},
    {"EClactate",                 1, 0,1,"","EC concentration of lactate in the medium (bottom)"},
    {"ECdrugA",                   0, 0,1,"","EC concentration of drug A in the medium (bottom)"},
    {"ECdrugAmet1",               0, 0,1,"","EC concentration of drug A metabolite 1 in the medium (bottom)"},
    {"ECdrugAmet2",               0, 0,1,"","EC concentration of drug A metabolite 2 in the medium (bottom)"},
    {"ECdrugB",                   0, 0,1,"","EC concentration of drug B in the medium (bottom)"},
    {"ECdrugBmet1",               0, 0,1,"","EC concentration of drug B metabolite 1 in the medium (bottom)"},
    {"ECdrugBmet2",               0, 0,1,"","EC concentration of drug B metabolite 2 in the medium (bottom)"},
    {"ICoxygen",                  1, 0,1,"","IC concentration of oxygen"},
    {"ICglucose",                 1, 0,1,"","IC concentration of glucose"},
    {"IClactate",                 1, 0,1,"","IC concentration of lactate"},
    {"ICpyruvate",                1, 0,1,"","IC concentration of pyruvate"},
    {"ICdrugA",                   0, 0,1,"","IC concentration of drug A"},
    {"ICdrugAmet1",               0, 0,1,"","IC concentration of drug A metabolite 1"},
    {"ICdrugAmet2",               0, 0,1,"","IC concentration of drug A metabolite 2"},
    {"ICdrugB",                   0, 0,1,"","IC concentration of drug B"},
    {"ICdrugBmet1",               0, 0,1,"","IC concentration of drug B metabolite 1"},
    {"ICdrugBmet2",               0, 0,1,"","IC concentration of drug B metabolite 2"},
    {"Medoxygen",                 0, 0,1,"","Average medium concentration of oxygen"},
    {"Medglucose",                0, 0,1,"","Average medium concentration of glucose"},
    {"Medlactate",                0, 0,1,"","Average medium concentration of lactate"},
    {"MeddrugA",                  0, 0,1,"","Average medium concentration of drug A"},
    {"MeddrugAmet1",              0, 0,1,"","Average medium concentration of drug A metabolite 1"},
    {"MeddrugAmet2",              0, 0,1,"","Average medium concentration of drug A metabolite 2"},
    {"MeddrugB",                  0, 0,1,"","Average medium concentration of drug B"},
    {"MeddrugBmet1",              0, 0,1,"","Average medium concentration of drug B metabolite 1"},
    {"MeddrugBmet2",              0, 0,1,"","Average medium concentration of drug B metabolite 2"},
    {"doublingtime",              0, 0,1,"","Average doubling time"},
    {"Grate",                     1, 0,1,"","Normalised glycolysis rate"},
    {"Prate",                     1, 0,1,"","Normalised pyruvate utilisation rate"},
    {"Arate",                     1, 0,1,"","Normalised ATP production rate"},
    {"Irate",                     1, 0,1,"","Normalised rate of production of anabolic intermediates"},
    {"f_G",                       0, 0,1,"","f_G"},
    {"f_P",                       0, 0,1,"","f_P"},
    {"HIF-1",                     1, 0,1,"","HIF-1"},
    {"PDK1",                      1, 0,1,"","PDK1"},
    {"dividerate",                0, 0,1,"","# divided/hour"},


// Profile plots
    {"MULTI",                     1, 0,1,"","Selected constituent on a line through the blob centre"},
//    {"CFSE",                      0, 0,1,"","Extracellular CFSE concentration on a line through the blob centre"},
    {"Oxygen",                    0, 0,1,"","Extracellular oxygen concentration on a line through the blob centre"},
    {"Glucose",                   0, 0,1,"","Extracellular glucose concentration on a line through the blob centre"},
//    {"Tracer",                    0, 0,1,"","Extracellular tracer concentration on a line through the blob centre"},
    {"Drug_A",                    0, 0,1,"","Extracellular drug A concentration on a line through the blob centre"},
    {"Drug_A_metab1",             0, 0,1,"","Extracellular drug A metabolite 1 concentration on a line through the blob centre"},
    {"Drug_A_metab2",             0, 0,1,"","Extracellular drug A metabolite 2 concentration on a line through the blob centre"},
    {"Drug_B",                    0, 0,1,"","Extracellular drug Bconcentration on a line through the blob centre"},
    {"Drug_B_metab1",             0, 0,1,"","Extracellular drug B metabolite 1 concentration on a line through the blob centre"},
    {"Drug_B_metab2",             0, 0,1,"","Extracellular drug B metabolite 2 concentration on a line through the blob centre"},
/*
    {"IC_MULTI",                  1, 0,1,"","Selected constituent on a line through the blob centre"},
    {"IC_Oxygen",                 0, 0,1,"","Intracellular oxygen concentration on a line through the blob centre"},
    {"IC_Glucose",                0, 0,1,"","Intracellular glucose concentration on a line through the blob centre"},
    {"IC_Drug_A",                 0, 0,1,"","Intracellular drug A concentration on a line through the blob centre"},
    {"IC_Drug_A_metab1",          0, 0,1,"","Intracellular drug A metabolite 1 concentration on a line through the blob centre"},
    {"IC_Drug_A_metab2",          0, 0,1,"","Intracellular drug A metabolite 2 concentration on a line through the blob centre"},
    {"IC_Drug_B",                 0, 0,1,"","Intracellular drug Bconcentration on a line through the blob centre"},
    {"IC_Drug_B_metab1",          0, 0,1,"","Intracellular drug B metabolite 1 concentration on a line through the blob centre"},
    {"IC_Drug_B_metab2",          0, 0,1,"","Intracellular drug B metabolite 2 concentration on a line through the blob centre"},
    {"IC_CFSE",                   0, 0,1,"","CFSE concentration on a line through the blob centre"},
    {"IC_growthrate",             0, 0,1,"","Cell growth rate on a line through the blob centre"},
    {"IC_cellvolume",             0, 0,1,"","Cell volume fraction on a line through the blob centre"},
    {"IC_O2byvolume",             0, 0,1,"","Cell volume fraction on a line through the blob centre"},
// Distribution plots
//    {"Oxygen",                    0, 0,1,"","Probability distribution of extracellular oxygen concentration"},
//    {"cellvolume",                0, 0,1,"","Probability distribution of cell volume fraction"}
*/

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}

    nInfolabel = sizeof(label_info)/sizeof(INFOSTRUCT);
    workingInfolabelList = new INFOSTRUCT[nInfolabel];
    for (int i=0; i<nInfolabel; i++) {
        workingInfolabelList[i] = label_info[i];
    }
    /*
    nInfocheckbox = sizeof(checkbox_info)/sizeof(INFOSTRUCT);
    workingInfocheckboxList = new INFOSTRUCT[nInfocheckbox];
    for (int i=0; i<nInfocheckbox; i++) {
        workingInfocheckboxList[i] = checkbox_info[i];
    }
    */
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}

void Params::get_labeltag(int i, QString *tag)
{
    *tag = workingInfolabelList[i].tag;
}

void Params::infoLabelInfo(QString tag, QString *info)
{
    for (int i=0; i<nInfolabel; i++) {
        if (tag == workingInfolabelList[i].tag) {
            *info = workingInfolabelList[i].info;
            return;
        } else {
            *info = "";
        }
    }
}

