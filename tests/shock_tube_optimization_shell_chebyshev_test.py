import sys, os
sys.path.append('../../') #get rid of this at some point with central test script or when package is built
os.chdir('../../')

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.optimization.opt_runner as opt
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
#import MSI.optimization.shock_tube_optimization_shell_six_param_fit as stMSIspf
import MSI.optimization.shock_tube_optimization_shell_chebyshev as stMSIcheb
import cantera as ct
import pandas as pd
import numpy as np
import MSI.utilities.plotting_script as plotter
import MSI.utilities.post_processor as post_processor




files_to_include = [['Hong_0.yaml'],
                    ['Hong_2.yaml'],
                    ['Hong_3.yaml'],
                    ['Hong_4.yaml','Hong_4_abs.yaml'],
                    ['Hong_1.yaml'],
                    ['Troe_4.yaml','Troe_4_abs.yaml'],
                    ['Troe_5.yaml','Troe_5_abs.yaml'],
                    ['Troe_6.yaml','Troe_6_abs.yaml'],
                    ['Troe_7.yaml','Troe_7_abs.yaml'],
                    ['Troe_8.yaml','Troe_8_abs.yaml'],            
                    ['Hong_5.yaml','Hong_5_abs.yaml']]
                    



#['Hong_6_high_temp_large_uncertainty.yaml', 'Hong_6_high_temp_abs.yaml']] 
                                                      

#files_to_include = [['Hong_0.yaml']]
                    



                                                      
numer_of_iterations = 2
cti_file = 'FFCM1_custom_cheb_extra_zeros_new.cti'
#cti_file = 'FFCM1_custom_cheb_extra_zeros_new_1_colliders.cti'


working_directory = 'MSI/data/chebyshev_data'
reaction_uncertainty_csv = 'FFCM1_reaction_uncertainty.csv'

master_reaction_equation_cti_name= 'master_reactions_FFCM1_cheb_extra_zeros_new.cti'

#rate_constant_target_value_data = 'burke_target_value_single_reactions.csv'

#this would be an empty string '' if you do not want to include it 
run_with_k_target_values = 'On'
master_equation_reactions = ['H2O2 + OH <=> H2O + HO2',
                             '2 HO2 <=> H2O2 + O2',
                             'HO2 + OH <=> H2O + O2',
                             '2 OH <=> H2O + O',
                             'CH3 + HO2 <=> CH4 + O2',
                             'CH3 + HO2 <=> CH3O + OH']

master_index = [2,3,4,5,6,7]

master_equation_uncertainty_df = pd.read_csv('MSI/data/test_data/six_parameter_fit_uncertainty_df.csv')

#rate_constant_target_value_data_for_plotting = 'FFCM1_target_reactions_1.csv'
#rate_constant_target_value_data = 'FFCM1_target_reactions_1.csv'
#rate_constant_target_value_data_extra = 'FFCM1_target_reactions_1.csv'



rate_constant_target_value_data_for_plotting = 'FFCM1_target_reactions_1_plotting.csv'
rate_constant_target_value_data = 'FFCM1_target_reactions_1.csv'
rate_constant_target_value_data_extra = 'FFCM1_target_reactions_extra_data.csv'


cheb_sensitivity_dict = {'2 HO2 <=> H2O2 + O2': [np.array([[-0.842718,0],
         [ 0.443424,0],
         [ 0.168757,0],
         [ 0.0322031,0],
         [-0.0057967,0],
         [-0.00719272,0],
         [-0.00394207,0],
         [-0.00143191,0],
         [ 6.63935e-06,0],
         [0.000493172,0],
         [0.000712372,0],
         [0.000523586,0],
         [0.000554322,0],
         [0.000186269,0],
         [0.000466622,0]]), np.array([[-2.74126,0],
         [-0.7876,0],
         [-0.350124,0],
         [-0.100554,0],
         [0.0124844,0],
         [ 0.0220203,0],
         [ 0.0161857,0],
         [0.00972219,0],
         [0.00522206,0],
         [0.00261056,0],
         [ 0.00115085,0],
         [0.000485739,0],
         [9.53231e-05,0],
         [5.45247e-05,0],
         [-0.000140707,0]]), np.array([[-2.46305,0],
         [-0.725156,0],
         [-0.301131,0],
         [-0.0742016,0],
         [ 0.00162991,0],
         [0.00893466,0],
         [0.00580178,0],
         [ 0.00251803,0],
         [0.000429041,0],
         [-0.000388348,0],
         [-0.000795072,0],
         [-0.000616906,0],
         [-0.00069521,0],
         [-0.000230894,0],
         [-0.000608478,0]]), np.array([[-0.141337,0],
         [-0.218581,0],
         [-0.138015,0],
         [-0.066208,0],
         [-0.0205548,0],
         [-0.00630092,0],
         [-0.00112983,0],
         [ 0.0005282,0],
         [ 0.000938152,0],
         [0.000826389,0],
         [0.000683619,0],
         [0.000431331,0],
         [ 0.000361901,0],
         [0.000127995,0],
         [0.000251668,0]]), np.array([[1.83004,0],
         [-1.4719,0],
         [-0.359648,0],
         [ 0.0334891,0],
         [0.0598449,0],
         [ 0.0354098,0],
         [0.0140112,0],
         [0.00245409,0],
         [-0.00294833,0],
         [-0.00405897,0],
         [-0.00431978,0],
         [-0.00296089,0],
         [-0.0028622,0],
         [-0.000978951,0],
         [-0.00225978,0]]), np.array([[ 0.0284025,0],
         [-0.0268788,0],
         [-0.00416215,0],
         [0.00197177,0],
         [0.00118763,0],
         [0.000558997,0],
         [0.000162112,0],
         [-2.1396e-05,0],
         [-9.46499e-05,0],
         [-9.73581e-05,0],
         [-9.16351e-05,0],
         [-6.05716e-05,0],
         [-5.57102e-05,0],
         [-1.92296e-05,0],
         [-4.23241e-05,0]]), np.array([[-0.539735,0],
         [-0.0221098,0],
         [0.0115172,0],
         [ 0.0230455,0],
         [0.0208784,0],
         [0.0124353,0],
         [0.0065188,0],
         [0.00312159,0],
         [0.00130829,0],
         [ 0.000468186,0],
         [5.36476e-05,0],
         [-4.53044e-05,0],
         [-0.000136469,0],
         [-4.21932e-05,0],
         [-0.000159041,0]])], '2 OH <=> H2O + O': [np.array([[-1.01166,0],
         [0.73383,0],
         [ 0.145363,0],
         [-0.00319219,0],
         [-0.00528572,0],
         [-0.00453066,0],
         [-0.0024436,0],
         [-0.000979143,0],
         [-0.000191809,0],
         [0.00013987,0],
         [0.000251601,0],
         [0.000234661,0],
         [0.000206912,0],
         [0.000135437,0],
         [0.000120021,0],
         [ 4.22412e-05,0],
         [8.93489e-05,0]]), np.array([[-1.19558,0],
         [-1.17259,0],
         [-0.369084,0],
         [-0.0687507,0],
         [0.00436644,0],
         [0.0123823,0],
         [0.0088778,0],
         [0.00483607,0],
         [0.00215561,0],
         [0.000725566,0],
         [3.46434e-05,0],
         [-0.000194751,0],
         [-0.000282099,0],
         [-0.00021042,0],
         [-0.000216972,0],
         [-7.46855e-05,0],
         [-0.000180641,0]]), np.array([[-0.987634,0],
         [-0.802963,0],
         [-0.112675,0],
         [ 0.0632623,0],
         [ 0.0395816,0],
         [ 0.0212316,0],
         [0.00918935,0],
         [0.00284805,0],
         [-4.22451e-05,0],
         [-0.00105373,0],
         [-0.0012725,0],
         [-0.0010703,0],
         [-0.000882971,0],
         [-0.000563633,0],
         [-0.000481131,0],
         [-0.000170513,0],
         [-0.000346272,0]]), np.array([[-0.102416,0],
         [-0.167246,0],
         [-0.101446,0],
         [-0.0414027,0],
         [-0.00418504,0],
         [ 0.00173931,0],
         [ 0.00211631,0],
         [0.0014449,0],
         [0.000812592,0],
         [0.000403986,0],
         [0.000171349,0],
         [5.96113e-05,0],
         [3.09839e-06,0],
         [-9.54566e-06,0],
         [-2.24133e-05,0],
         [-7.08954e-06,0],
         [-2.5505e-05,0]]), np.array([[ 1.97118,0],
         [-2.02716,0],
         [-0.0105087,0],
         [0.115295,0],
         [0.0172879,0],
         [ 0.0030159,0],
         [-0.000196345,0],
         [-0.000552792,0],
         [-0.000327369,0],
         [-0.000108848,0],
         [2.46866e-05,0],
         [6.8754e-05,0],
         [8.78688e-05,0],
         [6.4049e-05,0],
         [6.61824e-05,0],
         [ 2.25656e-05,0],
         [5.56654e-05,0]]), np.array([[ 0.215294,0],
         [-0.272304,0],
         [ 0.0509158,0],
         [0.00988242,0],
         [-0.000309211,0],
         [-0.00128376,0],
         [-0.000786195,0],
         [-0.0002875,0],
         [ 1.75342e-06,0],
         [0.000118814,0],
         [0.000151274,0],
         [0.000130815,0],
         [0.000111053,0],
         [7.17137e-05,0],
         [6.26303e-05,0],
         [2.2064e-05,0],
         [4.61163e-05,0]]), np.array([[-0.0233324,0],
         [-0.0324617,0],
         [-0.0138609,0],
         [-0.000993681,0],
         [0.00399712,0],
         [0.00289286,0],
         [0.00159349,0],
         [0.000741607,0],
         [ 0.000275427,0],
         [5.65035e-05,0],
         [-3.66328e-05,0],
         [-5.79734e-05,0],
         [-6.19733e-05,0],
         [-4.30751e-05,0],
         [-4.09011e-05,0],
         [-1.42749e-05,0],
         [-3.20888e-05,0]]), np.array([[-0.170542,0],
         [-0.304177,0],
         [-0.230334,0],
         [-0.145391,0],
         [-0.0723741,0],
         [-0.0343758,0],
         [-0.0150007,0],
         [-0.00568862,0],
         [-0.00148832,0],
         [ 0.000169741,0],
         [0.000727856,0],
         [0.000735619,0],
         [0.000661361,0],
         [0.000435364,0],
         [0.000383713,0],
         [0.000135782,0],
         [0.000282886,0]]), np.array([[-0.139339,0],
         [-0.243557,0],
         [-0.180303,0],
         [-0.108115,0],
         [-0.0460284,0],
         [-0.0196129,0],
         [-0.00740786,0],
         [-0.00205729,0],
         [9.34898e-05,0],
         [ 0.000766094,0],
         [ 0.000871709,0],
         [ 0.00071246,0],
         [0.000571344,0],
         [0.00036043,0],
         [0.000300675,0],
         [0.000107184,0],
         [0.000211316,0]]), np.array([[-0.0258964,0],
         [-0.0483299,0],
         [-0.039072,0],
         [-0.0274521,0],
         [-0.0165071,0],
         [-0.00849802,0],
         [-0.00396711,0],
         [-0.00163896,0],
         [-0.000526986,0],
         [-5.52984e-05,0],
         [ 0.000123945,0],
         [0.000151126,0],
         [0.000147034,0],
         [9.93398e-05,0],
         [9.05221e-05,0],
         [3.18676e-05,0],
         [6.86076e-05,0]])], 'CH3 + HO2 <=> CH3O + OH': [np.array([[ 0.988823,0],
         [0.00838564,0],
         [-0.0130431,0],
         [-0.00479895,0],
         [-0.00198649,0],
         [-0.000876941,0],
         [-0.000403419,0],
         [-0.000190643,0],
         [-9.2226e-05,0],
         [-4.49547e-05,0],
         [-2.26046e-05,0],
         [-1.08911e-05,0],
         [-6.00769e-06,0],
         [-2.19272e-06,0],
         [-2.42801e-06,0]]), np.array([[0.0903702,0],
         [0.098767,0],
         [0.0459482,0],
         [0.0169057,0],
         [0.00699799,0],
         [0.00308928,0],
         [0.00142116,0],
         [0.000671594,0],
         [0.000324893,0],
         [0.000158366,0],
         [7.96313e-05,0],
         [3.83672e-05,0],
         [2.11639e-05,0],
         [7.72447e-06,0],
         [8.55337e-06,0]])], 'CH3 + HO2 <=> CH4 + O2': [np.array([[-1.36438,0],
         [1.15445,0],
         [-0.000316787,0],
         [-0.000116556,0],
         [-4.82473e-05,0],
         [-2.12989e-05,0],
         [-9.79813e-06,0],
         [-4.63027e-06,0],
         [-2.23996e-06,0],
         [-1.09184e-06,0],
         [-5.49017e-07,0],
         [-2.64519e-07,0],
         [-1.45915e-07,0],
         [-5.32544e-08,0],
         [-5.89714e-08,0]]), np.array([[-1.96949,0],
         [-1.42416,0],
         [-0.112362,0],
         [-0.0413414,0],
         [-0.0171129,0],
         [-0.00755455,0],
         [-0.00347532,0],
         [-0.00164232,0],
         [-0.000794496,0],
         [-0.000387269,0],
         [-0.000194731,0],
         [-9.38236e-05,0],
         [-5.17542e-05,0],
         [-1.88895e-05,0],
         [-2.09165e-05,0]]), np.array([[0.0912312,0],
         [0.0976024,0],
         [0.0470496,0],
         [0.017311,0],
         [0.00716574,0],
         [0.00316334,0],
         [0.00145523,0],
         [0.000687693,0],
         [0.000332681,0],
         [0.000162162,0],
         [8.15402e-05,0],
         [3.9287e-05,0],
         [2.16712e-05,0],
         [7.90965e-06,0],
         [8.7584e-06,0]])], 'H2O2 + OH <=> H2O + HO2': [np.array([[-0.846524,0],
         [ 0.447649,0],
         [0.174235,0],
         [ 0.0271819,0],
         [-0.00729249,0],
         [-0.00731393,0],
         [-0.00345375,0],
         [-0.000804495,0],
         [0.000601499,0],
         [ 0.000950295,0],
         [ 0.001073,0],
         [0.000747398,0],
         [ 0.00074263,0],
         [ 0.00025223,0],
         [ 0.000599317,0]]), np.array([[-2.52546,0],
         [-0.948863,0],
         [-0.412186,0],
         [-0.0896899,0],
         [0.0252699,0],
         [ 0.027219,0],
         [ 0.016151,0],
         [ 0.00750078,0],
         [  0.0023407,0],
         [0.000126069,0],
         [-0.0009829634,0],
         [-0.000882211,0],
         [-0.00111511,0],
         [-0.000365497,0],
         [-0.00103303,0]]), np.array([[-2.32485,0],
         [-0.828851,0],
         [-0.330234,0],
         [-0.0655762,0],
         [ 0.00455499,0],
         [0.00944265,0],
         [ 0.00503916,0],
         [ 0.00138588,0],
         [-0.000704993,0],
         [-0.00128207,0],
         [-0.00151531,0],
         [-0.00106773,0],
         [-0.00108,0],
         [-0.000365284,0],
         [-0.000883316,0]]), np.array([[-0.172272,0],
         [-0.252333,0],
         [-0.145669,0],
         [-0.0613297,0],
         [-0.0143016,0],
         [-0.00298425,0],
         [ 7.40548e-05,0],
         [ 0.000667491,0],
         [ 0.000606362,0],
         [ 0.000416257,0],
         [ 0.000265347,0],
         [ 0.000148069,0],
         [ 9.28048e-05,0],
         [ 3.56271e-05,0],
         [ 4.20592e-05,0]]), np.array([[ 3.0698,0],
         [-2.61454,0],
         [-0.541824,0],
         [ 0.151565,0],
         [ 0.0936031,0],
         [ 0.0457696,0],
         [0.0134952,0],
         [-0.00179418,0],
         [-0.00801073,0],
         [-0.00826738,0],
         [-0.0078117,0],
         [-0.00517042,0],
         [-0.00476933,0],
         [-0.00164482,0],
         [-0.00363309,0]]), np.array([[ 0.0994763,0],
         [-0.103695,0],
         [-0.00802014,0],
         [ 0.0117249,0],
         [ 0.00296628,0],
         [ 0.00080273,0],
         [-7.36281e-05,0],
         [-0.000342559,0],
         [-0.000385284,0],
         [-0.000312284,0],
         [-0.000253033,0],
         [-0.00015862,0],
         [-0.00013451,0],
         [-4.71665e-05,0],
         [-9.52027e-05,0]]), np.array([[ 0.471246,0],
         [-0.133993,0],
         [-0.031926,0],
         [-0.00360874,0],
         [-0.00553848,0],
         [-0.00370957,0],
         [-0.00238204,0],
         [-0.00148553,0],
         [-0.000926082,0],
         [-0.000557261,0],
         [-0.000350964,0],
         [-0.000195797,0],
         [-0.00013786,0],
         [-4.97421e-05,0],
         [-8.01491e-05,0]]), np.array([[-0.440552,0],
         [ 0.207406,0],
         [ 0.0864837,0],
         [ 0.0129719,0],
         [0.00118811,0],
         [-0.000816893,0],
         [-0.000589815,0],
         [-0.00018314,0],
         [ 8.95157e-05,0],
         [ 0.000170824,0],
         [ 0.000208261,0],
         [ 0.000147914,0],
         [ 0.000152199,0],
         [ 5.12068e-05,0],
         [ 0.000126199,0]])], 'HO2 + OH <=> H2O + O2': [np.array([[-0.391894,0],
         [0.0883868,0],
         [ 0.0980892,0],
         [ 0.0418015,0],
         [ 0.00994407,0],
         [ 0.00179395,0],
         [-0.000508979,0],
         [-0.000929414,0],
         [-0.000833894,0],
         [-0.000598773,0],
         [-0.000425728,0],
         [-0.000252007,0],
         [-0.000189625,0],
         [-6.85477e-05,0],
         [-0.000117339,0]]), np.array([[-1.74785,0],
         [-0.573832,0],
         [-0.0614013,0],
         [ 0.0185055,0],
         [ 0.0381509,0],
         [ 0.0305318,0],
         [0.0211072 ,0],
         [0.0136315 ,0],
         [ 0.0086409,0],
         [ 0.00523805,0],
         [ 0.00330108,0],
         [ 0.00184156,0],
         [ 0.00128788,0],
         [ 0.000466313,0],
         [0.000739679,0]]), np.array([[-1.45469,0],
         [-0.49767,0],
         [-0.0836142,0],
         [-0.00423854,0],
         [0.0129732 ,0],
         [0.0111135 ,0],
         [0.00724386 ,0],
         [0.0042474,0],
         [0.00234722,0],
         [0.00123948 ,0],
         [0.000631547 ,0],
         [0.000305339 ,0],
         [0.000143759 ,0],
         [5.73754e-05 ,0],
         [ 2.76603e-05,0]]), np.array([[-0.13615,0],
         [-0.169447,0],
         [-0.0919602,0],
         [-0.041309,0],
         [-0.0124285,0],
         [-0.00311021,0],
         [0.000174357 ,0],
         [ 0.00107748,0],
         [0.00117704 ,0],
         [0.000932631 ,0],
         [ 0.000733807,0],
         [0.000454426,0],
         [ 0.000375098,0],
         [0.000132517,0],
         [0.00025804,0]]), np.array([[0.351208,0],
         [-0.159316,0],
         [-0.114981,0],
         [-0.0365852,0],
         [-0.00143334,0],
         [0.00358292,0],
         [0.00324047,0],
         [ 0.00210435,0],
         [ 0.00118683,0],
         [0.000616998,0],
         [ 0.000286945,0],
         [0.000127624,0],
         [3.5656e-05,0],
         [1.76008e-05,0],
         [-2.30869e-05,0]]), np.array([[0.0639567,0],
         [-0.0270925,0],
         [-0.0187749,0],
         [-0.00548315,0],
         [-4.48604e-05,0],
         [0.000628959,0],
         [0.000505092,0],
         [0.000299064,0],
         [0.00014746,0],
         [6.38485e-05,0],
         [1.7017e-05,0],
         [2.14152e-06,0],
         [-1.01244e-05,0],
         [-2.70077e-06,0],
         [-1.51413e-05,0]]), np.array([[-0.422838,0],
         [0.0847614,0],
         [0.068645,0],
         [0.0396959,0],
         [0.0114195,0],
         [ 0.00354108,0],
         [0.000772429,0],
         [-0.000114316,0],
         [-0.000345142,0],
         [-0.000320127,0],
         [-0.000266607,0],
         [-0.000168467,0],
         [-0.000139589,0],
         [-4.97216e-05,0],
         [-9.5359e-05,0]]), np.array([[-0.0466193,0],
         [ 0.00499562,0],
         [ 0.0106031,0],
         [0.00562057,0],
         [0.00183642,0],
         [0.000573735,0],
         [8.94612e-05,0],
         [-7.06421e-05,0],
         [-0.00011099,0],
         [-9.67588e-05,0],
         [-8.07271e-05,0],
         [-5.11309e-05,0],
         [-4.35466e-05,0],
         [-1.53141e-05,0],
         [-3.0813e-05,0]]), np.array([[-0.186478,0],
         [-0.177408,0],
         [-0.0865274,0],
         [-0.0373494,0],
         [-0.0117016,0],
         [-0.00340818,0],
         [-0.000431228,0],
         [ 0.000484195,0],
         [ 0.000679238,0],
         [ 0.000572451,0],
         [ 0.000463884,0],
         [0.000290479,0],
         [ 0.000241979,0],
         [8.55709e-05,0],
         [0.000167414,0]]), np.array([[-0.148795,0],
         [-0.108753,0],
         [-0.0395802,0],
         [-0.0107795,0],
         [0.000422726,0],
         [0.00175679,0],
         [0.00140655,0],
         [0.000883189,0],
         [0.000494368,0],
         [0.000257855,0],
         [0.000123077,0],
         [5.61476e-05,0],
         [1.91185e-05,0],
         [8.63376e-06,0],
         [-5.27643e-06,0]]), np.array([[-0.00837597,0],
         [-0.0145598,0],
         [-0.0104529,0],
         [-0.00638864,0],
         [-0.00326243,0],
         [-0.00149822,0],
         [-0.000593584,0],
         [-0.000176743,0],
         [5.3114e-06,0],
         [5.80772e-05,0],
         [7.57048e-05,0],
         [5.43745e-05,0],
         [5.43968e-05,0],
         [1.86323e-05,0],
         [4.36952e-05,0]]), np.array([[-0.0437584,0],
         [-0.0165163,0],
         [-0.00719935,0],
         [-0.00261478,0],
         [-0.000759719,0],
         [-0.000257808,0],
         [-9.95723e-05,0],
         [-4.6266e-05,0],
         [-2.74131e-05,0],
         [-1.75486e-05,0],
         [-1.32878e-05,0],
         [-8.14891e-06,0],
         [-7.22477e-06,0],
         [-2.45129e-06,0],
         [-5.46175e-06,0]]), np.array([[ 0.222559,0],
         [-0.15142,0],
         [0.0196283,0],
         [0.0223217,0],
         [0.0117266,0],
         [0.00597953,0],
         [0.00280908,0],
         [0.00120402,0],
         [0.000420594,0],
         [9.56523e-05,0],
         [-5.25402e-05,0],
         [-6.43195e-05,0],
         [-9.19336e-05,0],
         [-3.02674e-05,0],
         [-8.84464e-05,0]])]}



MSI_st_instance_one = stMSIcheb.MSI_shocktube_optimization_chebyshev(cti_file,
                                                   .01,
                                                   1,
                                                   1,
                                                   working_directory,
                                                   files_to_include,                 
                                                   reaction_uncertainty_csv,rate_constant_target_value_data,
                                                   master_equation_reactions = master_equation_reactions,
                                                   chebyshev_sensitivities = cheb_sensitivity_dict,
                                                   master_reaction_equation_cti_name = master_reaction_equation_cti_name,
                                                   master_index = master_index,
                                                   master_equation_uncertainty_df = master_equation_uncertainty_df,
                                                   chebyshev_fit_nominal_parameters_dict = None)
MSI_st_instance_one.one_run_shock_tube_optimization()





S_matrix_original = MSI_st_instance_one.S_matrix
exp_dict_list_original = MSI_st_instance_one.experiment_dictonaries
original_covariance = MSI_st_instance_one.covarience
X_one_itteration = MSI_st_instance_one.X
MSI_st_instance_one.deltaXAsNsEas
Ydf_original = MSI_st_instance_one.Y_data_frame
Zdf_original = MSI_st_instance_one.z_data_frame
Ydf_original = MSI_st_instance_one.Y_data_frame


##need to fix this and return _s_matrix and y_matrix
MSI_st_instance_two = stMSIcheb.MSI_shocktube_optimization_chebyshev(cti_file,
                                                   .01,
                                                   1,
                                                   1,
                                                   working_directory,
                                                   files_to_include,                 
                                                   reaction_uncertainty_csv,rate_constant_target_value_data,
                                                   master_equation_reactions = master_equation_reactions,
                                                   chebyshev_sensitivities = cheb_sensitivity_dict,
                                                   master_reaction_equation_cti_name = master_reaction_equation_cti_name,
                                                   master_index = master_index,
                                                   master_equation_uncertainty_df = master_equation_uncertainty_df,
                                                   chebyshev_fit_nominal_parameters_dict = None)

delta_X_list = MSI_st_instance_two.multiple_shock_tube_runs(numer_of_iterations)

##ALL OF THIS STUFF CAN PROBABLY GO INTO SOME SORT OF CLASS
#
#
deltaXAsNsEas = MSI_st_instance_two.deltaXAsNsEas
physical_obervable_updates_list = MSI_st_instance_two.physical_obervable_updates_list
absorbance_observables_updates_list = MSI_st_instance_two.absorbance_coef_update_dict
Ydf = MSI_st_instance_two.Y_data_frame
Zdf = MSI_st_instance_two.z_data_frame
experimental_dicts = MSI_st_instance_two.experiment_dictonaries
z_matrix = MSI_st_instance_two.z_matrix
s_matrix = MSI_st_instance_two.s_matrix
y = MSI_st_instance_two.y_matrix
Y_matrix = MSI_st_instance_two.Y_matrix
S_matrix = MSI_st_instance_two.S_matrix

X = MSI_st_instance_two.X
Xdf = MSI_st_instance_two.X_data_frame
covarience = MSI_st_instance_two.covarience
exp_dict_list_optimized = MSI_st_instance_two.experiment_dictonaries
parsed_yaml_list = MSI_st_instance_two.list_of_parsed_yamls
sigma = MSI_st_instance_two.sigma
X = MSI_st_instance_two.X
delta_X = MSI_st_instance_two.delta_X
molecular_parameter_updates = MSI_st_instance_two.delta_x_molecular_params_by_reaction_dict
original_diag = np.diag(original_covariance)




#target_value_rate_constant_csv = 'MSI/data/test_data/FFCM1_custom_target_value_test.csv'
original_cti_file = MSI_st_instance_two.data_directory +'/'+ MSI_st_instance_two.cti_file_name

experiment_dict_uncertainty = MSI_st_instance_two.experiment_dict_uncertainty_original
target_value_csv = MSI_st_instance_two.data_directory +'/'+ MSI_st_instance_two.k_target_values_csv
if run_with_k_target_values == 'On' or run_with_k_target_values == 'on':
    k_target_value_S_matrix = MSI_st_instance_two.k_target_values_for_S
 #   print('poop')
else:
    k_target_value_S_matrix = None


##########################################################################################################################
#PLOTTING##
##########################################################################################################################
#csv_file_sigma = MSI_st_instance_two.data_directory +'/'+'sigma_for_uncertainty_weighted_sensitivity_FFCM1.csv'
#csv_file_sigma =  MSI_st_instance_two.data_directory +'/'+'sigma_for_uncertainty_weighted_sensitivity_glarborg.csv'
csv_file_sigma = ''

plotting_instance = plotter.Plotting(S_matrix,
                                     s_matrix,
                                     Y_matrix,
                                     y,
                                     z_matrix,
                                     X,
                                     sigma,
                                     covarience,
                                     original_covariance,
                                     S_matrix_original,
                                     exp_dict_list_optimized,
                                     exp_dict_list_original,
                                     parsed_yaml_list,
                                     Ydf,
                                     target_value_rate_constant_csv= MSI_st_instance_two.data_directory +'/'+ rate_constant_target_value_data_for_plotting ,
                                     target_value_rate_constant_csv_extra_values = MSI_st_instance_two.data_directory +'/'+rate_constant_target_value_data_extra,
                                     k_target_value_S_matrix =k_target_value_S_matrix,
                                     k_target_values=run_with_k_target_values,
                                     working_directory = working_directory,
                                     sigma_uncertainty_weighted_sensitivity_csv=csv_file_sigma,
                                     cheby_sensitivity_dict = cheb_sensitivity_dict,
                                     mapped_to_alpha_full_simulation=MSI_st_instance_two.mapped_to_alpha_full_simulation)

#csv_file_sigma = MSI_st_instance_two.data_directory +'/'+'sigma_for_uncertainty_weighted_sensitivity_updated.csv'
observable_counter_and_absorbance_wl,length_of_experimental_data = plotting_instance.lengths_of_experimental_data()
sigmas_optimized,test = plotting_instance.calculating_sigmas(S_matrix,covarience)
sigmas_original,test2 = plotting_instance.calculating_sigmas(S_matrix_original,original_covariance)
plotting_instance.plotting_observables(sigmas_original = sigmas_original,sigmas_optimized= sigmas_optimized)
diag = plotting_instance.getting_matrix_diag(covarience)

plotting_instance.plotting_rate_constants(optimized_cti_file=MSI_st_instance_two.new_cti_file,
                                original_cti_file=original_cti_file,
                                initial_temperature=250,
                                final_temperature=2500,
                                master_equation_reactions = master_equation_reactions)



sensitivity, top_sensitivity = plotting_instance.sort_top_uncertainty_weighted_sens()
obs = plotting_instance.plotting_uncertainty_weighted_sens()





