'''
Minnestota River, Jordan, MN application
Normal distribution FDC analysis
Input

Script by Kensuke Naito
'''


# Import libraries #########################################
############################################################
import numpy as np


# Parameters ###############################################
############################################################
# set paths
path_csv = 'path to fdc csv' 

# setup filenames
out_filename0 = '_results_scenario0.csv'
out_filename1 = '_results_scenario1.csv'
out_filename2 = '_results_scenario2.csv'

'''
Between 1934-1954
Q_mean = 91.5 cms
Q_total = 668622 cms total of 20 yrs
Q_total = 33431 cms total of 1 yr
Q_max = 1781 cms
'''

# set parameters for each scenario
# universal base discharge
Qbase = 5.     # [m^3/s] base discharge
Qscaler = 33000
Tmu = 50
# base scenario (scenario 0)
Tsigma0 = 10
# scenario 1
Tsigma1 = 15
# scenario 2
Tsigma2 = 5
# parameters for FDC and hydrograph
Nfdc = 365    #[1] # of bin in FDC (100)
Nevent = 365     # total number of event

# space / time step
dxv = 1000.     # [m] vary step length
M = 10     # number of grids
dt = 0.1     # [yr]  time step
Nloop0 = 10001   # year of run for scenario 0 (spin-up)
Nloop12 = 20001   # year of run for scenario 1,2

# initial channel geometry
BbfI = 150.
HncI = 1.5
HcI = 1.5
HbfI = HncI + HcI
ScI = 0.0002164

# auxullary parameter
Czc = 13.          # [1] Coef of Chezy coef of channel
aEH = 0.3         # [1] Engelund-Hansen correction coefficient
au = 0.5          # [1]

# downstream boundary
Etab_d = 10.

# basic setup
Ds = 0.3 / 1000.         # [mm(m)] D50 of bed material
# Bmb = 2000.     # [m] meander belt width
sinu = 2.          # [1] channel sinuosity
Qtfeed = 0.22    # [Mt/yr] bed material supply

# Parameter for floodplain material rating curve
# afm = 3e-7
# bfm = 0.59
Cfm = 0.00015     # [1] Characteristic volumetric concentration of floodplain material

 # parameter for deposition / erosion calculation
A = 4.81          # [1] dimensionless coeff (Johannesson & Parker) (5)
phimb = 10.          # [1] ratio of meandering belt width to bankfull width (10)
phic = 1.          # [1] ratio of slump block size to cohesive layer thickness (0.5)
TsbR = 0.05          # [yr] characteristic life time of slump block (0.5)
c_vegR = 9.          # [m/yr] characteristic time for the bank encroachment
lpc = 0.3          # [1] porosity of floodplain material (0.3)
lp = 0.3          # [1] porosity of bank material (0.3)
phie = 0.06        # [1] fractional effectivity of floodplain deposition (0.5)
Dfm = 0.04 / 1000.          # [mm(m)] floodplain material size (0.05)

# global constants
g = 9.81          #[m/s^2]
Rr = 1.65          #[1]
day_in_sec = 60. * 60. * 24.          #[s]
year_in_sec = day_in_sec * 365.25          #[s]
converter_Qt = 1e6 / (Rr + 1.)

# unit conversion and further calculation
Xv = np.ones(M)     # River valley coordinates for plot
for j in range(M):
    Xv[j] *= j * dxv / 1000.
Qtfeed *= converter_Qt     # [Mt/yr -> m^3/yr]
Rep = np.sqrt(g * Rr * Dfm**3.) / 1e-6          # [1] particle Reynolds number
b1 = 2.89139447769084; b2 = 0.95296; b3 = 0.0568346711984055 
b4 = 0.00289204602084475; b5 = 0.00024464688411386
Rf = np.exp(-b1 + b2 * np.log(Rep) - b3 * np.log(Rep)**2. - b4 * np.log(Rep)**3. + b5 * np.log(Rep)**4.)
vfall = Rf * np.sqrt(g * Rr * Dfm)          # [m/s] fall velocity of floodplain material
phiN1 = (sinu - 1.) * A / phimb     # [1]
phiN3 = phic / ((1. - lpc) * TsbR)     # [1/yr]
