'''
Minnestota River, Jordan, MN application
4-regime analysis (1: 1935-1954, 2: 1955-1974, 3: 1975-1994, 4: 1995-2014)
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
fdc_filename0 = '_fdc_minnesotariver_jordan_1935_1954.csv'
fdc_filename1 = '_fdc_minnesotariver_jordan_1935_1954.csv'
fdc_filename2 = '_fdc_minnesotariver_jordan_1955_1974.csv'
fdc_filename3 = '_fdc_minnesotariver_jordan_1975_1994.csv'
fdc_filename4 = '_fdc_minnesotariver_jordan_1995_2014.csv'
out_filename0 = '_results_regime0.csv'
out_filename1 = '_results_regime1.csv'
out_filename2 = '_results_regime2.csv'
out_filename3 = '_results_regime3.csv'
out_filename4 = '_results_regime4.csv'

# space / time step
dxv = 2000.     # [m] vary step length
M = 20     # number of grids
dt = 0.01     # [yr]  time step
Nloop_0 = 100001   # year of run for regime 1 (spin-up)
Nloop_1234 = 2001   # year of run for regime 2
print('spin-up run [yr] = ' + str(int(dt * Nloop_0)))
print('regime 2-4 run [yr] = ' + str(int(dt * Nloop_1234)))

# initial channel geometry
BbfI = 100.
HncI = 2.
HcI = 2.
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

# parameters for FDC
Nfdc = 100    #[1] # of bin in FDC (100)

# Parameter for floodplain material rating curve
Cfm = 0.00015     # [1] Characteristic volumetric concentration of floodplain material

 # parameter for deposition / erosion calculation
A = 4.81          # [1] dimensionless coeff (Johannesson & Parker) (5)
phimb = 10.          # [1] ratio of meandering belt width to bankfull width (10)
phic = 1.          # [1] ratio of slump block size to cohesive layer thickness (0.5)
TsbR = 0.1          # [yr] characteristic life time of slump block (0.5)
c_vegR = 8.          # [m/yr] characteristic time for the bank encroachment
lpc = 0.3          # [1] porosity of floodplain material (0.3)
lp = 0.3          # [1] porosity of bank material (0.3)
phie = 0.2        # [1] fractional effectivity of floodplain deposition (0.5)
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
