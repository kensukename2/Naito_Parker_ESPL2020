'''
Minnestota River, Jordan, MN application
4-regime analysis (1: 1935-1954, 2: 1955-1974, 3: 1975-1994, 4: 1995-2014)
Main script

Script by Kensuke Naito
'''


# Import libraries #########################################
############################################################
import os 
import sys
menow = os.getcwd()
sys.path.append(menow)
import numpy as np
import csv  
from _main_input_4regimes import *


# Main routine #############################################
############################################################
def main():
    # delete old files
    delete_old_files()

    # zero the time
    Time = 0.

    # regime 1935-1954 (spin-up) run
    print('Run: regime 0')
    Stop, Time, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
    Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt \
          = regime_0(Time)

    # regime 1935-1954, 1955-1974, 1975-1994, 1995-2014 runs
    for regime in range(1, 5):
        if Stop == True:
            break
        else:
            print('Run: regime ' + str(int(regime)))
            Stop, Time, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
            Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt \
                  = regime_1234(regime, Time, Qbf, Hbf, Bbf, Sc,
                                Hc, Hnc, Hpb, Etab, Etac, Etanc, Etapb,
                                Cero, Cdep, Vdep, Qt)


# Main functions ###########################################
############################################################
# sub routines for regime 1 (spin-up run)
def regime_0(Time):
    # setup empty arrays
    TIME = []; QBF = []; HBF = []; BBF = []; SC = []
    HC = []; HNC = []; HPB = []
    ETAB = []; ETAC = []; ETANC = []; ETAPB = []
    CERO = []; CDEP = []; VDEP = []; QT = []

    # get FDC
    Qi, pi = get_fdc(fdc_filename1)

    # set spacial array for initial condition
    usbf = (g * HbfI * ScI)**(1./2.)
    QbfI = Czc * usbf * HbfI * BbfI
    HpbI = phiN1 * HbfI
    CeroI = outerbank_erosion_rate(QbfI, BbfI, HbfI, HncI, ScI, Qi, pi)
    CdepI = innerbank_deposition_rate(QbfI, BbfI, HbfI, ScI, Qi, pi)
    VdepI = overbank_deposition_rate(QbfI, Qi, pi)
    QtI = bedmaterial_load(QbfI, BbfI, HbfI, ScI, Qi, pi)

    # set spacial array for initial condition
    Qbf = QbfI * np.ones(M); Bbf = BbfI * np.ones(M)
    Hbf = HbfI * np.ones(M); Sc = ScI * np.ones(M)
    Hnc = HncI * np.ones(M); Hc = HcI * np.ones(M)
    Hpb = HpbI * np.ones(M)
    Etab = Etab_d * np.ones(M)
    for j in range(M):
        Etab[j] += ScI * sinu * dxv * (M - j - 1.)
    Etanc = Etab + HncI; Etac = Etab + HbfI
    Etapb = Etab + HpbI
    Cero = CeroI * np.ones(M); Cdep = CdepI * np.ones(M)
    Vdep = VdepI * np.ones(M); Qt = QtI * np.ones(M)

    # time loop
    for k in range(Nloop_0):
        # spatial calculation
        Stop, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
        Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt \
            = main_calculation(Qbf, Bbf, Sc, Hnc, Hpb, Etac, Etanc, Etab, Etapb,
                               Cero, Cdep, Vdep, Qt, Qtfeed, Qi, pi)
        # update lists and time
        if Stop == True:
            break
        else:
            TIME.append(Time); QBF.append(np.mean(Qbf)); HBF.append(np.mean(Hbf))
            BBF.append(np.mean(Bbf)); SC.append(np.mean(Sc))
            HC.append(np.mean(Hc)); HNC.append(np.mean(Hnc)); HPB.append(np.mean(Hpb))
            ETAB.append(np.mean(Etab)); ETAC.append(np.mean(Etac))
            ETANC.append(np.mean(Etanc)); ETAPB.append(np.mean(Etapb))
            CERO.append(np.mean(Cero)); CDEP.append(np.mean(Cdep))
            VDEP.append(np.mean(Vdep)); QT.append(np.mean(Qt))
            Time += dt

    # print out
    if Stop != True:
        out_csv(out_filename0, np.array(TIME), np.array(QBF), np.array(HBF),
                np.array(BBF), np.array(SC),
                np.array(HC), np.array(HNC), np.array(HPB),
                np.array(ETAB), np.array(ETAC), np.array(ETANC), np.array(ETAPB),
                np.array(CERO), np.array(CDEP), np.array(VDEP), np.array(QT))

    return Stop, Time, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
           Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt


# sub routines for regime 1-4
def regime_1234(regime, Time, Qbf, Hbf, Bbf, Sc,
                Hc, Hnc, Hpb, Etab, Etac, Etanc, Etapb,
                Cero, Cdep, Vdep, Qt):

    # setup empty arrays
    TIME = []; QBF = []; HBF = []; BBF = []; SC = []
    HC = []; HNC = []; HPB = []
    ETAB = []; ETAC = []; ETANC = []; ETAPB = []
    CERO = []; CDEP = []; VDEP = []; QT = []

    # set file names and parameters for each regime
    if regime == 1:
        fdc_filename = fdc_filename1
        out_filename = out_filename1
    elif regime == 2:
        fdc_filename = fdc_filename2
        out_filename = out_filename2
    elif regime == 3:
        fdc_filename = fdc_filename3
        out_filename = out_filename3
    else:
        fdc_filename = fdc_filename4
        out_filename = out_filename4

    # get FDC
    Qi, pi = get_fdc(fdc_filename)

    # time loop
    for k in range(Nloop_1234):

        # spatial calculation
        Stop, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
        Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt \
            = main_calculation(Qbf, Bbf, Sc, Hnc, Hpb, Etac, Etanc, Etab, Etapb,
                               Cero, Cdep, Vdep, Qt, Qtfeed, Qi, pi)

        # update lists and time
        if Stop == True:
            break
        else:
            TIME.append(Time); QBF.append(np.mean(Qbf)); HBF.append(np.mean(Hbf))
            BBF.append(np.mean(Bbf)); SC.append(np.mean(Sc))
            HC.append(np.mean(Hc)); HNC.append(np.mean(Hnc)); HPB.append(np.mean(Hpb))
            ETAB.append(np.mean(Etab)); ETAC.append(np.mean(Etac))
            ETANC.append(np.mean(Etanc)); ETAPB.append(np.mean(Etapb))
            CERO.append(np.mean(Cero)); CDEP.append(np.mean(Cdep))
            VDEP.append(np.mean(Vdep)); QT.append(np.mean(Qt))
            Time += dt

    # print out
    if Stop != True:
        out_csv(out_filename, np.array(TIME), np.array(QBF), np.array(HBF),
                np.array(BBF), np.array(SC),
                np.array(HC), np.array(HNC), np.array(HPB),
                np.array(ETAB), np.array(ETAC), np.array(ETANC), np.array(ETAPB),
                np.array(CERO), np.array(CDEP), np.array(VDEP), np.array(QT))

    return Stop, Time, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
           Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt


# Sub function ############################################
###########################################################
# main calculation for channel geometry
def main_calculation(Qbf, Bbf, Sc, Hnc, Hpb, Etac, Etanc, Etab, Etapb,
                     Cero, Cdep, Vdep, Qt, Qtfeed, Qi, pi):
    # create arrays
    d_Etab, d_Etanc = np.zeros(M), np.zeros(M)
    d_Etac, d_Bbf = np.zeros(M), np.zeros(M)

    # calculate update
    for j in range(M):
        # Qt at adjacent nodes
        Qtmid = Qt[j]
        if j == 0:
            Qtup = Qtfeed
            Qtdw = Qt[j + 1]
        elif j == M - 1:
            Qtup = Qt[j - 1]
            Qtdw = Qtmid
        else:
            Qtup = Qt[j - 1]
            Qtdw = Qt[j + 1]

        d_Etab[j] = (Hnc[j] * Cero[j] - Hpb[j] * Cdep[j]) / Bbf[j] \
                    - (au * (Qtmid - Qtup) + (1. - au) * (Qtdw - Qtmid)) \
                    / (dxv * sinu * (1. - lpc) * Bbf[j])
        d_Etanc[j] = (Etapb[j] - Etanc[j]) * sinu * Cdep[j] / (phimb * Bbf[j] - sinu * Bbf[j])
        d_Etac[j] = Vdep[j] / (1. - lpc) + sinu * (Cdep[j] * (Etapb[j] - Etanc[j]) - Cero[j] * (Etac[j] - Etanc[j])) \
                    / (phimb * Bbf[j] - sinu * Bbf[j])
        d_Bbf[j] = Cero[j] - Cdep[j]

    # update geometry
    Bbf = Bbf + dt * d_Bbf
    Etab = Etab + dt * d_Etab
    Etanc = Etanc + dt * d_Etanc
    Etac = Etac + dt * d_Etac
    Hnc = Etanc - Etab
    Hc = Etac - Etanc
    Hbf = Hnc + Hc
    Hpb = phiN1 * Hbf
    Etapb = Etab + Hpb

    # update the rest
    for j in range(M):
        if j == 0:
            Sc[j] = (Etab[0] - Etab[1]) / (sinu * dxv)
        elif j == M - 1:
            Sc[j] = (Etab[-2] - Etab[-1]) / (sinu * dxv)
        else:
            Sc[j] = (Etab[j - 1] - Etab[j + 1]) / (2. * sinu * dxv)

        # detect NAN and negative geometry
        if Bbf[j] != Bbf[j] or Bbf[j] != Bbf[j] or Qt[j] != Qt[j]:
            Stop = True
            print('[!] Got NAN')
            break
        if Hc[j] <= 0 or Bbf[j] <= 0 or Sc[j] <= 0:
            Stop = True
            print('[!] Negative geometry')
            break
        else:
            Stop = False
            # bankfull discharge
            usbf = (g * Hbf[j] * Sc[j]) ** (1. / 2.)
            Qbf[j] = Czc * usbf * Hbf[j] * Bbf[j]
            # overbank deposition rate and channel migration rate
            Cero[j] = outerbank_erosion_rate(Qbf[j], Bbf[j], Hbf[j], Hnc[j], Sc[j], Qi, pi)  # [m/yr]
            Cdep[j] = innerbank_deposition_rate(Qbf[j], Bbf[j], Hbf[j], Sc[j], Qi, pi)  # [m/yr]
            Vdep[j] = overbank_deposition_rate(Qbf[j], Qi, pi)  # [m/yr]
            # bed material transport rate
            Qt[j] = bedmaterial_load(Qbf[j], Bbf[j], Hbf[j], Sc[j], Qi, pi)  # [m]

    return Stop, Qbf, Hbf, Bbf, Sc, Hc, Hnc, Hpb, \
           Etab, Etac, Etanc, Etapb, Cero, Cdep, Vdep, Qt


# function to delete old files
def delete_old_files():
    os.chdir(path_csv)
    for file in os.listdir(path_csv):
        if file.endswith('.csv'):
            os.remove(file)
    os.chdir(menow)

# function to dump results into csv file
def out_csv(out_filename, TIME, QBF, HBF, BBF, SC, HC, HNC, HPB,
            ETAB, ETAC, ETANC, ETAPB, CERO, CDEP, VDEP, QT):
    os.chdir(path_csv)
    Nrow = len(TIME)
    QT /= converter_Qt  # [m^3/s -> Mt/yr]
    with open(out_filename, 'wb') as ff:
        writer = csv.writer(ff)
        for kk in range(Nrow):
            # [0: TIME, 1:QBF, 2:HBF, 3:BBF, 4:SC,
            #  5:ETAB, 6:ETAC, 7:ETANC, 8:ETAPB,
            #  9:HC, 10:HNC, 11:HPB, 12:CERO, 13:CDEP, 14:VDEP, 15:QT]
            writer.writerow([
                '{:.1f}'.format(TIME[kk]), '{:.3f}'.format(QBF[kk]),
                '{:.3f}'.format(HBF[kk]), '{:.3f}'.format(BBF[kk]),
                '{:.8f}'.format(SC[kk]),
                '{:.3f}'.format(ETAB[kk]), '{:.3f}'.format(ETAC[kk]),
                '{:.3f}'.format(ETANC[kk]), '{:.3f}'.format(ETAPB[kk]),
                '{:.3f}'.format(HC[kk]), '{:.3f}'.format(HNC[kk]),
                '{:.3f}'.format(HPB[kk]),
                '{:.3f}'.format(CERO[kk]), '{:.3f}'.format(CDEP[kk]),
                '{:.10f}'.format(VDEP[kk]), '{:.3f}'.format(QT[kk])])
    os.chdir(menow)


# function to read FDC
def get_fdc(fdc_filename):
    # read file
    with open(fdc_filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        Q_raw = [];
        counter = 0
        for row in csv_reader:
            if str(row[3]).isdigit() and str(row[3]) > 0.:
                Q_raw.append(float(row[3]) * (12. ** 3.) * (2.54 ** 3.) / (100. ** 3.))  # [ft^3/s -> m^3/s]
            counter += 1
    Q_raw = np.array(Q_raw)
    number_of_data = len(Q_raw)
    Q_sort, F_sort = np.zeros(number_of_data), np.zeros(number_of_data)
    rank = np.linspace(1, number_of_data + 1, number_of_data, dtype='int')
    Q_sort[:] = -1.0 * np.sort(-Q_raw[:])
    F_sort[0] = float(rank[0]) / (float(number_of_data) + 1.)
    storage = 0
    for j in range(1, number_of_data):
        if Q_sort[j] < Q_sort[j - 1]:
            rank[j] = rank[j - 1] + 1 + storage
            storage = 0
        elif Q_sort[j] == Q_sort[j - 1]:
            rank[j] = rank[j - 1]
            storage += 1
        F_sort[j] = float(rank[j]) / (float(number_of_data) + 1.)
    # get pdf
    Qi, pi = np.zeros(Nfdc), np.zeros(Nfdc)
    nn, bins = np.histogram(Q_sort, Nfdc, density=True)
    for k in range(Nfdc):
        Qi[k] = 0.5 * (bins[k] + bins[k + 1])
        pi[k] = nn[k] * (bins[k + 1] - bins[k])
    return Qi, pi


# calculate overbank deposition rate
def overbank_deposition_rate(Qbf, Qi, pi):
    v_dep = 0.
    for i in range(Nfdc):
        # below bankfull flow
        if Qi[i] > Qbf:
            v_dep += phie * vfall * year_in_sec * Cfm * pi[i]
    return v_dep


# calculate outerbank erosion rate
def outerbank_erosion_rate(Qbf, Bbf, Hbf, Hnc, Sc, Qi, pi):
    if Hnc < 0.:
        Hnc = 0.
    tausbf = Hbf * Sc / (Rr * Ds)
    c_ero = 0.
    for i in range(Nfdc):
        # below bankfull flow
        if Qi[i] <= Qbf:
            H = (Qi[i] ** 2. / (Czc ** 2. * g * Sc * Bbf ** 2.)) ** (1. / 3.)
            tausi = H * Sc / (Rr * Ds)
            c_ero += 0.5 * phiN3 * Hnc * (tausi / tausbf) * pi[i]
        # above bankfull flow
        else:
            c_ero += 0.5 * phiN3 * Hnc * 1. * pi[i]
    return c_ero


# calculate innerbank deposition rate
def innerbank_deposition_rate(Qbf, Bbf, Hbf, Sc, Qi, pi):
    tausbf = Hbf * Sc / (Rr * Ds)
    c_dep = 0.
    for i in range(Nfdc):
        # below bankfull flow
        if Qi[i] <= Qbf:
            H = (Qi[i] ** 2. / (Czc ** 2. * g * Sc * Bbf ** 2.)) ** (1. / 3.)
            tausi = H * Sc / (Rr * Ds)
            c_dep += 0.5 * c_vegR / (1. - lp) * (1. - tausi / tausbf) * pi[i]
    return c_dep


# calculate bed material load
def bedmaterial_load(Qbf, Bbf, Hbf, Sc, Qi, pi):
    tausbf = Hbf * Sc / (Rr * Ds)
    Einstein_no = 0.
    for i in range(Nfdc):
        # below overbank flow
        if Qi[i] <= Qbf:
            H = (Qi[i] ** 2. / (Czc ** 2. * g * Sc * Bbf ** 2.)) ** (1. / 3.)
            tausi = H * Sc / (Rr * Ds)
            Einstein_no += Bbf * aEH * 0.05 * Czc ** 2. * tausi ** (5. / 2.) * pi[i]
        # above overbank flow
        else:
            Einstein_no += Bbf * aEH * 0.05 * Czc ** 2. * tausbf ** (5. / 2.) * pi[i]
    Qt = Einstein_no * (g * Rr * Ds ** 3.) ** (1. / 2.) * year_in_sec  # [m^3/yr]
    return Qt


# run the script ###########################################
############################################################
if __name__ == '__main__':
    main()
