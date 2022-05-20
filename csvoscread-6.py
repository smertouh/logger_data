# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 17:58:19 2019

@author: sotnikov
"""

import csv
import numpy as np
import glob, os
from scipy.interpolate import interp1d
import sys

def column(matrix, i):
    return [row[i] for row in matrix]



#xxxx.insert(0,x)

namefolder=[]
j=2

for file in glob.glob("**/chany0.txt", recursive=True):
    namefolder.append(file.rsplit("chany0",1)[0])
for snf in namefolder:
    try:
        names = ["time"]
        coefs = [0.001]
        zeroval = [0]
        xxxx = []
        xxxxt = []
        xxx = []
        xxxt = []
        if snf == 'chany0.txt':
            namefolder[0] = ''
        namefile = []
        for file in glob.glob(snf + "chany*.txt"):
            namefile.append(file.split(".")[0])

        for s in namefile:
            name = str(s) + '.txt'
            # name2 = "param" + name
            name2 = name.rsplit("chany", 1)[0] + "paramchany" + name.rsplit("chany", 1)[1]
            # открываем файл осц и конфигурации
            with open(name) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=';')
                a = [row for row in csv_reader]
            with open(name2) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=';')
                a1 = [row for row in csv_reader]
            # читаем название и коэффициент
            namei = "noname"
            newnamei=["noname","Utot"]#2022-03-16
            coefi = 1
            zerovali = 0
            for row in a1:

                try:
                    s = str(row).split("=")[0].rsplit("'", 2)[1]
                    if s == 'label':
                        namei = (str(row).split("=")[1].rsplit("'", 2)[0])
                    if s == 'display_unit':
                        coefi = (float(str(row).split("=")[1].rsplit("'", 2)[0]))
                    if s == 'zero_start':
                        try:
                            zerovali = (float(str(row).split("=")[1].rsplit("'", 2)[0]))
                            zerovali = 1
                            if 'T0'in namei:
                                zerovali=0
                            elif namei in newnamei: #2022-03-16
                                zerovali = 0 #2022-03-16
                        except:
                            zerovali = 0
                except:
                    print("no label")
            names.append(namei)
            coefs.append(coefi)
            zeroval.append(zerovali)
            xxx = []
            xx = []
            xx = (column(a, 1))
            xx = str(xx)
            xx = (xx.split('\''))
            # вычитываем столбец значений
            for ss in xx:
                while ss is not float:
                    try:
                        sss = float(ss.replace(',', '.'))
                        xxx.append(sss)
                        break
                    except ValueError:
                        break
            # если получилось зачитать запоминаем
            if xxx != []:
                xxxx.insert(j, xxx)
            # вычитываем стобец времен
            xxxt = []
            xxt = []
            xxt = (column(a, 0))
            xxt = str(xxt)
            xxt = (xxt.split('\''))
            for ss in xxt:
                while ss is not float:
                    try:
                        sss = float(ss.replace(',', '.'))
                        xxxt.append(sss)
                        break
                    except ValueError:
                        break
                # k=k+1
                # if s.isdigit():
                # xxx.append(s)
            # если получилось зачитать запоминаем
            if xxxt != []:
                xxxxt.insert(j, xxxt)
            j = j + 1

        xxx = []
        xx = []
        # вычитываем время
        with open(snf + "chany0.txt") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=';')
            a = [row for row in csv_reader]
        xxxt = []
        xxt = []
        xxt = (column(a, 0))
        xxt = str(xxt)
        xxt = (xxt.split('\''))
        for ss in xxt:
            while ss is not float:
                try:
                    sss = float(ss.replace(',', '.'))
                    ssss=(round(float(sss),0))#2021-06-03 new
                    if ssss in xxxt: #2021-06-03 new
                        print("повтор1!")
                        pass #2021-06-03 new
                    else: #2021-06-03 new
                        #print(ssss)
                        xxxt.append(ssss)#2021-06-03 sss
                    break
                except ValueError:
                    break
            # k=k+1
            # if s.isdigit():
            # xxx.append(s)
        for z in xxxxt:
            for zz in z:
                ssss = (round(float(zz), 0))
                if ssss in xxxt:
                    #print("повтор2!")
                    pass
                else:
                    xxxt.append(ssss)
                    #print(zz)

        xxxt = sorted(xxxt)

        if xxxt != []:
            xxxx.insert(0, xxxt)

        coef = coefs
        names.append("Ib")
        Ib = np.size(names) - 1
        names.append("PRF")
        PRF = np.size(names) - 1

        i = 0
        # ищем канал перехвата и второго тока
        Iac = 1
        IAG = 1
        RF_UG1 = 1
        RF_UA1 = 1
        Cath1_C = 1
        S_C1A = 1
        for s in names:

            if s == 'Iac_2':
                Iac = i
            elif s == 'IAG':
                IAG = i
            elif s == 'RF_UG1':
                RF_UG1 = i
            elif s == 'RF_UA1':
                RF_UA1 = i
            elif s == 'Cath1_C':
                Cath1_C = i
            elif s == 'S_C1(A)':
                S_C1A = i

            i+= 1

        # z=np.zeros((2245,1))

        xxxx.append(xxxt)
        xxxxt.append(xxxt)
        zeroval.append(0)
        coef.append(1)
        xxxx.append(xxxt)
        xxxxt.append(xxxt)
        zeroval.append(0)
        coef.append(1)
        name = names

        # интерполируем

        n = 0
        while (n < np.size(names) - 1):
            f = interp1d(xxxxt[n], xxxx[n + 1], kind='nearest',
                         fill_value="extrapolate")  # fill_value=(0,0), bounds_error=False
            xxxx[n + 1] = f(xxxx[0])
            n = n + 1

        xxxxx = np.transpose(xxxx)
        xxxxxt = np.transpose(xxxxt)

        # выставляем 0 и умножаем на коэффициент
        n = 0
        while (n < np.size(names) - 1):
            xxxxx[:, n] = (xxxxx[:, n] * coef[n])
            if zeroval[n]:
                xxxxx[:, n] = xxxxx[:, n] - xxxxx[1, n]
            n = n + 1
        # щитаем Ib
        for i in np.arange(0, np.size(xxxt)):
            xxxxx[i, Ib] = xxxxx[i, Iac] - xxxxx[i, IAG]

        # щитаем PRF
        for i in np.arange(0, np.size(xxxt)):
            if  abs(xxxxx[i,RF_UG1])>77:
                alfa=np.arccos((-480 + 403) / xxxxx[i,RF_UG1])
                xxxxx[i, PRF] = (alfa - np.sin(alfa) * np.cos(alfa)) / (np.sin(alfa) - alfa * np.cos(alfa)) * 0.5 * \
                                xxxxx[i, RF_UA1] * (xxxxx[i, Cath1_C] - xxxxx[i, S_C1A])
            else:
                alfa=0
                xxxxx[i, PRF] = 0

        # записываем в файл
        if snf!='':
            wn=snf.rsplit("\\ADC_")[0]+".sotnikov"
        else:
            wn="osc.sotnikov"
        with open((snf + wn), 'w') as csv_file2:
            print(wn)
            wr = csv.writer(csv_file2)  # ,dialect='excel')
            wr.writerow(name)
            wr.writerows(xxxxx)
    except:
        print("папка "+snf+" проклята")




















