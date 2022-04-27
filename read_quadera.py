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
import time


namefolder=[]

def get_nearest_value(value, llist):
    list_of_difs = [10000000000000000000]
    xx=0
    for x in llist:
        if x>value:
            break
        if abs(value - x)>list_of_difs[-1]:
            break
        else:
            list_of_difs.append(abs(value - x))
            xx=x

    return xx#list_of_difs[-1]



def get_nearest_value2(value, llist):
    xx=llist[0]

    if value>llist[-1]:
        xx=llist[-1]
        ind=-1
        return xx, ind
    for x in llist:
        if x>value:
            if xx in llist:
                ind=llist.index(xx)
            else:
                ind = 0
            break
        xx=x

    return xx, ind



def mass_value(m, mass, ttime,dt=0):

    t = time.mktime(time.strptime(ttime, '%d.%m.%Y %H:%M:%S'))-dt
    [ntime,n]= get_nearest_value2(t, m[mass+'Time'])
    ncurrent= m[mass][n]
    return ntime, ncurrent


def read_quadera_dict():
    for file in glob.glob("*.asc", recursive=True):
        namefolder.append(file.rsplit("chany0", 1)[0])
    name = namefolder[0]
    with open(name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        a = [row for row in csv_reader]

    Main_dict = {}
    mass_list = []
    for s in a:
        if len(s) > 2:
            if not mass_list:
                for ss in s:
                    if ss != '':
                        mass_list.append(ss)
            else:
                pass
                """
                try:
                    t=time.strptime(s[0],'%d.%m.%Y %H:%M:%S')
                    print(t)
                except:
                    passtime.localtime
                """
            if len(s) % 3 == 0:
                j = 0
                for mass in mass_list:
                    if j+2<len(s)-1:
                        if ((mass + 'Time') in Main_dict):
                            Main_dict[(mass + 'Time')].append(time.mktime(time.strptime(s[j + 0], '%d.%m.%Y %H:%M:%S')))
                            Main_dict[mass].append(s[j + 2])
                        else:
                            Main_dict.update({mass + 'Time': []})
                            Main_dict.update({mass: []})

                    j += 3
    return Main_dict

Main_dict=read_quadera_dict()

print('end')

rest=time.strftime('%d.%m.%Y %H:%M:%S',time.localtime(mass_value(Main_dict,"18","2.05.2021 18:43:22",20)[0]))
resm=mass_value(Main_dict,"18","2.05.2021 18:43:22",20)[1]
print(rest, resm)