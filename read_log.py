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






def mass_value(m, atr, ttime):
    for t in m['time']:
        if t==ttime:
            n=m['time'].index(t)
            return m[atr][n], t
    return '', ''




def read_log_dict():
    for file in glob.glob("*-*-*.sotnikov", recursive=False):
        namefolder.append(file.rsplit("chany0", 1)[0])
    l0=100
    for name0 in namefolder:
        if len(name0)<l0:
            l0=len(name0)
            name = name0
    with open(name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        a = [row for row in csv_reader]

    Main_dict = {}
    atr_list = []
    for s in a:
        if len(s) > 2:
            if not atr_list:
                for ss in s:
                    atr_list.append(ss)
                    Main_dict.update({ss:[]})
            else:
                j=0
                for atr in atr_list:
                    Main_dict[atr].append(s[j])
                    j+=1


    return Main_dict


Main_dict=read_log_dict()

print('end')

rest=(mass_value(Main_dict,"Iex","15:53:31")[0])
resm=mass_value(Main_dict,"Iex","15:53:31")[1]
print(rest, resm)
