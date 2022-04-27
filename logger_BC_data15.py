
import csv
import numpy as np
import glob, os
from scipy.interpolate import interp1d
import time
import math
import sys
try:
    import read_quadera #2020-06-09
except:
    pass
try:
    import read_log  # 2020-06-10
except:
    pass

dT=0.15 #полуширина метки
#пытаемся прочитать данные с призмы
try:
    MQD=read_quadera.read_quadera_dict()
except:
    pass
try:
    AD=read_log.read_log_dict()
except:
    pass

"""
rest=time.strftime('%d.%m.%Y %H:%M:%S',time.localtime(read_quadera.mass_value(MQD,"18","2.05.2021 18:43:22",20)[0]))
resm=read_quadera.mass_value(MQD,"18","2.05.2021 18:43:22",20)[1]
print(rest, resm)
"""
print(1)
def column(matrix, i):
    return [row[i] for row in matrix]
def colminus(matrix, i,j):

    return [float(row[i])-float(row[j]) for row in matrix]
def matrixmin(matrix, i):
    d=[]
    for s in matrix[:,i]:
        d.append(float(s))
    return np.min(d)

def matrixmax(matrix, i):
    d=[]
    for s in matrix[:,i]:
        d.append(float(s))
    return np.max(d)

def DeltaT(T1, T2):
    return (T2.tm_hour - T1.tm_hour) * 3600 + (T2.tm_min - T1.tm_min) * 60 + T2.tm_sec - T1.tm_sec

def colzerminus(matrix):
    size1=np.size(matrix)
    n=100
    delta=0#(float(matrix[-1])-float(matrix[n]))/(size1-n)
    s=0#float(matrix[n])
    for i in range(size1):
        matrix[i]=(float(matrix[i])-delta*(i-n)-s)
        if (float(matrix[i]) < 0.03):
            if float(matrix[i]) > -0.03:
                matrix[i]=0
    return 1

def colinterpmat(matrix, ii):

    matf=[]
    timevec=np.arange(-180,0,0.1)
    for _i in range(np.size(ii)+1):
        matf.append(timevec)
    vect=[]
    for d in matrix[:, 0]:
        vect.append(float(d))
    j=1
    for i in ii:
        vec=[]
        for d in matrix[:,i]:
            vec.append(float(d))
        f = interp1d(vect, vec, kind='nearest',
                     fill_value="extrapolate")  # fill_value=(0,0), bounds_error=False
        matf[j] = f(timevec)
        j+=1
    return matf


class signaloscs:
    """класс для чтения осциллограм"""

    def __init__(self, nm=[[0, 1],[2,4]], Time='00:00:00', ColumnNames={}, Columntype={},Longtime='01.01.2021 00:00:00'):
        if ColumnNames:
            pass
        else:

            self.Columntype = {
                'IPG': 'ADC',
                'UPG': 'ADC',
                'Iex': 'ADC',
                'Uex': 'ADC',
                'Iac_2': 'ADC',
                'IAG': 'ADC',
                'Utot': 'miniFC',
                'UHV': 'slowADC',
                'T01': 'slowADC',
                'T02': 'slowADC',
                'T03': 'slowADC',
                'T04': 'slowADC',
                'Flow': 'constvalue',
                'IHVo': 'slowADC',
                'Ip': 'ADC',
                'Q21': 'Bcal',
                'Q43': 'Bcal',
                'Qtot': 'Bcal',
                'Ib': 'calc',
                'Ie': 'Ie',
                'transp1': 'transp',
                'transp2': 'transp',
                'Ibmax':'Ibmax',
                'RF_UG1':'ADC',
                "RF_UA1":'ADC',
                "Cath1_C":'ADC_no_zero',
                'S_C1(A)':'ADC',
                'PRF':'PRF',
                'Ipmax':'Ipmax',
                'IFC':'ADC',
                'ifc1':'miniFC',
                'ifc2':'miniFC',
                'ifc3': 'miniFC',
                'ifc4': 'miniFC',
                'ifc5': 'miniFC',
                'quadera2':'Quadera',
                'quadera18':'Quadera',
                'Ret':'Log',
                'Tcst': 'Log',
                'OvenUpU': 'Log',
                'OvenUpI': 'Log',
                'PCsup':'PCs',
                'FCa':'Log',
                'FCb': 'Log',
                'FCl': 'Log',
                'PH2Pa':'Log',
                'I15kV': 'Log',
                'U15kV': 'Log',
                'I110kV': 'Log',
                'U110kV': 'Log',
                'ProtStatus': 'ProtStatus',
                'mark_time': 'mark_time',
                'Tcsb': 'Log',
                'SED1_B':'miniFC',
                'SED1_R': 'miniFC',
                'SED1_T': 'miniFC',
                'SED1_L': 'miniFC'

            }




            self.ColumnNames = self.Columntype.keys()

        self.ColumnNamesDict = {}

        for s in self.ColumnNames:
            self.ColumnNamesDict.update({s: ""})
        """    
            if not s in self.ColumnNamesDict.keys():
                name=s
                self.ColumnNamesDict.update({name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})
        """
        _i = 0
        if not nm:
            nm = [[1, 2], [3, 4]]
        for row in nm[0]:
            # row=el
            if (row != []):
                if row in self.ColumnNamesDict.keys():
                    self.ColumnNamesDict.update({row: signalosc(nm, _i, self.Columntype[row], signame=row)})



            _i += 1
        self.BC = []

        self.Time = Time
        self.Longtime=Longtime

        self.calcvalues()



        # print(self.ColumnNamesDict)  # self.ColumnNamesDict[row])

    def __del__(self):
        class_name = self.__class__.__name__
        #print('{} уничтожен'.format(class_name))
    def BCcal(self,name):
        try:
            if not self.BC:
                self.BC = self.cmdt(self.ColumnNamesDict['T01'], self.ColumnNamesDict['T02'],
                                    self.ColumnNamesDict['T03'], self.ColumnNamesDict['T04'],
                                    self.ColumnNamesDict['Flow'])
                print(self.BC)
            self.ColumnNamesDict[name].setvalue(self.BC[name])
        except:
            pass
            # print(self.BC)
    def Ibmax(self,name):
        try:
            z = (self.ColumnNamesDict['Iac_2']) - (self.ColumnNamesDict['IAG'])
            zmax = np.max(z.nm)
        except:
            zmax = 0
        self.ColumnNamesDict[name].setvalue(zmax)

    def Ipmax(self,name):
        try:
            z = (self.ColumnNamesDict['Ip'])
            #z.set_t1(6.75)
            #zmax = z.readvalue(t2=6.75)
            zampl=np.max(z.nm)-np.min(z.nm)
        except:
            zampl = 0
        self.ColumnNamesDict[name].setvalue(zampl)

    def PRF(self,name):
        try:
            alfa=np.arccos((-480+403)/self.ColumnNamesDict['RF_UG1'].readvalue())
            URF=self.ColumnNamesDict['RF_UA1'].readsinglevalue()
            C_C=self.ColumnNamesDict['Cath1_C'].readvalue()
            S_C1=self.ColumnNamesDict['S_C1(A)'].readvalue()
            pRF=(alfa-np.sin(alfa)*np.cos(alfa))/(np.sin(alfa)-alfa*np.cos(alfa))*0.5*URF*(C_C-S_C1)
            self.ColumnNamesDict[name].setvalue(pRF)
        except:
            pass
    #2021-06-07
    def Ie(self,name):
        print('Ie')
        try:
            I_ex=self.ColumnNamesDict['Iex'].readvalue()
            I_ac = self.ColumnNamesDict['Iac_2'].readvalue()
            I_ag = self.ColumnNamesDict['IAG'].readvalue()
            Ielectron=I_ex-I_ac+I_ag
            self.ColumnNamesDict[name].setvalue(Ielectron)
        except:
            pass
    def PCs(self,name):
        try:
            if name=='PCsup':
                U=float(self.ColumnNamesDict['OvenUpU'].readvalue())
                I =float(self.ColumnNamesDict['OvenUpI'].readvalue())
            else:
                pass
            P=U*I
            self.ColumnNamesDict[name].setvalue(P)
        except:
            pass

    def ProtStatus(self,name):#2022-01-11
        try:
            p1=(self.ColumnNamesDict['I15kV'].readvalue())
            p2=(self.ColumnNamesDict['U15kV'].readvalue())
            p3 = (self.ColumnNamesDict['I110kV'].readvalue())
            p4 = (self.ColumnNamesDict['U110kV'].readvalue())
            if p1=='True'or p2=='True' or p3=='True' or p4 == 'True':
                P='True'
            else:
                P='False'
            #print(P)
            self.ColumnNamesDict[name].setvalue(P)
        except:
            pass

    def mark_time(self,name):#2022-01-17
        try:
            P = float(sys.argv[1])
            if (self.ColumnNamesDict['ProtStatus'].readvalue())=='False':

                self.ColumnNamesDict[name].setvalue(P)
            else:
                Times=np.arange(P-2.4,P+2.4,0.3)
                #values=[]
                vMax=0
                time_without_BD = []
                for t in Times:
                    if round((self.ColumnNamesDict['Uex'].readvalue(np.nan, t)), 1)==vMax:
                        #values.append(vMax)
                        time_without_BD.append(P-t)
                    elif round((self.ColumnNamesDict['Uex'].readvalue(np.nan, t)), 1)>vMax:
                        vMax=round((self.ColumnNamesDict['Uex'].readvalue(np.nan, t)), 1)
                        #values=[vMax]
                        time_without_BD=[P-t]
                #values=[v0,v1, v2, v3, v4,v5]
                _i=time_without_BD.index(min(time_without_BD))
                P=round(P-time_without_BD[_i],1)
                self.ColumnNamesDict[name].setvalue(round(P,1))
            #print(P)
            return P
        except:
            pass

    def Quadera(self,name):
        try:
            quad=read_quadera.mass_value(MQD, name.rsplit('quadera')[-1], self.Longtime, 20)[1]
            self.ColumnNamesDict[name].setvalue(quad)
        except:
            pass
    def R_Log(self,name):
        try:
            quad=read_log.mass_value(AD,name,self.Time)[0]
            self.ColumnNamesDict[name].setvalue(quad)
        except:
            pass
    def calcvalues(self):


        """
        needednames=['mark_time', 'Uex','ProtStatus','I15kV','U15kV','I110kV','U110kV']
        for name in needednames:
            self.ColumnNamesDict.update({name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})
        self.markT = self.mark_time("mark_time")
        print(self.markT)
        """

        #print('начинаю читать')
        for name in self.ColumnNames:
            if self.ColumnNamesDict[name] == '':


                if self.Columntype[name] == 'Log': #2021-06-10
                    self.ColumnNamesDict.update({name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})
                    self.R_Log(name)
                elif self.Columntype[name] == 'ProtStatus': #2022-01-11
                    self.ColumnNamesDict.update({name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})
                    self.ProtStatus(name)
                elif self.Columntype[name] == 'mark_time':#2022-01-17
                    self.ColumnNamesDict.update({name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})
                    self.markT=self.mark_time(name)
                    #print(self.markT)




        #print(self.markT)
        for name in self.ColumnNames:
            if True:#self.ColumnNamesDict[name] == '':
                try:
                    ##self.ColumnNamesDict.update({name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})
                    if type(self.ColumnNamesDict[name]) is str:
                        self.ColumnNamesDict.update(
                            {name: signalosc(np.nan, np.nan, self.Columntype[name], signame=name)})

                    self.ColumnNamesDict[name].setT2(self.markT)
                    self.ColumnNamesDict[name].setType(self.Columntype[name])
                    self.ColumnNamesDict[name].setSigname(name)
                    # print(self.markT)
                    if self.ColumnNamesDict[name].type == 'Bcal':
                        self.BCcal(name)
                    elif self.ColumnNamesDict[name].type == 'Ibmax':
                        self.Ibmax(name)
                    elif self.ColumnNamesDict[name].type == 'PRF':
                        self.PRF(name)
                    elif self.ColumnNamesDict[name].type == 'Ipmax':
                        self.Ipmax(name)
                    elif self.ColumnNamesDict[name].type == 'Ie':  # 2021-06-07
                        self.Ie(name)
                    elif self.ColumnNamesDict[name].type == 'Quadera':  # 2021-06-09
                        self.Quadera(name)
                    elif self.ColumnNamesDict[name].type == 'Log':  # 2021-06-10
                        self.R_Log(name)
                    elif self.ColumnNamesDict[name].type == 'PCs':  # 2021-06-10
                        self.PCs(name)
                    elif self.ColumnNamesDict[name].type == 'ProtStatus':  # 2022-01-11
                        self.ProtStatus(name)
                    elif self.ColumnNamesDict[name].type == 'mark_time':  # 2022-01-17
                        self.markT = self.mark_time(name)
                        #print(self.markT)
                    else:
                        pass

                        #self.ColumnNamesDict.update({name: '0'})
                except:
                    pass
        #print("Прочёл!")

        """
        if self.markT!=float(sys.argv[1]):#2022-01-17
            for name in self.ColumnNames:  # 2022-01-17
                try:
                    self.ColumnNamesDict[name].setT2(self.markT)
                except:
                    print("проблемы с ", name)
        """






    def cmdt(self, T1, T2, T3, T4, Flow):
        nm = [T1.nmt, T1.nm, T2.nm, T3.nm, T4.nm, Flow.nm]
        matrix3 = self.colinterpmat(nm, [1, 2, 3, 4, 5])
        colzerminus(matrix3[1])
        colzerminus(matrix3[2])
        colzerminus(matrix3[3])
        colzerminus(matrix3[4])
        # считаем дельту
        dT21 = matrix3[2] - matrix3[1]
        dT32 = matrix3[3] - matrix3[2]
        dT41 = matrix3[4] - matrix3[1]
        dT43 = matrix3[4] - matrix3[3]
        Flowgpm = matrix3[5]
        q21 = 0
        q32 = 0
        q41 = 0
        q43 = 0

        wn2 = 'asd.osc'
        na = np.zeros([np.size(matrix3[0]) + 1, 5])
        with open((wn2), 'w') as csv_file2:

            name = ["time", "q21", "q32", "q41", "q43"]
            wr = csv.writer(csv_file2)  # ,dialect='excel')
            wr.writerow(name)
            _j = 0
            for s in matrix3[0]:
                na[_j][0] = s
                na[_j][1] = matrix3[1][_j]
                na[_j][2] = matrix3[2][_j]
                na[_j][3] = matrix3[3][_j]
                na[_j][4] = matrix3[4][_j]
                _j += 1
            wr.writerows(na)


        for i in np.arange(1, (np.size(dT21) - 2)):
            try:
                q21 = q21 + dT21[i] * float(Flowgpm[i]) * (0.1 * 4183 * 3.79 / 60 / 1000)
                q32 = q32 + dT32[i] * float(Flowgpm[i]) * (0.1 * 4183 * 3.79 / 60 / 1000)
                q41 = q41 + dT41[i] * float(Flowgpm[i]) * (0.1 * 4183 * 3.79 / 60 / 1000)
                q43 = q43 + dT43[i] * float(Flowgpm[i]) * (0.1 * 4183 * 3.79 / 60 / 1000)
            except:
                print(np.size(dT43), " ", np.size(dT43))
        if q21>0:
            if q43>0:
                qtot=q21+q43
            else:
                qtot=q21
        else:
            if q43>0:
                qtot=q43
            else:
                qtot=''
        return {'Q21': q21, 'Q43': q43, 'Qtot': qtot}

    def colinterpmat(self, nm, ii):

        # ii=[T01,T02,T03,T04,Flow]
        matf = []
        timevec = np.arange(-180, 0, 0.1)
        for _i in range(np.size(ii) + 1):
            matf.append(timevec)
        vect = nm[0]

        j = 1
        for i in ii:
            vec = nm[i]
            f = interp1d(vect, vec, kind='nearest',
                         fill_value="extrapolate")  # fill_value=(0,0), bounds_error=False
            matf[j] = f(timevec)
            j += 1
        return matf

    def logdata(self):
        # z=[]
        z = {}
        for name in self.ColumnNames:
            # z.append(str(self.ColumnNamesDict[name]))
            z.update({name: str(self.ColumnNamesDict[name])})
        return {self.Time: z}


class signalosc:
    def columnfl(self, matrix, i):
        return [float(row[i]) for row in matrix]

    def __init__(self, nm, signal, type='ADC', t1=4.2, t2=np.nan, signame=""):   #2021-05-31 t2=5.5 #2022-01-17 t2=np.nan
        if t2 is np.nan:
            try:
                t2 = float(sys.argv[1])  # new 2021-06-07
                # print(sys.argv[1])
            except:
                print("Нет времени")
                t2 = 5.0  # new 2021-06-09
                pass

        if signame=='IPG' or signame=='UPG' or signame=='RF_UG1' or  signame ==  "RF_UA1" or signame == "Cath1_C" or signame == 'S_C1(A)':
            t1=0
        if not signal is np.nan:
            matrix2 = np.delete(nm.copy(), 0, axis=0)
            self.nmt = self.columnfl(matrix2, 0)
            self.nm = self.columnfl(matrix2, signal)
        else:
            self.value = np.nan
        self.t1 = t1
        self.t2 = t2
        self.slowT=np.nan
        self.type = type
        self.signame = signame
    def __add__(self,a):
        _j = 0
        na = np.zeros([np.size(self.nmt) + 1, 2])
        for v in self.nm:
            na[_j, 1] = (float(v) + float(a.nm[_j]))
            na[_j, 0] = self.nmt[_j]
            _j += 1
        nm = na
        b = signalosc(nm, 1, type=self.type, t1=self.t1, t2=self.t2, signame=self.signame + "+" + a.signame)
        return b
    def __sub__(self,a):
        _j=0
        na=np.zeros([np.size(self.nmt)+1,2])
        for v in self.nm:
            na[_j,1]=(float(v)-float(a.nm[_j]))
            na[_j,0]=self.nmt[_j]
            _j += 1
        nm=na
        b=signalosc(nm,1,type=self.type,t1=self.t1,t2=self.t2,signame=self.signame+"-"+a.signame)
        return b
    def __mul__(self,a):
        _j=0
        na=np.zeros([np.size(self.nmt)+1,2])
        for v in self.nm:
            na[_j,1]=(float(v)*float(a.nm[_j]))
            na[_j,0]=self.nmt[_j]
            _j += 1
        nm=na
        b=signalosc(nm,1,type=self.type,t1=self.t1,t2=self.t2,signame=self.signame+"*"+a.signame)
        return b
    def __div__(self,a):
        _j=0
        na=np.zeros([np.size(self.nmt)+1,2])
        for v in self.nm:
            if a.nm[_j]!=0:
                na[_j,1]=(float(v)/float(a.nm[_j]))
            else:
                na[_j,1]=9999999
            na[_j,0]=self.nmt[_j]
            _j += 1
        nm=na
        b=signalosc(nm,1,type=self.type,t1=self.t1,t2=self.t2,signame=self.signame+"//"+a.signame)
        return b
    def __str__(self):
        # print("asf")

        try:
            s = str(self.readvalue(self.t1, self.t2))
        except:
            s= " " #2022-01-19
        return s

    def __repr__(self):
        # print("asf")
        s = str(self.readvalue(self.t1, self.t2))
        return s

    def readsinglevalue(self, t=np.nan):
        if t is np.nan:
            t=self.t2
        SS = 0
        j = 0
        jj = 0
        for s in self.nmt:
            if t + dT > float(s) and float(s) > t - dT:
                SS += float(self.nm[jj])
                j += 1
            jj += 1
        if j != 0:
            SS = SS / j
        # print(SS)
        return SS

    def setvalue(self, value):
        self.value = value

    def setT2(self, value):
        self.t2=value
    def setType(self, value):
        self.type=value
    def setSigname(self, value):
        self.signame=value


    def readvalue(self, t1=np.nan, t2=np.nan):
        if t1 is np.nan:
            t1=self.t1
        if t2 is np.nan:
            t2=self.t2
        # print("2124")
        if self.type == 'ADC' or self.type == 'calc':
            d = (self.readsinglevalue(t2) - self.readsinglevalue(t1))
        elif self.type== 'ADC_no_zero':
            d = self.readsinglevalue(t2)
        elif self.type == 'slowADC':
            d = np.max(self.nm) - np.min(self.nm)
        elif self.type == 'Ibmax':
            d = self.value
        elif self.type == 'constvalue':
            d = np.average(self.nm)
        elif self.type == 'Bcal':
            d = self.value
        elif self.type == 'PRF':
            d = self.value
        elif self.type == 'Ipmax':
            d = self.value
        elif self.type == 'Ie': #2021-06-07
            d = self.value
        elif self.type == 'PCs': #2021-06-07
            d = self.value
        elif self.type == 'Log': #2021-06-07
            d = self.value
        elif self.type == 'Quadera': #2021-06-09
            d = self.value
        elif self.type == 'miniFC':
            """ 2021-05-31
            if math.isnan (self.slowT):
                i=0
                for f in reversed(self.nm):

                    if float(f)**2 >0.01:
                        self.slowT=self.nmt[i]
                        break
                    i-=1
            """
            try:
                self.slowT=float(sys.argv[1]) #new 2021-05-31 #new 2021-06-07
                #print(sys.argv[1])
            except:
                self.slowT=5.0
                print('Нет времени')
            try:
                d = self.readsinglevalue(self.slowT) #2021-05-31  d = self.readsinglevalue(self.slowT-1.25)
            except:
                d=0
        elif self.type == 'ProtStatus': #2022-11-01
            d = self.value

        else:
            try:
                d = self.value
                return d
            except:
                print('DEBUG!', self.signame)
                return 0
        return d
    def set_t1(self,t1):
        self.t1=t1
    def set_t2(self,t2):
        self.t2=t2

class LOG:
    def __init__(self,osc):
        self.Log = {}
        self.TimeLog = []
        self.columnnames=['Time']+list(osc.ColumnNames)
        self.newcolumnnames=['Time','IPG','UPG','Iex','Uex','Iac_2','IAG','Utot','UHV','Flow','IHVo','Ip','Ipmax','Q21','Q43','Qtot','Ib','Ie','Ibmax','transp1','transp1_1','transp1_2','transp2','PRF','IFC','ifc1','ifc2','ifc3','ifc4','ifc5','quadera2','quadera18','Ret','Tcst','PCsup','ProtStatus', 'PH2Pa','mark_time','Tcsb','SED1_B','SED1_R','SED1_L','SED1_T']
    def update(self, osc, Time):
        self.Log.update(osc.logdata())
        self.TimeLog.append(Time)
    def buildtable(self,names=''):
        nm=[]
        if not names:
            names=self.columnnames
        #for s in self.TimeLog:
        for s in self.T:
            na = [s]
            z = self.Log[s]
            for s1 in names:
                if s1!='Time':
                    try:
                        na.append(z[s1])
                    except:
                        na.append('')
            nm.append(na)
        return nm

    def BC3min(self):
        T=[]
        for s in self.TimeLog:
            T.append(time.strptime(s,'%H:%M:%S'))
        T = sorted(T)
        #t0=time.mktime(T[0])
        t0=T[0]
        dT=[]
        for t in T:
            deltat = DeltaT(t0, t)
            dT.append(deltat)

            if (175 < deltat) and (deltat <185):
                self.BCshift(time.strftime('%H:%M:%S',t0),time.strftime('%H:%M:%S',t))
                t0 = time.strftime('%H:%M:%S', t0)
            else:
                #print('УДАЛЯЮ СТРОКУ ',time.strftime('%H:%M:%S',t0),' была неправильная пауза')
                t0 = time.strftime('%H:%M:%S', t0)
                self.Log[t0].update({'Q21':''})
                self.Log[t0].update({'Q43': ''})
                self.Log[t0].update({'Qtot': ''})

            #print ("получилось", t0, " Q21= ",self.Log[t0]["Q21"]," Q43= ",self.Log[t0]["Q43"]," Qtot= ", self.Log[t0]["Qtot"])

            # считаю транспортировку на 10м
            self.transcalc(t,t0)
            t0=t
        T = [time.strftime('%H:%M:%S',x) for x in T]
        #dT = [time.strftime('%H:%M:%S', x) for x in dT]
        #print (T)
        #print (dT)
        self.T=T
    def transcalc(self,t,t0):

        try:
            self.Log[t0].update({"transp2": (float(self.Log[t0]['Qtot']) / float(self.Log[t0]['Ib']) / 2.5 / (
                    float(self.Log[t0]['Utot']) + float(self.Log[t0]['UHV'])))})
        except:
            pass
        t0 = time.strftime('%H:%M:%S', t)
        try:
            self.Log[t0].update({"transp1": (float(self.Log[t0]['IHVo']) / float(self.Log[t0]['Ibmax']))})
        except:
            pass
        try:
            self.Log[t0].update({"transp1_1": (float(self.Log[t0]['Ip']) / float(self.Log[t0]['Ib']))})
        except:
            pass
        try:
            self.Log[t0].update({"transp1_2": (float(self.Log[t0]['Ipmax']) / float(self.Log[t0]['Ibmax']))})
        except:
            pass
    def minuscheck(self,t0,t,Q):
        if not '-' in self.Log[t][Q]:
            #print(self.Log[t][Q]," ", Q,">0 двигаю значения", t," в ",t0)
            try:
                self.Log[t0].update({Q: float(self.Log[t][Q])})
            except:
                pass

        else:
            #self.Log[t0].update({'Qtot': float(self.Log[t0]['Qtot'])-float(self.Log[t0][Q])})
            #print(self.Log[t][Q]," ", Q,"<0 удаляю значения в ", t0)
            self.Log[t0].update({Q: ''})

    def BCshift(self,t0,t):

        self.minuscheck(t0,t,'Q21')
        self.minuscheck(t0, t, 'Q43')
        self.minuscheck(t0, t, 'Qtot')
        #print("получилось", t0, " Q21= ", self.Log[t0]["Q21"], " Q43= ", self.Log[t0]["Q43"], " Qtot= ",              self.Log[t0]["Qtot"])
        #print("удаляю", t, " Q21= ", self.Log[t]["Q21"], " Q43= ", self.Log[t]["Q43"], " Qtot= ",              self.Log[t]["Qtot"])
        self.Log[t].update({'Q21': ''})
        self.Log[t].update({'Q43': ''})
        self.Log[t].update({'Qtot': ''})

    def save(self,wn):
        with open((wn), 'w') as csv_file2:
            try:
                self.BC3min()
            except:
                pass
            nm1 = self.buildtable(self.newcolumnnames)
            wr = csv.writer(csv_file2)  # ,dialect='excel')
            if not self.newcolumnnames:
                wr.writerow(self.columnnames)
            else:
                wr.writerow(self.newcolumnnames)
            wr.writerows(nm1)

            #self.BCshift()


class LOGFC(LOG):
    def __init__(self,osc):
        self.Log = {}
        self.TimeLog = []
        self.columnnames=['Time']+list(osc.ColumnNames)
        #self.newcolumnnames=['Time','IFC','ifc1','ifc2','ifc3','ifc4','ifc5','a','b','L'] #2022-01-11
        self.newcolumnnames = ['Time', 'IFC', 'ifc1', 'ifc2', 'ifc3', 'ifc4', 'ifc5', 'FCa', 'FCb', 'FCl','Uex','ProtStatus','Ib','Utot']  # 2022-01-11
    def selfT(self):
        T=[]
        for s in self.TimeLog:
            T.append(time.strptime(s,'%H:%M:%S'))
        T = sorted(T)
        T = [time.strftime('%H:%M:%S',x) for x in T]
        self.T=T

    def BC3min(self):
        self.selfT()
        self.FCposABL()
        pass
    def FCposABL(self):
        with open('abl.txt') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=';')
            a = [row for row in csv_reader]
        print(a)
        for s in a:
            try:
                Tfc = time.strptime(s[0], '%H:%M')
            except:
                Tfc = time.strptime(s[0], '%H:%M:%S')
            for t in self.T:
                if time.strptime(t, '%H:%M:%S')>=Tfc:
                    print(s[1],' ',s[2],' ',s[3])
                    #t0 = time.strftime('%H:%M:%S', t)
                    self.Log[t].update({'FCa': s[1]})#a
                    self.Log[t].update({'FCb': s[2]})#b
                    self.Log[t].update({'FCl': s[3]})#L

                    pass
    def save(self,wn,Uex='1'):
        #print('safasfgasgaga')
        Uex1=6.2
        Uex2=6.8

        try: #2022-01-14
            with open((wn), 'w') as csv_file2:
                try:
                    self.BC3min()
                except:
                    pass
                nm1 = self.buildtable(self.newcolumnnames)
                wr = csv.writer(csv_file2)  # ,dialect='excel')
                if not self.newcolumnnames:
                    wr.writerow(self.columnnames)
                else:
                    wr.writerow(self.newcolumnnames)
                wr.writerows(nm1)
        except:
            pass





        if not Uex:
            with open((wn), 'w') as csv_file2:
                try:
                    self.BC3min()
                except:
                    pass
                nm1 = self.buildtable(self.newcolumnnames)
                wr = csv.writer(csv_file2)  # ,dialect='excel')
                if not self.newcolumnnames:
                    wr.writerow(self.columnnames)
                else:
                    wr.writerow(self.newcolumnnames)
                wr.writerows(nm1)

        else:
            Uex = [5.5, 5.8, 6.2, 6.5, 6.8, 7.2, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12]
            for Uex0 in Uex:
                Uex1=Uex0-0.1
                Uex2=Uex0+0.1

                if os.path.exists('FCtable_Uex'):
                    pass
                else:
                    os.mkdir('FCtable_Uex')

                wn='FCtable_Uex\\'+str(Uex0)+'.osc'
                with open((wn), 'w') as csv_file2:
                    try:
                        self.BC3min()
                    except:
                        pass
                    nm1 = self.buildtable(self.newcolumnnames)
                    wr = csv.writer(csv_file2)  # ,dialect='excel')
                    if not self.newcolumnnames:
                        wr.writerow(self.columnnames)
                    else:
                        wr.writerow(self.newcolumnnames)
                    # wr.writerows(nm1)
                    for row in nm1:
                        if Uex2 > float(row[10]) > Uex1:
                            wr.writerow(row)




osc=signaloscs()
Log=LOG(osc)
LogFC=LOGFC(osc)
namefile=[]
rowT=[]
rowT2=[]
for file in glob.glob("**/*.sotnikov", recursive=True):
    namefile.append(file)
for namefiles in namefile:
    print(namefiles)
    xxxx = []
    name=str(namefiles)
    with open(name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        a=[row for row in csv_reader]
    s=str(a[0]).split(',')

    line_count = 0
    nm = []
    i=0
    for row in a:
        na = []
        try:
            for s in row:
                ss = str(s).split(',')
                for ss2 in ss:
                    sss = str(ss2).rsplit("'")
                    ssss = str(sss).split("'")[1]

                    if (ssss != "="):
                        if (ssss != ""):
                            if ssss != []:
                                na.append(str(ssss))
        except:
            break
        if na!=[]:
            nm.append(na)
            #print(na)



    try:
        # Считываем импульс
        nfs = name.rsplit('.', 1)[0].rsplit('\\', 1)[-1].rsplit('_', 1)[1]
        Time = nfs[0] + nfs[1] + ":" + nfs[2] + nfs[3] + ":" + nfs[4] + nfs[5]
        newString = name.rsplit('\\', 1)[1].rsplit('.sotnikov')[0]
        Longtime=newString[8]+newString[9]+'.'+newString[5]+newString[6]+'.'+newString[0]+newString[1]+newString[2]+newString[3]+" "+newString[-6]+newString[-5]+":"+newString[-4]+newString[-3]+":"+newString[-2]+newString[-1]
        osc = signaloscs(nm, Time,Longtime=Longtime)
        try:
            osc.cmdt(osc.ColumnNamesDict['T01'], osc.ColumnNamesDict['T02'], osc.ColumnNamesDict['T03'],
                     osc.ColumnNamesDict['T04'], osc.ColumnNamesDict['Flow'])
        except:
            print('нет измерений калориметром')
        # Записываем ЛОГ
        Log.update(osc, Time)
        LogFC.update(osc, Time)
        # print(osc.ColumnNamesDict)
        print(osc.logdata())
        # print(Log.Log)
        #print(Log.TimeLog)
    except:
        pass

if os.path.exists('table'):
    pass
else:
    os.mkdir("table")



try:
    wn = "table\\"+newString.rsplit("_")[0]+"-caldata-"+sys.argv[1]+".osc"#"caldata"+sys.argv[1]+".osc"
    Log.save(wn)
    print('измерения мощности на калориметр ок')
except:
    print('нет данных калориметра')

try:
    wn2= "table\\"+"IFC"+sys.argv[1]+".osc"
    LogFC.save(wn2)
    print('измерение профиля пучка ОК')
except:
    print('ЦФ не двигали')