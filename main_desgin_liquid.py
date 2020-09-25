# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 07:22:22 2018

@author: TJ
"""

import numpy as np
import pandas as pd
import os
import json
from datetime import datetime
from tqdm import tqdm
from scipy import optimize
from cea_post import Read_datset

class Cal_excond:
    """ Class for calculate combustion parameter to determin experimental condition
    
    Attributes
    -------
    cond: dictionary of float and string
        Dictionary of experimental and calculating condition
        "d": float, [m] port diamter
        "N": int, [-] the number of port
        "Df": float, [m] outer diameter of fuel cylinder
        "eta": float, [-] efficiency of characteristics exhaust velocity
        "rho_f": float, [kg/m3] fuel density
        "Rm": float, [J/kg/K] gas constant
        "Tox": float, [K] oxidizer gas temperature
        "mu": float, [Pa-s] oxidizer dynamic viscocity
        "a": float, [-] fuel filling rate
        "cea_path": string, folder path of CEA csv data-base 
    model_const: dictionary of float and string
        Dictionary of model constants which are used in Vf experimental formula
        "mode": string, mode selection. 
                "C1C2" mode uses the following experimental formula: Vf = (C1/Vox+C2)P^n, 
                "EXP" mode uses the following experimental formula: Vf = beta*P^n*Vox^m,
                "PROP" mode uses the following experimental formula: Vf = alpha*P^n
        "n": float, pressure exponent
        "C1": float, optional, SI-unit experimental constant
        "C2": float, optional, SI-unit experimental constant
        "beta": float, optional, SI-unit experimental constant
        "m": float, optional, exponent of oxidizer port velocity
        "alpha": float, optional, SI-unit experimental constant
    """
    def __init__(self, cond, const, **kwargs):
        """constructer

        Parameter
        -------------
        cond: dictionary of float and string
            Dictionary of experimental and calculating condition
            "d": float, [m] port diamter
            "N": int, [-] the number of port
            "Df": float, [m] outer diameter of fuel cylinder
            "eta": float, [-] efficiency of characteristics exhaust velocity
            "rho_f": float, [kg/m3] fuel density
            "Rm": float, [J/kg/K] gas constant
            "Tox": float, [K] oxidizer gas temperature
            "mu": float, [Pa-s] oxidizer dynamic viscocity
            "a": float, [-] fuel filling rate
            "cea_path": string, folder path of CEA csv data-base 
        model_const: dictionary of float and string
            Dictionary of model constants which are used in Vf experimental formula
            "mode": string, mode selection. 
                    "C1C2" mode uses the following experimental formula: Vf = (C1/Vox+C2)P^n, 
                    "EXP" mode uses the following experimental formula: Vf = beta*P^n*Vox^m,
                    "PROP" mode uses the following experimental formula: Vf = alpha*P^n
            "n": float, pressure exponent
            "C1": float, optional, SI-unit experimental constant
            "C2": float, optional, SI-unit experimental constant
            "beta": float, optional, SI-unit experimental constant
            "m": float, optional, exponent of oxidizer port velocity
            "alpha": float, optional, SI-unit experimental constant
        """
        self.cond = cond
        self.model_const = const
        self.d = self.cond["d"]
        self.N = self.cond["N"]
        self.Df = cond["Df"]
        self.eta = cond["eta"]
        self.rho_f = cond["rho_f"]
        self.Rm = cond["Rm"]
        self.mu = cond["mu"]
        self.Tox = cond["Tox"]
        self.cea_path = cond["cea_path"]
        self.C = cond["C"]
        self.be = cond["be"]
        self.d_o = cond["d_o"]
        self.rho_ox = cond["rho_ox"]
        self.a = 1 - self.N*np.power(self.d, 2)/np.power(self.Df, 2)
        if self.model_const["mode"] == "C1C2":
            self.C1 = self.model_const["C1"]
            self.C2 = self.model_const["C2"]
            self.n = self.model_const["n"]
        elif self.model_const["mode"] == "EXP":
            self.m = self.model_const["m"]
            self.beta = self.model_const["beta"]
            self.n = self.model_const["n"]
        elif self.model_const["mode"] == "PROP":
            self.alpha = self.model_const["alpha"]
            self.n = self.model_const["n"]
        elif self.model_const["mode"] == "LIQUID":
            self.C11 = self.model_const["C11"]
            self.C12 = self.model_const["C12"]
            self.C21 = self.model_const["C21"]
            self.C22 = self.model_const["C11"]
            self.m1 = self.model_const["m1"]
            self.m2 = self.model_const["m2"]
            self.n = self.model_const["n"]
        else:
            print("This program does't have such a mode: \"{}\".".format(self.model_const["mode"]))
            print("Please assign a mode from the following: \"C1C2\", \"EXP\" or \"PROP\". ")

    def _iterat_func_Pc_(self, Pc, mox, Dt):
        """ function which output error of Pc for Pc iteration
        
        Parameters
        ----------
        Pc : float
            [Pa] chamber pressure
        mox : float
            [kg/s] oxidzer mass flow rate
        Dt : float
            [m] nozzle throat diameter

        Return
        ----------
        error: float
            relative error of cmaber pressure between assumed and calculated value
        """
        Vox = func_Vox(Pc, mox, self.Rm, self.Tox, self.Df, self.a)
        if self.model_const["mode"] == "C1C2":
            Vf = func_Vf(Vox, Pc, self.C1, self.C2, n=self.n)
        elif self.model_const["mode"] == "EXP":
            # Vf = func_Vf_exp(Vox, Pc, self.beta, self.m, n=self.n)
            pass
        else:
            # Vf = func_Vf_prop(Vox, Pc, self.alpha, n=self.n)
            pass
        of = func_of(Vox, Vf, self.rho_f, self.a, Pc, self.Rm, self.Tox)
        func_cstr = gen_func_cstr(self.cea_path)
        Pc_cal = func_Pc_cal(of, Pc, mox, func_cstr, Dt, self.eta)
        diff = Pc_cal - Pc
        error = diff/Pc_cal
        return error

    def get_excond(self, mox, Dt, Pc_init=1.0e+6):
        """function for calculate combustion parameter from assigned experimental condition
        
        Parameters
        ----------
        mox : float
            [kg/s] oxidizer port velocity
        Dt : float
            [m] nozzle throat diameter
        Pc_init : float, optional
            [Pa] initial guess of chamber pressure, by default 1.0e+6
        
        Returns
        -------
        dic_excond: dictionary of float
            dictionary of calculated combustion parameter
            "mox": float, [kg/s] oxidizer mass flow rate. This is assingned value
            "Dt": float, [m] nozzle throat diameter. This parameter is assigned value
            "Pc": float, [Pa] chamber pressure
            "Vox": float, [m/s] oxidzer port velocity
            "Vf": float, [m/s] axial fuel regression rate
            "of": float, [-] oxidizer to fuel mass ratio
            "cstr_ex": float, [m/s] experimental characteristic exhaust velocity
            "Re": float, [-] Reynoldz number of oxidizer at each port
            "ustr_lam": float, [m/s] friction velocity when the flow is laminar
            "ustr_turb": float, [m/s] friction velocity when the flow is turbulence
        """
        func_cstr = gen_func_cstr(self.cea_path)
        try:
            Pc = optimize.newton(self._iterat_func_Pc_, Pc_init, maxiter=10, tol=1.0e-3, args=(mox, Dt))
        except (RuntimeError, RuntimeWarning):
            try:
                Pc = optimize.brentq(self._iterat_func_Pc_, 0.01e+6, 50e+6, maxiter=100, xtol=1.0e-3, args=(mox, Dt), full_output=False)
            except (RuntimeError, ValueError, RuntimeWarning):
                Pc = np.nan

        Vox = func_Vox(Pc, mox, self.Rm, self.Tox, self.Df, self.a)
        if self.model_const["mode"] == "C1C2":
            Vf = func_Vf(Vox, Pc, self.C1, self.C2, n=self.n)
        elif self.model_const["mode"] == "EXP":
            # Vf = func_Vf_exp(Vox, Pc, self.beta, self.m, n=self.n)
            pass
        else:
            # Vf = func_Vf_prop(Vox, Pc, self.alpha, n=self.n)
            pass
        of = func_of(Vox, Vf, self.rho_f, self.a, Pc, self.Rm, self.Tox)
        cstr = self.eta* func_cstr(of, Pc)
        Re = func_Re(Pc, Vox, self.Tox, self.Rm, self.d, self.mu)
        ustr_lam = func_ustr_lam(Pc, Vox, self.Tox, self.Rm, self.d, self.mu)
        ustr_turb = func_ustr_turb(Pc, Vox, self.Tox, self.Rm, self.d, self.mu)
        dic_excond = {"mox": mox,
                      "Dt": Dt,
                      "Pc": Pc,
                      "Vox": Vox,
                      "Vf": Vf,
                      "of": of,
                      "cstr_ex": cstr,
                      "Re": Re,
                      "ustr_lam": ustr_lam,
                      "ustr_turb": ustr_turb
                      }
        return dic_excond

    def _iterat_func_Pc_liquid_(self, Pc, Pup, Dt):
        """ function which output error of Pc for Pc iteration
        
        Parameters
        ----------
        Pc : float
            [Pa] chamber pressure
        Pup : float
            [Pa] reservoir pressure. In this calculation, Pup is the same as orifice up pressure.
        Dt : float
            [m] nozzle throat diameter

        Return
        ----------
        error: float
            relative error of cmaber pressure between assumed and calculated value
        """
        mox = func_mox_liquid(Pup, Pc, self.C, self.be, self.d_o, self.rho_ox)
        Vox = func_Vox_liquid(mox, self.rho_ox, self.Df, self.a)
        Vf = func_Vf_liquid(Vox, Pc, self.d, self.C11, self.C12, self.C21, self.C22, self.m1, self.m2, n=self.n)
        of = func_of_liquid(mox, Vf, self.rho_f, self.a, self.Df)
        func_cstr = gen_func_cstr(self.cea_path)
        Pc_cal = func_Pc_cal(of, Pc, mox, func_cstr, Dt, self.eta)
        diff = Pc_cal - Pc
        # error = diff/Pc_cal
        error = diff/Pc
        return error


    def get_excond_liquid(self, Pup, Dt, Pc_init=0.1e+6):
        """function for calculate combustion parameter from assigned experimental condition
        
        Parameters
        ----------
        mox : float
            [kg/s] oxidizer port velocity
        Dt : float
            [m] nozzle throat diameter
        Pc_init : float, optional
            [Pa] initial guess of chamber pressure, by default 1.0e+6
        
        Returns
        -------
        dic_excond: dictionary of float
            dictionary of calculated combustion parameter
            "mox": float, [kg/s] oxidizer mass flow rate. This is assingned value
            "Dt": float, [m] nozzle throat diameter. This parameter is assigned value
            "Pc": float, [Pa] chamber pressure
            "Vox": float, [m/s] oxidzer port velocity
            "Vf": float, [m/s] axial fuel regression rate
            "of": float, [-] oxidizer to fuel mass ratio
            "cstr_ex": float, [m/s] experimental characteristic exhaust velocity
            "Re": float, [-] Reynoldz number of oxidizer at each port
            "ustr_lam": float, [m/s] friction velocity when the flow is laminar
            "ustr_turb": float, [m/s] friction velocity when the flow is turbulence
        """
        func_cstr = gen_func_cstr(self.cea_path)
        try:
            # Pc = optimize.newton(self._iterat_func_Pc_liquid_, Pc_init, maxiter=100, tol=1.0e-3, args=(Pup, Dt))
            Pc = optimize.newton(self._iterat_func_Pc_liquid_, Pup*0.99, maxiter=100, tol=1.0e-3, args=(Pup, Dt))
        except (RuntimeError, RuntimeWarning):
            # Pc = optimize.brentq(self._iterat_func_Pc_liquid_, 0.1e+6, Pup*0.99 , maxiter=100, xtol=1.0e-3, args=(Pup, Dt), full_output=False)
            try:
                Pc = optimize.brentq(self._iterat_func_Pc_liquid_, 0.1e+6, Pup*0.99 , maxiter=100, xtol=1.0e-3, args=(Pup, Dt), full_output=False)
            except (RuntimeError, ValueError, RuntimeWarning):
                Pc = np.nan
        mox = func_mox_liquid(Pup, Pc, self.C, self.be, self.d_o, self.rho_ox)
        Vox = func_Vox_liquid(mox, self.rho_ox, self.Df, self.a)
        Vf = func_Vf_liquid(Vox, Pc, self.d, self.C11, self.C12, self.C21, self.C22, self.m1, self.m2, n=self.n)
        of = func_of_liquid(mox, Vf, self.rho_f, self.a, self.Df)
        func_cstr = gen_func_cstr(self.cea_path)
        cstr = self.eta* func_cstr(of, Pc)
        Re = self.rho_ox*Vox*self.d/self.mu
        # ustr_lam = func_ustr_lam(Pc, Vox, self.Tox, self.Rm, self.d, self.mu)
        # ustr_turb = func_ustr_turb(Pc, Vox, self.Tox, self.Rm, self.d, self.mu)
        dic_excond = {"mox": mox,
                      "Dt": Dt,
                      "Pc": Pc,
                      "Vox": Vox,
                      "Vf": Vf,
                      "of": of,
                      "cstr_ex": cstr,
                      "Re": Re,
                      "mox": mox
                    #   "ustr_lam": ustr_lam,
                    #   "ustr_turb": ustr_turb
                      }
        return dic_excond
    



class Gen_excond_table(Cal_excond):
    def __init__(self, cond, const, row_range=None, Dt_range=None, **kwargs):
        super().__init__(cond, const)
        self.row_range = row_range
        self.Dt_range = Dt_range

    def _input_range_(self):
        """ function to input the range of oxidizer mass flow rate and nozzle throat diameter
        
        Returns
        -------
        row_range: 1d-ndarray of float
            [kg/s] or [Pa] the range of oxidizer mass flow rate or reservoir pressure
        Dt_range: 1d-ndarray of float
            [m] the range of nozzle throat diameter
        """
        if self.cond["mode_liquid"]:
            print("\n\nInput the range of Pup [MPa], orifice up pressure, where you want to generate the table." )
            print("\ne.g. If the range is 1.0 to 5.0 MPa and the interval is 0.1 MPa\n1.0 5.0 0.1")
            tmp = list(map(lambda x: float(x) ,input().split()))
            row_range = np.arange(tmp[0], tmp[1]+tmp[2]/2, tmp[2])*1.0e+6 # generate range and convert [MPa] to [Pa]
        else:
            print("\n\nInput the range of mox [g/s], oxidizer mass flow rate, where you want to generate the table." )
            print("\ne.g. If the range is 10.0 to 20.0 g/s and the interval is 1.0 g/s\n10.0 20.0 1.0")
            tmp = list(map(lambda x: float(x) ,input().split()))
            row_range = np.arange(tmp[0], tmp[1]+tmp[2]/2, tmp[2])*1.0e-3 # generate range and convert [g/s] to [kg/s]

        print("\n\nInput the range of Dt [mm], nozzle throat diameter, where you want to generate the table." )
        print("\ne.g. If the range is 5.0 to 10.0 mm and the interval is 1.0 mm\n5.0 10.0 1.0")
        tmp = list(map(lambda x: float(x) ,input().split()))
        Dt_range = np.arange(tmp[0], tmp[1]+tmp[2]/2, tmp[2])*1.0e-3 # generate range and convert [mm] to [m]
        return row_range, Dt_range

    def gen_table(self):
        """ function to generate combustion parameter
        
        Returns
        -------
        dic_excond: dictionary of pandas.DataFrame
            dictionary of calculated combustion parameter data table
            "mox": pandas.DataFrame, [kg/s] oxidizer mass flow rate. This is assingned value
            "Dt": pandas.DataFrame, [m] nozzle throat diameter. This parameter is assigned value
            "Pc": pandas.DataFrame, [Pa] chamber pressure
            "Vox": pandas.DataFrame, [m/s] oxidzer port velocity
            "Vf": pandas.DataFrame, [m/s] axial fuel regression rate
            "of": pandas.DataFrame, [-] oxidizer to fuel mass ratio
            "cstr_ex": pandas.DataFrame, [m/s] experimental characteristic exhaust velocity
            "Re": pandas.DataFrame, [-] Reynoldz number of oxidizer at each port
            "ustr_lam": pandas.DataFrame, [m/s] friction velocity when the flow is laminar
            "ustr_turb": pandas.DataFrame, [m/s] friction velocity when the flow is turbulence
        """
        if self.row_range is None or self.Dt_range is None:
            self.row_range , self.Dt_range = self._input_range_()
        df_base = pd.DataFrame({}, index=self.row_range)
        df_pc = df_base.copy(deep=True)
        df_vox = df_base.copy(deep=True)
        df_vf = df_base.copy(deep=True)
        df_of = df_base.copy(deep=True)
        df_cstr = df_base.copy(deep=True)
        df_re = df_base.copy(deep=True)
        df_ulam = df_base.copy(deep=True)
        df_uturb = df_base.copy(deep=True)
        df_mox = df_base.copy(deep=True)
        for Dt in tqdm(self.Dt_range, desc="Dt loop", leave=True):
            Pc = np.array([])
            Vox = np.array([])
            Vf = np.array([])
            of = np.array([])
            cstr_ex = np.array([])
            Re = np.array([])
            ustr_lam = np.array([])
            ustr_turb = np.array([])
            mox = np.array([])
            for row_val in tqdm(self.row_range, desc="mox or Pup loop", leave=False):
                if self.cond["mode_liquid"]:
                    tmp = self.get_excond_liquid(row_val, Dt, Pc_init=0.1e+6)
                else:
                    tmp = self.get_excond(row_val, Dt, Pc_init=0.1e+6)
                Pc = np.append(Pc, tmp["Pc"])
                Vox = np.append(Vox, tmp["Vox"])
                Vf = np.append(Vf, tmp["Vf"])
                of = np.append(of, tmp["of"])
                cstr_ex = np.append(cstr_ex, tmp["cstr_ex"])
                Re = np.append(Re, tmp["Re"])
                if self.cond["mode_liquid"]:
                    ustr_lam = np.append(ustr_lam, np.nan)
                    ustr_turb = np.append(ustr_turb, np.nan)
                    mox = np.append(mox, tmp["mox"])
                else:
                    ustr_lam = np.append(ustr_lam, tmp["ustr_lam"])
                    ustr_turb = np.append(ustr_turb, tmp["ustr_turb"])
                    mox = np.append(mox, row_val)
            df_pc[Dt] = Pc # insert calculated Pc to data frame, df.
            df_vox[Dt] = Vox # insert calculated Vox to data frame, df.
            df_vf[Dt] = Vf # insert calculated Vf to data frame, df.
            df_of[Dt] = of # insert calculated of to data frame, df.
            df_cstr[Dt] = cstr_ex # insert calculated cstr_ex to data frame, df.
            df_re[Dt] = Re # insert calculated Re to data frame, df.
            df_ulam[Dt] = ustr_lam # insert calculated ustr_lam to data frame, df.
            df_uturb[Dt] = ustr_turb # insert calculated ustr_turb to data frame, df.
            df_mox[Dt] = mox
        self.excond_table = {"Pc": df_pc,
                        "Vox": df_vox,
                        "Vf": df_vf,
                        "of": df_of,
                        "cstr_ex": df_cstr,
                        "Re": df_re,
                        "mox": df_mox,
                        "ustr_lam": df_ulam,
                        "ustr_turb": df_uturb
                        }
        return self.excond_table

    def output(self, fldname=None):
        """ function to output the calculation result to csv file
        
        Parameters
        ----------
        fldname : string, optional
            folder name which contains calculated result as csv file, by default None
        """
        if fldname is None:
            fldname = datetime.now().strftime("%Y_%m%d_%H%M%S")   # folder name which contain animation and figure of calculation result
        os.mkdir(fldname)
        dic_json = {"CALCOND": self.cond,
                    "MODELCONST": self.model_const
                    }
        with open(os.path.join(fldname, "cond.json"), "w") as f:
            json.dump(dic_json, f, ensure_ascii=False, indent=4)
        for param, table in self.excond_table.items():
            table.to_csv(os.path.join(fldname, param+".csv"), na_rep="NaN")


def func_mox_liquid(Pup, Pc, C, beta, do, rho_ox):
    A = np.power(do, 2)*np.pi/4
    mox = C*A/np.sqrt(1-np.power(beta, 4)) *np.sqrt(2*rho_ox*(Pup-Pc))
    return mox

def func_Vox(Pc, mox, Rm, Tox, Df, a):
    Af = np.pi*np.power(Df, 2)/4
    Vox = mox*Rm*Tox/(Pc*Af*(1-a)) #[m/s]
    return(Vox)

def func_Vox_liquid(mox, rho_ox, Df, a):
    Af = np.pi*np.power(Df, 2)/4
    Vox = mox/(rho_ox*Af*(1-a)) #[m/s]
    return(Vox)

def func_Vf(Vox, Pc, C1, C2, n=1.0):
    Vf = (C1/Vox + C2)*np.power(Pc, n) #[m/s]
    return(Vf)

def func_Vf_liquid(Vox, Pc, d, C11, C12, C21, C22, m1, m2, n=1.0):
    C1 = C11*np.power(d, m1) + C12
    C2 = C21*np.power(d, m2) + C22
    Vf = (C1/Vox + C2)*np.power(Pc, n) #[m/s]
    return(Vf)

def func_of(Vox, Vf, rho_f, a, Pc, Rm, Tox):
    rho_ox = Pc/(Rm*Tox)
    of = (Vox*rho_ox*(1-a))/(rho_f*a*Vf)
    return(of)

def func_of_liquid(mox, Vf, rho_f, a, Df):
    Af = np.pi*np.power(Df, 2)/4
    of = mox/(a*Af*Vf*rho_f)
    return of

def func_Re(P, u, T, Rm, d, mu):
    rho = P/(Rm*T)
    Re = rho*u*d/mu
    return Re

def func_ustr_lam(P, u, T, Rm, d, mu):
    rho = P/(Rm*T)
    grad = 4*u/d
    tau = mu*grad
    ustr = np.sqrt(tau/rho)
    return ustr

def func_ustr_turb(P, u, T, Rm, d, mu):
    rho = P/(Rm*T)
    nu = mu/rho
    lmbd = 0.3164*np.power(u*d/nu, -1/4)
    tau = lmbd*rho*np.power(u, 2)/8
    ustr = np.sqrt(tau/rho)
    return ustr

def gen_func_cstr(cea_path):
    func = Read_datset(cea_path).gen_func("CSTAR")
    return(func)

def func_Pc_cal(of, Pc, mox, func_cstr, Dt, eta):
    cstr = func_cstr(of,Pc)
    At = np.pi*np.power(Dt, 2)/4
    Pc_cal = eta*cstr*mox*(1 + 1/of)/At
    return(Pc_cal)
   

if __name__ == "__main__":
    CALCOND = {"d": 0.5e-3,      # [m] port diameter
               "N": 132,        # [-] the number of port
               "Df": 40.7e-3,     # [m] fuel outer diameter
               "eta": 1.0,      # [-] efficiency of characteristic exhaust velocity
               "rho_f": 1191,   # [kg/m3] density of fuel
               "Rm": 259.8,     # [J/kg/K] gas constant
               "Tox": 300,      # [K] oxidizer tempereture
               "mu": 20.3e-6,  # [Pa-s] oxidizer dynamic viscosity
               "mode_liquid": True, # mode selection; using liquid oxidizer or using gasous oxidizer.
               "C": 0.61,       # discharge coefficient of orifice
               "be": 0.25,      # diamter ratio between pipe and orifice diameter
               "d_o": 1.7e-3,   # [m] diameter of orifice
               "rho_ox": 1190,  # [kg/m3] density of liquid oxidizer
            #    "cea_path": os.path.join("cea_db", "GOX_CurableResin", "csv_database")   # relative folder path to CEA .csv data-base
               "cea_path": os.path.join("cea_db", "LOX_CurableResin", "csv_database")   # relative folder path to CEA .csv data-base
               }
    CONST_VF = {"mode": "LIQUID",     # mode selectioin.
                                    # "C1C2" mode uses the following experimental formula: Vf = (C1/Vox+C2)P^n, 
                                    # "EXP" mode uses the following experimental formula: Vf = beta*P^n*Vox^m,
                                    # "PROP" mode uses the following experimental formula: Vf = alpha*P^n
                                    # "LIQUID" mode for liquid oxidizer uses the following experimental fromaula:
                                    #  Vf = (C1(d)/Vox+C2(d))P^n, C1(d) = (C11*d^m1 + C12), C2(d) = (C21*d^m2 + C22)
                "n": 1.0,           # pressure exponent
                "C1": 9.34e-8,      # SI-unit, optional. experimental constant
                # "C1": 13.9e-8,      # SI-unit, optional. experimental constant
                "C2": 2.46e-9,      # SI-unit, optional. experimental constant
                # "C2": 1.61e-9,      # SI-unit, optional. experimental constant
                "beta": None,       # SI-unit, optional. experimental constant
                "m": None,          # [-] optional. exponent of oxidizer port velocity
                "alpha": None,      # SI-unit, optional. experimental constant
                "C11": -104e-9,
                "C12": 3.55e-9,
                "m1": 0.715,
                "C21": 43.4e-9,
                "C22": -0.905e-9,
                "m2": 0.426
                }

    # inst = Cal_excond(CALCOND, CONST_VF)
    # out = inst.get_excond(mox=15.0e-3, Dt=7.0e-3, Pc_init=0.3e+6)
    # print(out)

    # MOX_RANGE = np.arange(2.0e-3, 8.1e-3, 0.1e-3)
    # DT_RANGE = np.arange(5.5e-3, 6.6e-3, 0.1e-3)
    # inst = Gen_excond_table(CALCOND, CONST_VF, mox_range=MOX_RANGE, Dt_range=DT_RANGE)

    PUP_RANGE = np.arange(1.5, 10.1, 0.1)*1e+6
    DT_RANGE = np.arange(0.5, 4.1, 0.25)*1e-3
    inst = Gen_excond_table(CALCOND, CONST_VF, row_range=PUP_RANGE, Dt_range=DT_RANGE)


    # inst = Gen_excond_table(CALCOND, CONST_VF)

    out = inst.gen_table()
    inst.output()
    print("Suceeded!")