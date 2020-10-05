# Condition file for main_design.py
import os

PARAM_EX = {"d": 0.3e-3,     # [m] port diameter
            "N": 432,        # [-] the number of port
            "Df": 38.0e-3,   # [m] fuel outer diameter
            "eta": 0.83,     # [-] efficiency of characteristic exhaust velocity
            "rho_f": 1191,   # [kg/m3] density of fuel
            "Rm": 259.8,     # [J/kg/K] gas constant
            "Tox": 300,      # [K] oxidizer tempereture
            "mu": 20.3e-6    # [Pa-s] oxidizer dynamic viscosity
            }

PARAM_LIQUID = {"mode_liquid": False, # mode selection; using liquid oxidizer or using gasous oxidizer.
                "C": 0.61,       # discharge coefficient of orifice
                "be": 0.25,      # diamter ratio between pipe and orifice diameter
                "d_o": 1.7e-3,   # [m] diameter of orifice
                "rho_ox": 1190   # [kg/m3] density of liquid oxidizer
                }

PARAM_CEA = {"cea_path": os.path.join("cea_db", "sample", "csv_database"),   # relative folder path to CEA .csv data-base
             "wt_other": 0.0        # [%] mass fraction of other chemical species, which are not oxidizer or fuel.
                                    # you should assign the cea data base which consider the mixture of other chemical species if "wt_other" is not zero.
            }

CONST_VF = {"mode": "C1C2",     # mode selectioin.
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
            "C11": -104e-9,     # SI-unit, optional. experimental constant
            "C12": 3.55e-9,     # SI-unit, optional. experimental constant
            "m1": 0.715,        # [-] optional. exponent of port diameter
            "C21": 43.4e-9,     # SI-unit, optional. experimental constant
            "C22": -0.905e-9,   # SI-unit, optional. experimental constant
            "m2": 0.426         # [-] optional. exponent of port diameter
            }