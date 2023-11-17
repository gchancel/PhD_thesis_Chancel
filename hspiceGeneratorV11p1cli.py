
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 12:42:15 2023.

This is the version 11.1 of a HSPICE netlists generator which creates simulable ICs based on
standard cells segments models from Geoffrey CHANCEL researches about BBI.
This generator automatically connects together a set of elementary netlists
to form a complete IC, with explorable silicon substrate and power distribution network.
The power distribution networks is connected to the IC from the outside through
Metal TOP layer external pins, laid out as an external ring, which are equipotentials.
The IC generated uses a mix of Triple-Well and Dual-Well elementary blocks. They are
connected to each other in a similar way of a real STM32F439VIT6 die. However, as it is
extremely heavy to simulate a real sized IC, which equals approximately 4.5 mm x 5.5 mm,
it is proposed here to reduce its total size, according to <mul> value, which is to be
set by the user.

@author: Geoffrey Chancel alias geoff2.
"""


import os
import sys
import gzip
import time
import argparse
# import colorama
import numpy as np
# import defLibGenV11 as const
import HSPICE_gen_test_twordw_v11p1 as pgen

from progress.bar import Bar
from defLibGenV11 import (ProbeError, SUBCKT_BEGIN_TW, BLOC_TW, SUBCKT_BEGIN_DW,
                          BLOC_DW, GEN_OK, ELEM_BLOCKS_BEGIN, ELEM_BLOCKS_END, GEN_ERROR,
                          WSEG, WSEG6, HSEG, BBI_END_V11_mixed, BBI_END_V11_TW, BBI_END_V11_DW)

VER = "V11.1"
VINT = 11.1


def dual_or_triple(x_l, y_l, mul_ic, ext = False, s_x = None, s_y = None):
    """
    Determine if a dual well or triple well block is to be used.

    Parameters
    ----------
    X : list
        X coordinates of block.
    Y : list
        Y coordinates of block.
    m: float
        Ratio of circuit size.
    ext : boolean
        Use the extended function.

    Returns
    -------
    elem_bloc: str
        The corresponding elementary block string.

    """
    global tex
    global tey

    ext = True
    mul_ic = 1.0

    x_l = [float(i.replace('p', '.')) for i in x_l]
    y_l = [float(i.replace('p', '.')) for i in y_l]
    if not ext:
        if (x_l[0] >= 0 and x_l[0] < int(60 * mul_ic)):  # LEFT VERTICAL BAND
            return "elem_blocDW"
        elif (x_l[0] >= int(840 * mul_ic) and x_l[0] <= int(900 * mul_ic)):  # RIGHT VBAND
            return "elem_blocDW"
        elif (y_l[0] >= int(0 * mul_ic) and y_l[0] < int(60 * mul_ic)):  # BOTTOM V BAND
            return "elem_blocDW"
        elif (y_l[0] >= int(1040 * mul_ic) and y_l[0] <= int(1100 * mul_ic)):  # BOTTOM VBAND
            return "elem_blocDW"
        elif (y_l[0] >= 0 and y_l[0] <= int(340 * mul_ic)
              and x_l[0] >= int(330 * mul_ic) and x_l[0] < int(570 * mul_ic)):
            return "elem_blocDW"
        elif (y_l[0] > int(340 * mul_ic) and y_l[0] < int(405 * mul_ic)):
            return "elem_blocDW"
        elif (y_l[0] >= int(700 * mul_ic) and y_l[0] < int(735 * mul_ic)
              and x_l[0] >= int(60 * mul_ic) and x_l[0] < int(360 * mul_ic)):
            return "elem_blocDW"
        elif (y_l[0] >= int(700 * mul_ic) and y_l[0] < int(1100 * mul_ic)
              and x_l[0] >= int(330 * mul_ic) and x_l[0] < int(360 * mul_ic)):
            return "elem_blocDW"
        else:
            return "elem_blocTW"

    else:
        if s_x is None or s_y is None:
            raise ValueError("If ext is True, sx and sy must be float or integers")

        mulx = tex / 900
        muly = tey / 1100

        m60, m840, m900, m1040 = 60 * mulx, 840 * mulx, 900 * mulx, 1040 * muly
        m1100, m700 = 1100 * muly, 700 * muly
        m330, m540, m340, m405 = 330 * mulx, 570 * mulx, 340 * muly, 405 * muly
        m735, m700, m1100, m360 = 735 * muly, 700 * muly, 1100 * muly, 360 * mulx
        if m60 < 30:
            m60 = 30
        if tex - m840 < 60:
            m840 = m840 - 60 + (m900 - m840)
        if tey - m1040 < 35:
            m1040 = m1040 - 35 + (m1100 - m1040)
        # if m*(m840 - m570 )
        # if m

        if (x_l[0] >= 0 and x_l[0] < int(m60)):  # LEFT VERTICAL BAND
            return "elem_blocDW"
        elif (x_l[0] >= int(m840) and x_l[0] <= int(m900)):  # RIGHT VERTICAL BAND
            return "elem_blocDW"
        elif (y_l[0] >= int(0 * mul_ic) and y_l[0] < int(m60)):  # BOTTOM VERTICAL BAND
            return "elem_blocDW"
        elif (y_l[0] >= int(m1040) and y_l[0] <= int(m1100)):  # BOTTOM VERTICAL BAND
            return "elem_blocDW"
        elif (y_l[0] >= 0 and y_l[0] <= int(m340)
              and x_l[0] >= int(m330) and x_l[-1] < int(m540)):
            return "elem_blocDW"
        elif (y_l[0] > int(m340) and y_l[0] < int(m405)):
            return "elem_blocDW"
        elif (y_l[0] >= int(m700) and y_l[0] < int(m735)
              and x_l[0] >= int(m60) and x_l[0] < int(m360)):
            return "elem_blocDW"
        elif (y_l[0] >= int(m700) and y_l[0] < int(m1100)
              and x_l[0] >= int(m330) and x_l[0] < int(m360)):
            return "elem_blocDW"
        else:
            return "elem_blocTW"


class SubstrateGeneration(object):
    """
    Class to generate the substrate netlist.

    ...

    Attributes
    ----------
    subType : str
        The desired substrate type, either 'TW', 'DW' or 'mixed'

    Methods
    -------
    get_status():
        Return the status of substrate generation.
    """

    def __init__(self, sub_type):
        global bbi_gen
        global nH
        global tSub
        if sub_type == "mixed":
            bbi_gen.write(f"\n*elementary bloc {tSub}umTW")
            bbi_gen.write(SUBCKT_BEGIN_TW)
            self.__output__(sub_type)
            bbi_gen.write(BLOC_TW)
            bbi_gen.write(f"\n*elementary bloc {tSub}umDW")
            bbi_gen.write(SUBCKT_BEGIN_DW)
            self.__output__(sub_type)
            bbi_gen.write(BLOC_DW)
            bbi_gen.write("\n")
            self.__return__ = GEN_OK
        elif sub_type == "TW":
            bbi_gen.write(f"\n*elementary bloc {tSub}umTW")
            bbi_gen.write(SUBCKT_BEGIN_TW)
            self.__output__(sub_type)
            bbi_gen.write(BLOC_TW)
            bbi_gen.write("\n")
            self.__return__ = GEN_OK
        elif sub_type == "DW":
            bbi_gen.write(f"\n*elementary bloc {tSub}umDW")
            bbi_gen.write(SUBCKT_BEGIN_DW)
            self.__output__(sub_type)
            bbi_gen.write(BLOC_DW)
            bbi_gen.write("\n")
            self.__return__ = GEN_OK
        else:
            self.__return__ = GEN_ERROR

    def get_status(self):
        """
        Return the status of substrate generation.

        Parameters
        ----------
        None.

        Returns
        -------
        int
            The value of the generation status.

        """
        return self.__return__

    def __output__(self, sub_type):
        """
        ss.

        Parameters
        ----------
        sub_type : str
            The desired substrate type (TW, DW or mixed).

        Returns
        -------
        None.

        """
        global bbi_gen
        global nH
        bbi_gen.write("+")
        for prof in range(1, nH + 1, 1):
            if prof < nH + 1 - 1:
                bbi_gen.write(f"V_L_{prof} ")
            else:
                bbi_gen.write(f"V_L_{prof}\n")
        bbi_gen.write("+")
        for prof in range(1, nH + 1, 1):
            if prof < nH + 1 - 1:
                bbi_gen.write(f"V_R_{prof} ")
            else:
                bbi_gen.write(f"V_R_{prof}\n")
        for h in range(1, nH + 1, 1):
            for i in range(1, 7, 1):
                if i == 1:
                    bbi_gen.write("+")
                bbi_gen.write(f"V_F{i}_{h} ")
            for i in range(1, 7, 1):
                bbi_gen.write(f"V_RE{i}_{h} ")
            bbi_gen.write("\n")
        bbi_gen.write("+VDOWN_1 VDOWN_2 VDOWN_3 VDOWN_4 VDOWN_5 VDOWN_6\n\n")
        # Internal signals generation
        NG = []
        for h in range(1, nH + 1, 1):
            NGloc = []
            if h < nH:
                for u in range(0, 6, 1):
                    NGloc.append((u + 1) + 6 * (h - 1))
                NG.append(NGloc)
            bbi_gen.write(f"XX{h} ")
            if h == nH:
                bbi_gen.write("VDOWN_1 VDOWN_2 VDOWN_3 VDOWN_4 VDOWN_5 VDOWN_6 ")
            else:
                for ng in NGloc:
                    bbi_gen.write(f"NG{ng} ")
            for i in range(1, 7, 1):
                bbi_gen.write(f"V_F{i}_{h} ")
            bbi_gen.write(f"V_L_{h} V_R_{h} ")
            for i in range(1, 7, 1):
                bbi_gen.write(f"V_RE{i}_{h} ")
            if h == 1:
                for i in range(6):
                    bbi_gen.write("vepi ")
            else:
                for ng in NG[h - 1 - 1]:
                    bbi_gen.write(f"NG{ng} ")
            bbi_gen.write(f"VSUBC{h} elementary_blocx6\n")


def main():
    """
    Execute the main code.

    Raises
    ------
    ValueError
        DESCRIPTION.
    ProbeError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # =============================================================================
    # Global variables
    # =============================================================================
    global bbi_gen
    global nH
    global tSub
    global nC
    global nL
    global tex
    global tey

    # =============================================================================
    # CLI INTERFACE
    # =============================================================================
    desc = '''Generate a netlists of standard-cells according to user settings.
    One can specify the required substrate type, the model type (bad gnd, good gnd, good gnd
                                                                  + impedance matching).
    The user has to specify the IC size (see the parameters). When choosing a size,
    the program might resize the IC according to elemetary models size, as they are
    unbreakable. Please also note that 1:1 ratios might often produce impossible results.
    If the generator does not manage to generate what you request, it will notify you and abort
    the generation process, asking you to change the IC size.
    In addition to this,, <mul> can be specified to change the global size of the IC.
    Note that some configurations may not work properly.
    '''

    parser = argparse.ArgumentParser(prog = "python ./hspiceGeneratorV11cli.py",
                                     description = desc,
                                     epilog = "gchancel2023")

    parser.add_argument("-s", dest = "subType",
                        help = "Substrate type, can take 3 values: <mixed>, <tw> or <dw>",
                        required = True)

    parser.add_argument("-m", dest = "model",
                        help = "Model type, can take 3 values: <0>, <1>, or <2>",
                        default = "2",
                        type = int,
                        required = False)

    parser.add_argument("-t", dest = "tSub",
                        help = '''Substrate thickness, integer between 10 and infinity,
                                multiple of 10''',
                        required = True, type = int)

    parser.add_argument("-eT", dest = "elemThickness",
                        help = '''Thickness of the elementary substrate block to choose.
                                A lower value brings more precision, a higher value decreases
                                precision but increases simulation speed for a constat tSub.
                                Can take any floating point value
                                between 10 and 1000 (in µm).''',
                        required = True, type = int)

    parser.add_argument("--mul",
                        help = '''Resize factor (float)''',
                        default = "1.0",
                        type = float,
                        required = False)

    parser.add_argument("-v", dest = "vpUU",
                        help = '''Amplitude of the applied pulse (in V)''',
                        type = float,
                        required = True)

    parser.add_argument("-p", dest = "pW",
                        help = '''Pulse width (in ns)''',
                        type = float,
                        required = True)

    parser.add_argument("-tr", dest = "tFR",
                        help = '''Rise/fall time (ns)''',
                        type = float,
                        required = True)

    parser.add_argument("--step", dest = "sim_step",
                        help = '''Simulation time step (specifiy the unit without space):
                            example: --step 50ps''',
                        required = True)

    parser.add_argument("--time", dest = "sim_time",
                        help = '''Simulation duration (specify the unit without space):
                            example: --time 40ns''',
                        required = True)

    parser.add_argument("--auto", dest = "auto",
                        help = '''Execute the generator without user prompt.''',
                        required = False,
                        default = "False",
                        type = str)

    parser.add_argument("--res", dest = "res",
                        help = '''Display the result graphically for visual inspection.''',
                        required = False,
                        default = "True",
                        type = str)

    parser.add_argument("-tex", dest = "tex",
                        help = '''X size''',
                        type = float,
                        required = True)

    parser.add_argument("-tey", dest = "tey",
                        help = '''Y size''',
                        type = float,
                        required = True)

    parser.add_argument("--psInt", dest = "plotSubInt",
                        help = '''Output all substrate interior signals into output file''',
                        type = str,
                        default = "False",
                        required = False)

    parser.add_argument("--compress", dest = "compress",
                        help = '''Compress output file and remove non compressed file.''',
                        type = str,
                        default = "False",
                        required = False)

    parser.add_argument("--af", dest = "appendFname",
                        help = '''Append custom string to the end of filename.''',
                        type = str,
                        default = "",
                        required = False)

    args = parser.parse_args()
    print("\nThe settings you entered are the following:")
    print(f"SubType = {args.subType}")
    print(f"Model = {args.model}")
    print(f"Thickness = {args.tSub} µm")
    print(f"Elementary thickness = {args.elemThickness} µm")
    print(f"Mult = {args.mul}")
    print(f"Vpulse = {args.vpUU} V")
    print(f"Pulse width = {args.pW} ns")
    print(f"Rise/fall time = {args.tFR} ns")
    print(f"Simulation time step = {args.sim_step}")
    print(f"Simulation duration = {args.sim_time}")
    print(f"Auto mode = {args.auto}")
    print(f"Display mode = {args.res}")
    print(f"X size = {args.tex}")
    print(f"Y size = {args.tey}")
    print(f"Plot Substrate Interior = {args.plotSubInt}")

    if not args.auto == "True":
        input("\n!!! PRESS ENTER TO CONTINUE !!!\n")
        locked = True
        msg = '''Are you sure to continue? Please answer with 'y'
        to continue or with 'n' to abort.'''
        while locked:
            print(msg)
            choice = input().lower()
            if choice == 'n':
                sys.exit(0)
            elif choice == 'y':
                locked = False

    print("Netlist generation in progess")

    # =============================================================================
    # IC SIZE SECTION (USER CUSTOMIZATION)
    # =============================================================================
    # mul = 0.35  # Set the scale used to reduce or increase original IC size.
    mul = args.mul  # Set the scale used to reduce or increase original IC size.
    # tex = mul * 900  # Final x size of the IC. DO NOT CHANGE ON ANY SITUATION !!!
    # tey = mul * 1100  # Final y size of the IC. DO NOT CHANGE ON ANY SITUATION !!!
    tex = args.tex * mul  # Final x size of the IC. DO NOT CHANGE ON ANY SITUATION !!!
    tey = args.tey * mul  # Final y size of the IC. DO NOT CHANGE ON ANY SITUATION !!!

    nC = int(tex / 30.0)  # According to the previously computed size, final number of columns.
    nL = int(tey / 5.0)  # According to the previously computed size, final number of lines.
    tex = int(nC * 30.0)
    tey = int(nL * 5.0)
    tSub = args.tSub  # Substrate thickness in µm.
    eT = args.elemThickness  # Thickness of an elementary substrate block.
    nH = int(tSub / eT)  # (Number of silicon substrate layers).
    print(f"\033[1;32mNumber of substrate layers = {nH}")
    resultingTSUB = nH * eT
    print(f"Resulting Substrate chosen = {resultingTSUB} µm\033[1;0m")

    if args.subType == "mixed":
        mixedDWTW = True
        TWorDW = True
    elif args.subType == "tw":
        mixedDWTW = False
        TWorDW = True
    elif args.subType == "dw":
        mixedDWTW = False
        TWorDW = False

    bar = Bar("Generation", max = nL * nC)

    model = args.model  # Can take 3 different values : 0, 1 or 2.
    '''
    Model 0 corresponds to the first BBI attempts, with bad grounding and impedance mismatch.
    Model 1 corresponds to Good grounding and impedance mismatch.
    Model 2 corresponds to Good Grounding + Impedance matching.
    '''

    vpUU = args.vpUU  # Amplitude of the voltage pulse.
    pW = args.pW  # Pulse width (in ns)
    tFR = args.tFR  # Fall/Rise time (in ns)
    sim_step = args.sim_step  # Time step of the simulation.
    sim_time = args.sim_time  # Duration of the simulation.

    Rm1 = '26'
    Rmup = '5'
    Cdecp = '2.25f'
    Ra = '2'
    Cgp = '35.2f'
    Cgn = '25.2f'
    Rp = '9.57'
    Rn = '5.3'
    Cnw = '20f'
    RcontactTW = '3.1k'
    RcontactDW = '3.1k'

    probe_count = 0

    # =============================================================================
    # BBI PROBE SECTION: !!! DO NOT TOUCH, NOT FINALIZED AND HIGHLY PRONE TO PRODUCE
    # UNEXPECTED RESULTS. !!! (BETA TEST, INCONSISTENT)
    # =============================================================================
    probe_center = True
    c_probe = True  # BBI PROBE CONNECTION (ONLY PARAMETER ALLOWED TO BE MODIFIED WITH CARE).
    if probe_center:
        # prbX = tex / 2.0
        # prbY = tey / 2.0
        prbX = int((nC * 30.0) / 2.0)
        prbY = int((nL * 5.0) / 2.0)
        # print(f"PROBE X = {prbX}")
        # print(f"PROBE Y = {prbY}")
    else:  # Custom probe position (may leads to inconsistencies, use with care).
        prbX = 210 * mul
        prbY = 890 * mul
    # DEFAULT PROBE IS 30 µm * 30 µm, not anything else. Other sizes will be added in
    # future versions of this script. Default location is the center of the IC, no choice.

    # =============================================================================
    # Development test variables, do not touch. (ALPHA TEST, unstable)
    # =============================================================================
    # cgnd = 0
    # cvdd = 0

    # xvdd = []
    # xgnd = []

    x_trace = []
    y_trace = []

    if mul == 1:
        ext_loc = False
    else:
        ext_loc = True

    ring = False
    if ring:
        print('RING POWER SUPPLY')
    else:
        print('TOP/BOTTOM POWER SUPPLIES (VDD AND GND TOGETHER)')
    print('X size =', tex, 'µm')
    print('Y size =', tey, 'µm')
    print('Nb of columns =', nC)
    print('Nb of lines =', nL)
    print('Depth of the IC (FIXED VALUE) =', int(nH * 10), 'µm')
    print('Probe connected =', c_probe)
    if c_probe:
        print('!CENTRAL PROBE, square 30 µm * 30 µm, not editable! (for now at least)')
    print('Mixed TW and DW cells =', mixedDWTW)
    if not mixedDWTW:
        if TWorDW:
            print('Exclusive TRIPLE-WELL')
        else:
            print('Exclusive DUAL-WELL')

    # =============================================================================
    # BEGINNING OF THE SCRIPT
    # =============================================================================
    filename = "./GENERATED_NETLISTS/"
    if not c_probe:
        if mixedDWTW:
            filename += f'./bbiGenVer{VINT}_T{int(tSub)}eT{eT}x'
            filename += f'{int(tex)}y{int(tey)}mixDwTwOp{args.appendFname}.sp'
        else:
            if TWorDW:
                filename += f'./bbiGenVer{VINT}_T{int(tSub)}eT{eT}x'
                filename += f'{int(tex)}y{int(tey)}excTwOp{args.appendFname}.sp'
            else:
                filename += f'./bbiGenVer{VINT}_T{int(tSub)}eT{eT}x'
                filename += f'{int(tex)}y{int(tey)}excDwOp{args.appendFname}.sp'
    else:
        filename += f'./bbiGenVer{VINT}_T{int(tSub)}eT{eT}'
        filename += f'x{int(tex)}y{int(tey)}Px{prbX}Py{prbY}'
        if mixedDWTW:
            filename += f"mixDwTwTranVp{int(vpUU)}Pw{pW}tFR{tFR}M{model}{args.appendFname}.sp"
        else:
            if TWorDW:
                filename += f"ExcTwVp{int(vpUU)}Pw{pW}tFR{tFR}M{model}{args.appendFname}.sp"
            else:
                filename += f"ExcDwVp{int(vpUU)}Pw{pW}tFR{tFR}M{model}{args.appendFname}.sp"

    bbi_gen = open(filename, 'w', encoding = 'utf-8')

    bbi_gen.write(f'BBI GENERATOR {VER}\n')
    bbi_gen.write('\n.param vpUU=' + str(vpUU) + '\n')
    bbi_gen.write('\n.param RM1=' + Rm1 + ' Rmup=' + Rmup + ' Cdecp=' + Cdecp + ' Ra=' + Ra)
    bbi_gen.write(' Cgp=' + Cgp + ' Cgn=' + Cgn + ' RP=' + Rp + ' RN=' + Rn + ' CNW=' + Cnw)
    bbi_gen.write(' RcontactTW=' + RcontactTW + ' RcontactDW=' + RcontactDW)
    bbi_gen.write('\n\n')

    bbi_gen.write('vBBI_source VDOWN_BBI_source GND pwl(')
    bbi_gen.write('0 0 3n 0 ' + str(3 + tFR) + 'n vpUU ' + str(3 + pW) + 'n vpUU ')
    bbi_gen.write(str(3 + pW + tFR) + 'n 0)\n\n')

    bbi_gen.write("Rvbbi VDOWN_BBI_source VDOWN_BBI_trans '100e9*(TIME < 3n && TIME > ")
    bbi_gen.write(str(3 + pW + tFR) + "n) + 0.00001'\n\n")

    bbi_gen.write("Tx1 VDOWN_BBI_trans GND VDOWN_BBI GND z0=50 td=3n\n\n")

    if model == 2:
        bbi_gen.write("* Impedance Matching Enabled\n")
        bbi_gen.write("R50 VDOWN_BBI GND 50\n\n")

    bbi_gen.write("vBBI VDOWN_BBI VDOWN_BBI_STM32 dc=0\n\n")

    bbi_gen.write("V_ALIM VDD GND dc=1.2\n")
    bbi_gen.write("VmeasGND GNDmeasSTM 0 dc=0\n\n")
    if model == 0:
        bbi_gen.write("R150_0 GNDm GNDmeasSTM 150\n")
        bbi_gen.write("R150_1 GNDg 0 150\n\n")
    elif model == 1 or model == 2:
        bbi_gen.write("R150_0 GNDm GNDmeasSTM 0.01\n")
        bbi_gen.write("R150_1 GNDg 0 0.01\n\n")

    bbi_gen.write(f'* xSize {tex} µm\n')
    bbi_gen.write(f'* ySize {tey} µm\n')
    bbi_gen.write(f'* tSize {tSub} µm\n\n')

    if mixedDWTW:
        bbi_gen.write('* MIXED TRIPLE-WELL AND DUAL-WELL\n')
    else:
        if TWorDW:
            bbi_gen.write('* EXCLUSIVE TRIPLE-WELL\n')
        else:
            bbi_gen.write('* EXCLUSIVE DUAL-WELL\n')

    bbi_gen.write(f"* MODEL {model}\n")

    # subElemResDepth10u = 2000
    # subElemResOther10u =

    rho = 0.01  # Ohm*meter
    ResVert = (rho * ((eT * 1e-6) / 2)) / (5e-6 * 5e-6)
    ResOther = (rho * 2.5e-6) / (5e-6 * (eT * 1e-6))

    bbi_gen.write(ELEM_BLOCKS_BEGIN)
    bbi_gen.write(f"""R1 U N001 {ResVert}
R2 N001 D {ResVert}
R3 Re N001 {ResOther}
R4 N001 F {ResOther}
R5 N001 L {ResOther}
R6 R N001 {ResOther}""")
    bbi_gen.write(ELEM_BLOCKS_END)

#     bbi_gen.write(ELEM_BLOCKS_BEGIN)
#     bbi_gen.write(f"""R1 U N001 {subElemResDepth10u * (eT / 10)}
# R2 N001 D {subElemResDepth10u * (eT / 10)}
# R3 Re N001 {subElemResOther10u * (eT / 10)}
# R4 N001 F {subElemResOther10u * (eT / 10)}
# R5 N001 L {subElemResOther10u * (eT / 10)}
# R6 R N001 {subElemResOther10u * (eT / 10)}""")
#     bbi_gen.write(ELEM_BLOCKS_END)

    # rho = 0.01  # Ohm*meter
    # subRT = 2750
    # subRB = 2750
    # subRL = 2750
    # subRR = 2750
    # subRRe = 2750
    # subRF = 2750

#     bbi_gen.write(ELEM_BLOCKS_BEGIN)
#     bbi_gen.write(f"""R1 U N001 {subRT}
# R2 N001 D {subRT}
# R3 Re N001 {subRT}
# R4 N001 F {subRT}
# R5 N001 L {subRT}
# R6 R N001 {subRT}""")
#     bbi_gen.write(ELEM_BLOCKS_END)

    if mixedDWTW:
        gen_status = SubstrateGeneration("mixed").get_status()
    elif not mixedDWTW and TWorDW:
        gen_status = SubstrateGeneration("TW").get_status()
    elif not mixedDWTW and not TWorDW:
        gen_status = SubstrateGeneration("DW").get_status()
    else:
        gen_status = SubstrateGeneration("error").get_status()

    if gen_status == GEN_ERROR:
        raise ValueError

    print("Main generation beginning")
    for le in range(nL):
        for co in range(nC):
            # ================================================================================
            # Generating first the local X and Y coordinates lists needed in nets names.
            # ================================================================================
            X = [co * WSEG + 0 * (WSEG6 / 2),
                 co * WSEG + 1 * (WSEG6 / 2),
                 co * WSEG + 3 * (WSEG6 / 2),
                 co * WSEG + 5 * (WSEG6 / 2),
                 co * WSEG + 7 * (WSEG6 / 2),
                 co * WSEG + 9 * (WSEG6 / 2),
                 co * WSEG + 11 * (WSEG6 / 2),
                 co * WSEG + 12 * (WSEG6 / 2)]
            X = [float(i) for i in X]  # To be sure they are all floats.
            Xf = X.copy()
            X = [str(i).replace('.', 'p') for i in X]

            Y = [le * HSEG,
                 le * HSEG + 1 / 2 * HSEG,
                 (le + 1) * HSEG]
            Y = [float(i) for i in Y]  # To be sure they are all floats.
            Yf = Y.copy()
            Y = [str(i).replace('.', 'p') for i in Y]

            x_trace.append(Xf)
            y_trace.append(Yf)

            # ================================================================================
            # VDD nets (Metal TOP and Metal 1)
            # ================================================================================
            vic = 'V_' + X[1] + '_' + Y[2] + '_0'
            vic1 = 'V_' + X[6] + '_' + Y[2] + '_0'
            vi1c = 'V_' + X[1] + '_' + Y[0] + '_0'
            vi1c1 = 'V_' + X[6] + '_' + Y[0] + '_0'
            vlt = 'VM_' + X[0] + '_' + Y[1] + '_0'
            vrt = 'VM_' + X[7] + '_' + Y[1] + '_0'

            # ================================================================================
            # GND nets (Metal TOP and Metal 1)
            # ================================================================================
            gic = 'G_' + X[1] + '_' + Y[2] + '_0'
            gic1 = 'G_' + X[6] + '_' + Y[2] + '_0'
            gi1c = 'G_' + X[1] + '_' + Y[0] + '_0'
            gi1c1 = 'G_' + X[6] + '_' + Y[0] + '_0'
            glb = 'GM_' + X[0] + '_' + Y[1] + '_0'
            grb = 'GM_' + X[7] + '_' + Y[1] + '_0'

            # ================================================================================
            # LEFT NETS
            # ================================================================================
            VL = []
            for left in range(1, nH + 1, 1):
                VL.append('VLR_' + X[0] + '_' + Y[1] + '_' + str(left))

            # ================================================================================
            # RIGHT NETS
            # ================================================================================
            VR = []
            for right in range(1, nH + 1, 1):
                VR.append('VLR_' + X[7] + '_' + Y[1] + '_' + str(right))

            # ================================================================================
            # REAR NETS
            # Will be a 3D numpy array,
            # ================================================================================
            VRE = np.empty((nH, 6), dtype = object)
            for rear in range(1, nH + 1, 1):
                for floc in range(6):
                    VRE[rear - 1, floc] = 'VFRE' + str(floc + 1) + '_'
                    VRE[rear - 1, floc] += X[floc + 1] + '_'
                    VRE[rear - 1, floc] += Y[2] + '_' + str(rear)

            # ================================================================================
            # FRONT NETS
            # ================================================================================
            VFR = np.empty((nH, 6), dtype = object)
            for front in range(1, nH + 1, 1):
                for floc in range(6):
                    VFR[front - 1, floc] = 'VFRE' + str(floc + 1) + '_'
                    VFR[front - 1, floc] += X[floc + 1] + '_'
                    VFR[front - 1, floc] += Y[0] + '_' + str(front)

            # ================================================================================
            # EXTREME BOTTOM NETS
            # ================================================================================
            VDOWN = []
            for bot in range(6):
                VDOWN.append('VDOWN_' + str(bot + 1) + '_' + X[bot + 1] + '_' + Y[1])

            # ================================================================================
            # EVALUATING AND CONNECTINF EXTERNAL POWER SUPPLY, GROUND, AND BBI PROBE.
            # ================================================================================
            if ring:
                alimplace = False
                if np.isclose(Xf[0], 0.0):  # LEFT SIDE
                    vic = 'VDD'
                    vi1c = 'VDD'
                    gic = 'GNDm'
                    gi1c = 'GNDm'
                    alimplace = True
                elif np.isclose(Xf[-1], tex):  # RIGHT SIDE
                    vic1 = 'VDD'
                    vi1c1 = 'VDD'
                    gic1 = 'GNDm'
                    gi1c1 = 'GNDm'
                    alimplace = True
                elif np.isclose(Yf[0], 0.0):  # FRONT
                    vi1c = 'VDD'
                    vi1c1 = 'VDD'
                    gi1c = 'GNDm'
                    gi1c1 = 'GNDm'
                    alimplace = True
                elif np.isclose(Yf[-1], tey):  # REAR
                    vic = 'VDD'
                    vic1 = 'VDD'
                    gic = 'GNDm'
                    gic1 = 'GNDm'
                    alimplace = True
                else:
                    pass  # CAS GÉNÉRAL
            else:
                alimplace = False
                if np.isclose(Yf[0], 0.0):  # FRONT
                    vi1c = 'VDD'
                    vi1c1 = 'VDD'
                    gi1c = 'GNDm'
                    gi1c1 = 'GNDm'
                    alimplace = True
                    # print('alimplace')
                elif np.isclose(Yf[-1], tey):  # REAR
                    vic = 'VDD'
                    vic1 = 'VDD'
                    gic = 'GNDm'
                    gic1 = 'GNDm'
                    alimplace = True
                    # print('alimplace')
                else:
                    pass  # CAS GÉNÉRAL

            probe_here = False
            # print(Xf[0], tex / 2.0)
            # print(Yf[1], (tey / 2.0) + 15.0)
            # print(Yf[1], (tey / 2.0) - 15.0)
            # print(c_probe
            #       and np.isclose(Xf[0], tex / 2.0)
            #       and Yf[1] <= (tey / 2.0) + 15.0
            #       and Yf[1] >= (tey / 2.0) - 15.0)
            if probe_center:
                if (c_probe
                   and (np.isclose(Xf[0], prbX) or np.isclose(Xf[2], prbX))
                   and (Yf[0] <= (prbY) + 15.0)
                   and (Yf[0] > (prbY) - 15.0)):
                    VDOWN = ['VDOWN_BBI_STM32' for i in VDOWN]
                    # print(Xf[0], Yf[0], 'PROBE')
                    probe_here = True
                    probe_count += 1
                else:
                    pass
            else:
                if (c_probe
                   and (np.isclose(Xf[0], prbX) or np.isclose(Xf[2], prbX))
                   and Yf[0] <= prbY + 15.0
                   and Yf[0] > prbY - 15.0):
                    VDOWN = ['VDOWN_BBI_STM32' for i in VDOWN]
                    # print(Xf[0], Yf[0], 'PROBE')
                    probe_here = True
                    probe_count += 1
                else:
                    pass

            # ================================================================================
            # WRITING INTO FINAL FILE
            # ================================================================================
            if not probe_here:
                if not alimplace:
                    if mixedDWTW:
                        dualORtriple = dual_or_triple(X, Y, mul,
                                                      ext = ext_loc, s_x = tex, s_y = tey)
                        if dualORtriple == 'elem_blocTW':
                            bbi_gen.write('XblocTW_x' + X[0] + '_y' + Y[0] + ' ')
                        elif dualORtriple == 'elem_blocDW':
                            bbi_gen.write('XblocDW_x' + X[0] + '_y' + Y[0] + ' ')
                    else:
                        if TWorDW:
                            bbi_gen.write('XblocTW_x' + X[0] + '_y' + Y[0] + ' ')
                        else:
                            bbi_gen.write('XblocDW_x' + X[0] + '_y' + Y[0] + ' ')
                    # bbi_gen.write('Xbloc_x' + X[0] + '_y' + Y[0] + ' ')
                else:
                    if mixedDWTW:
                        dualORtriple = dual_or_triple(X, Y, mul,
                                                      ext = ext_loc, s_x = tex, s_y = tey)
                        if dualORtriple == 'elem_blocTW':
                            bbi_gen.write('XblocTWALM_x' + X[0] + '_y' + Y[0] + ' ')
                        elif dualORtriple == 'elem_blocDW':
                            bbi_gen.write('XblocDWALM_x' + X[0] + '_y' + Y[0] + ' ')
                    else:
                        if TWorDW:
                            bbi_gen.write('XblocTWALM_x' + X[0] + '_y' + Y[0] + ' ')
                        else:
                            bbi_gen.write('XblocDWALM_x' + X[0] + '_y' + Y[0] + ' ')
                    # bbi_gen.write('Xbloc_x' + X[0] + '_y' + Y[0] + ' ')
            else:
                # print('xblocpro')
                if not alimplace:
                    if mixedDWTW:
                        dualORtriple = dual_or_triple(X, Y, mul,
                                                      ext = ext_loc, s_x = tex, s_y = tey)
                        if dualORtriple == 'elem_blocTW':
                            bbi_gen.write('XblocTWPRB_x' + X[0] + '_y' + Y[0] + ' ')
                        elif dualORtriple == 'elem_blocDW':
                            bbi_gen.write('XblocDWPRB_x' + X[0] + '_y' + Y[0] + ' ')
                    else:
                        if TWorDW:
                            bbi_gen.write('XblocTWPRB_x' + X[0] + '_y' + Y[0] + ' ')
                        else:
                            bbi_gen.write('XblocDWPRB_x' + X[0] + '_y' + Y[0] + ' ')
                else:
                    if mixedDWTW:
                        dualORtriple = dual_or_triple(X, Y, mul,
                                                      ext = ext_loc, s_x = tex, s_y = tey)
                        if dualORtriple == 'elem_blocTW':
                            bbi_gen.write('XblocTWALMPRB_x' + X[0] + '_y' + Y[0] + ' ')
                        elif dualORtriple == 'elem_blocDW':
                            bbi_gen.write('XblocDWALMPRB_x' + X[0] + '_y' + Y[0] + ' ')
                    else:
                        if TWorDW:
                            bbi_gen.write('XblocTWALMPRB_x' + X[0] + '_y' + Y[0] + ' ')
                        else:
                            bbi_gen.write('XblocDWALMPRB_x' + X[0] + '_y' + Y[0] + ' ')

            bbi_gen.write(vic + ' ')
            bbi_gen.write(vic1 + ' ')
            bbi_gen.write(vi1c + ' ')
            bbi_gen.write(vi1c1 + ' ')

            bbi_gen.write('\n+')

            bbi_gen.write(gic + ' ')
            bbi_gen.write(gic1 + ' ')
            bbi_gen.write(gi1c + ' ')
            bbi_gen.write(gi1c1 + ' ')

            bbi_gen.write('\n+')

            bbi_gen.write(vlt + ' ')
            bbi_gen.write(vrt + ' ')
            bbi_gen.write(glb + ' ')
            bbi_gen.write(grb + ' ')

            bbi_gen.write('\n+')

            for line in VL:
                bbi_gen.write(line + ' ')

            bbi_gen.write('\n+')

            for line in VR:
                bbi_gen.write(line + ' ')

            bbi_gen.write('\n+')

            for depth in range(nH):
                for line in VFR[depth]:
                    bbi_gen.write(line + ' ')
                bbi_gen.write('\n+')
                for line in VRE[depth]:
                    bbi_gen.write(line + ' ')

            bbi_gen.write('\n+')

            for line in VDOWN:
                bbi_gen.write(line + ' ')

            if mixedDWTW:
                bbi_gen.write(dual_or_triple(X, Y, mul, ext = ext_loc, s_x = tex, s_y = tey))
            else:
                if TWorDW:
                    bbi_gen.write("elem_blocTW")
                else:
                    bbi_gen.write("elem_blocDW")

            bbi_gen.write('\n\n')
            bar.next()

    # print(f"\nPROBE COUNT = {probe_count}")

    bbi_gen.write('\n.tran ' + sim_step + ' ' + sim_time)
    if mixedDWTW:
        bbi_gen.write(BBI_END_V11_mixed)
    elif TWorDW:
        bbi_gen.write(BBI_END_V11_TW)
    else:
        bbi_gen.write(BBI_END_V11_DW)

    if args.plotSubInt == "True":
        bbi_gen.write('.print v(*.VSUBC*)\n')
        bbi_gen.write('.print v(VDOWN*)\n')
    bbi_gen.write("\n.end")
    bbi_gen.close()

    try:
        if probe_count != 6:
            raise ProbeError
        if args.res == "True":
            print("\033[1;32m\nGeneration ended without errors.")
            print("Drawing resulting netlist top-view, please wait...\033[1;0m")
            bbi_gtest = open(filename, 'r')
            bbi_results = pgen.plotGEN(bbi_gtest.read())
            bbi_gtest.close()
            if bbi_results:
                pass
            else:
                print("Plot error\033[1;0m")
        else:
            print("\033[1;32m\nGeneration ended without errors, exiting generator.\033[1;0m")
    except ProbeError:
        os.remove(filename)
        print("\033[1;31mProbe incorrectly placed.")
        print("Please choose another IC size.\033[1;0m")

    if args.compress == "True":
        print("\033[1;33mCompressing text file, please wait...\033[1;0m")
        with open(filename, 'r') as bbi_gen, gzip.open(f'{filename}.gz', 'wt') as bbi_gen_gz:
            bbi_gen_gz.write(bbi_gen.read())
        print("\033[1;34mFile compressed successfully\033[1;0m")
        print("\033[1;35mRemoving uncompressed .sp file!\033[1;0m")
        os.remove(filename)


if __name__ == '__main__':
    startTime = time.time()
    main()
    endTime = time.time()
    print(f"Elasped time for execution: {round(endTime - startTime, 2)} s")
