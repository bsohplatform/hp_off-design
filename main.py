from HP_dataclass import *
from COMP_module import *
from HX_module import *
from CoolProp.CoolProp import PropsSI

class VCHP_off:
    def __init__(self, evap_BC, cond_BC):
        self.evap_BC = evap_BC
        self.cond_BC = cond_BC
        
    def Input_Processing(self, InCond, OutCond, InEvap, OutEvap):
        if InCond.p <= 0.0:
            InCond.p = 101300.0
            
        if OutCond.p <= 0.0:
            OutCond.p = 101300.0
        
        if InCond.T <= 0.0:
            noCond = 0
        
        if OutCond.T <= 0.0:
            noCond = 1
        
        if InCond.m <= 0.0 and OutCond.m <= 0:
            noCond = 2
        else:
            if InCond.m == 0:
                InCond.m = OutCond.m
            else:
                OutCond.m = InCond.m
        
        if InEvap.p <= 0.0:
            InEvap.p = 101300.0
            
        if OutEvap.p <= 0.0:
            OutEvap.p = 101300.0
        
        if InEvap.T <= 0.0:
            noEvap = 0
        
        if OutEvap.T <= 0.0:
            noEvap = 1
        
        if InEvap.m <= 0.0 and OutEvap.m <= 0:
            noEvap = 2
        else:
            if InEvap.m == 0:
                InEvap.m = OutEvap.m
            else:
                OutEvap.m = InEvap.m
        
        if noCond == 0:
            OutCond.c = Aux_fn.PropCal(OutCond, 'C', 'T', 'P')
            OutCond.v = Aux_fn.PropCal(OutCond, 'V', 'T', 'P')
            OutCond.pr = Aux_fn.PropCal(OutCond, 'Prandtl', 'T', 'P')
            OutCond.l = Aux_fn.PropCal(OutCond, 'L', 'T', 'P')
            OutCond.d = Aux_fn.PropCal(OutCond, 'D', 'T', 'P')
        elif noCond == 1:
            InCond.c = Aux_fn.PropCal(InCond, 'C', 'T', 'P')
            InCond.v = Aux_fn.PropCal(InCond, 'V', 'T', 'P')
            InCond.pr = Aux_fn.PropCal(InCond, 'Prandtl', 'T', 'P')
            InCond.l = Aux_fn.PropCal(InCond, 'L', 'T', 'P')
            InCond.d = Aux_fn.PropCal(InCond, 'D', 'T', 'P')            
        else:
            InCond.c = Aux_fn.PropCal(InCond, 'C', 'T', 'P')
            InCond.v = Aux_fn.PropCal(InCond, 'V', 'T', 'P')
            InCond.pr = Aux_fn.PropCal(InCond, 'Prandtl', 'T', 'P')
            InCond.l = Aux_fn.PropCal(InCond, 'L', 'T', 'P')
            InCond.d = Aux_fn.PropCal(InCond, 'D', 'T', 'P')
            OutCond.c = Aux_fn.PropCal(OutCond, 'C', 'T', 'P')
            OutCond.v = Aux_fn.PropCal(OutCond, 'V', 'T', 'P')
            OutCond.pr = Aux_fn.PropCal(OutCond, 'Prandtl', 'T', 'P')
            OutCond.l = Aux_fn.PropCal(OutCond, 'L', 'T', 'P')
            OutCond.d = Aux_fn.PropCal(OutCond, 'D', 'T', 'P')
        
        if noEvap == 0:
            OutEvap.c = Aux_fn.PropCal(OutEvap, 'C', 'T', 'P')
            OutEvap.v = Aux_fn.PropCal(OutEvap, 'V', 'T', 'P')
            OutEvap.pr = Aux_fn.PropCal(OutEvap, 'Prandtl', 'T', 'P')
            OutEvap.l = Aux_fn.PropCal(OutEvap, 'L', 'T', 'P')
            OutEvap.d = Aux_fn.PropCal(OutEvap, 'D', 'T', 'P')
        elif noEvap == 1:
            InEvap.c = Aux_fn.PropCal(InEvap, 'C', 'T', 'P')
            InEvap.v = Aux_fn.PropCal(InEvap, 'V', 'T', 'P')
            InEvap.pr = Aux_fn.PropCal(InEvap, 'Prandtl', 'T', 'P')
            InEvap.l = Aux_fn.PropCal(InEvap, 'L', 'T', 'P')
            InEvap.d = Aux_fn.PropCal(InEvap, 'D', 'T', 'P')
        else:
            InEvap.c = Aux_fn.PropCal(InEvap, 'C', 'T', 'P')
            InEvap.v = Aux_fn.PropCal(InEvap, 'V', 'T', 'P')
            InEvap.pr = Aux_fn.PropCal(InEvap, 'Prandtl', 'T', 'P')
            InEvap.l = Aux_fn.PropCal(InEvap, 'L', 'T', 'P')
            InEvap.d = Aux_fn.PropCal(InEvap, 'D', 'T', 'P')
            OutEvap.c = Aux_fn.PropCal(OutEvap, 'C', 'T', 'P')
            OutEvap.v = Aux_fn.PropCal(OutEvap, 'V', 'T', 'P')
            OutEvap.pr = Aux_fn.PropCal(OutEvap, 'Prandtl', 'T', 'P')
            OutEvap.l = Aux_fn.PropCal(OutEvap, 'L', 'T', 'P')
            OutEvap.d = Aux_fn.PropCal(OutEvap, 'D', 'T', 'P')
            
        return (InCond, OutCond, InEvap, OutEvap, noCond, noEvap)
    
    def OffDesign_Solver(self, InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Cycle_Inputs, Comp_Inputs, Cond_Inputs, Evap_Inputs, Outputs, noCond, noEvap):
        pcr = PropsSI('PCRIT','',0,'',0,InCond_REF.Y)
        InCond_REF.pcr = pcr
        OutCond_REF.pcr = pcr
        InEvap_REF.pcr = pcr
        OutEvap_REF.pcr = pcr
        
        comp = COMP_module()
        cond = HX_module(hx_type=Cond_Inputs.htype, cor = Cond_Inputs.cor, Inputs=Cond_Inputs)
        evap = HX_module(hx_type=Evap_Inputs.htype, cor = Evap_Inputs.cor, Inputs=Evap_Inputs)

        cond_p_ub = pcr
        if noCond == 0:
            cond_p_lb = PropsSI('P','T',OutCond.T, 'Q', 1.0, OutCond_REF.Y)
        else:
            cond_p_lb = PropsSI('P','T',InCond.T, 'Q', 1.0, InCond_REF.Y)
            
        cond_a = 1    
        
        while cond_a:
            InCond_REF.p = 0.5*(cond_p_ub + cond_p_lb)
            InCond_REF.Ts = PropsSI('T','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )
            InCond_REF.hl = PropsSI('H','P',InCond_REF.p,'Q',0.0, InCond_REF.Y )
            InCond_REF.hg = PropsSI('H','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )  
                
            if noEvap == 0:
                    evap_p_ub = PropsSI('P','T',OutEvap.T, 'Q', 1.0, OutEvap_REF.Y)
            else:
                evap_p_ub = PropsSI('P','T',InEvap.T, 'Q', 1.0, InEvap_REF.Y)
            evap_p_lb = 101300.0
            
            evap_a = 1
            
            while evap_a:
                OutEvap_REF.p = 0.5*(evap_p_ub+evap_p_lb)
                OutEvap_REF.Ts = PropsSI('T','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                OutEvap_REF.hl = PropsSI('H','P',OutEvap_REF.p,'Q',0.0, OutEvap_REF.Y )
                OutEvap_REF.hg = PropsSI('H','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                
                (OutEvap_REF, InCond_REF, Outputs.comp_W, Outputs.comp_eff_isen, Outputs.DSH, evap_a) = comp.Off(OutEvap_REF, InCond_REF, Comp_Inputs, DSH = Cycle_Inputs.DSH)
                OutCond_REF.m = InCond_REF.m
                InEvap_REF.m = OutEvap_REF.m
                
                if Cond_Inputs.htype == 'phx':                    
                    (InCond_REF, OutCond_REF, InCond, OutCond, cond_Q, cond_rho, Outputs.cond_T_pp, err_cond_index)=cond.PHX('cond',InCond_REF, OutCond_REF, InCond, OutCond, noCond)
                
                if Evap_Inputs.htype == 'phx':
                    (InEvap_REF, OutEvap_REF, InEvap, OutEvap, evap_Q, evap_rho, Outputs.evap_T_pp, err_evap_index)=evap.PHX('evap',InEvap_REF, OutEvap_REF, InEvap, OutEvap, noEvap)
                    
                if InEvap_REF.h == 0:
                    if err_evap_index == 1:
                        evap_p_lb = OutEvap_REF.p
                    else:
                        evap_p_ub = OutEvap_REF.p
                    
                    err_evap_p = 1
                else:
                    if self.evap_BC == 'q':
                        err_evap_p = (Cycle_Inputs.evap_Q - evap_Q)/Cycle_Inputs.evap_Q_cal
                    else:    
                        err_evap_p = (InEvap_REF.h - OutCond_REF.h)/OutCond_REF.h
                        
                    if err_evap_p < 0:
                        evap_p_lb = OutEvap_REF.p
                    else:
                        evap_p_ub = OutEvap_REF.p
                    
                if abs(err_evap_p) < 1.0e-3:
                    evap_a = 0
                elif evap_p_ub - evap_p_lb < 0.1:
                    evap_a = 0
            
            if self.cond_BC == 'dsc':
                cal_h = PropsSI('H','T',OutCond_REF.Ts-Cycle_Inputs.DSC,'P',OutCond_REF.p,OutCond_REF.Y)
                err_cond_p = (cal_h - OutCond_REF.h)/cal_h
            elif self.cond_BC == 'x':
                cal_h = (OutCond_REF.hg - OutCond_REF.hl)*Cycle_Inputs.cond_x+OutCond_REF.hl
                err_cond_p = (cal_h - OutCond_REF.h)/cal_h
            elif self.cond_BC == 'q':
                err_cond_p = (cond_Q - Cycle_Inputs.cond_Q)/Cycle_Inputs.cond_Q
            elif self.cond_BC == 'm':
                try:
                    dl = PropsSI('D','T',OutCond_REF.T,'P',OutCond_REF.p,OutCond_REF.Y)
                except:
                    dl = PropsSI('D','P',OutCond_REF.p,'Q', 0.0, OutCond_REF.Y)
                    
                M_ref_cal = (cond.V*cond_rho+evap.V*evap_rho+Cycle_Inputs.V_rec*dl)
                err_cond_p = (M_ref_cal - Cycle_Inputs.M_ref)/Cycle_Inputs.M_ref
            
            if err_cond_p < 0:
                cond_p_lb = InCond_REF.p
            else:
                cond_p_ub = InCond_REF.p    
            
            if abs(err_cond_p) < 1.0e-3:
                cond_a = 0
            elif cond_p_ub - cond_p_lb < 0.1:
                cond_a = 0
        try:
            dl = PropsSI('D','T',OutCond_REF.T,'P',OutCond_REF.p,OutCond_REF.Y)
        except:
            dl = PropsSI('D','P',OutCond_REF.p,'Q', 0.0, OutCond_REF.Y)
        Outputs.M_ref = (cond.V*cond_rho+evap.V*evap_rho+Cycle_Inputs.V_rec*dl)
        Outputs.cond_x = (OutCond_REF.h - OutCond_REF.hl)/(OutCond_REF.hg - OutCond_REF.hl)
        Outputs.DSC = OutCond_REF.Ts - OutCond_REF.T
        Outputs.evap_Q = evap_Q
        Outputs.cond_Q = cond_Q
        Outputs.COP_c = evap_Q/Outputs.comp_W
        Outputs.COP_h = cond_Q/Outputs.comp_W
        
        return (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Outputs)
    
    def Result_summary(self, InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Outputs):
        print('--------------------- Primary Results ---------------------')
        print('Refrigerant: %s'%(InCond_REF.Y))
        print('COP: %5.2f'%(Outputs.COP_h))
        print('mass flow rate: %5.3f[kg/s]'%(InCond_REF.m)) 
        print('InCond_REF: %5.1f[℃]/ %5.1f[bar]'%(InCond_REF.T-273.15, InCond_REF.p/1.0e5))
        print('OutCond_REF: %5.1f[℃]/ %5.1f[bar]'%(OutCond_REF.T-273.15, OutCond_REF.p/1.0e5))
        print('InEvap_REF: %5.1f[℃]/ %5.1f[bar]'%(InEvap_REF.T-273.15, InEvap_REF.p/1.0e5))
        print('OutEvap_REF: %5.1f[℃]/ %5.1f[bar]'%(OutEvap_REF.T-273.15, OutEvap_REF.p/1.0e5))
        print('Total Inventory: %5.5f[kg]'%(Outputs.M_ref))
        print('Condenser outlet quality: %5.3f'%(Outputs.cond_x))
        print('--------------------- Components Results ---------------------')
        print('condQ: %5.2f [kW]'%(Outputs.cond_Q/1.0e3))
        print('evapQ: %5.2f [kW]'%(Outputs.evap_Q/1.0e3))
        print('compW: %5.2f [kW]'%(Outputs.comp_W/1.0e3))
        print('cond_Tpp: %5.1f [℃]'%(Outputs.cond_T_pp))
        print('evap_Tpp: %5.1f [℃]'%(Outputs.evap_T_pp))
        print('cond_dp: %5.2f [-]'%((InCond_REF.p-OutCond_REF.p)/InCond_REF.p*100.0))
        print('evap_dp: %5.2f [-]'%((InEvap_REF.p-OutEvap_REF.p)/OutEvap_REF.p*100.0))
        print('comp_eff: %5.1f [℃]'%(Outputs.comp_eff_isen))
        print('DSC: %5.1f [℃]'%(Outputs.DSC))
        print('--------------------- Secondary Results ---------------------')
        print('Cond fluid: %s'%(InCond.Y))
        print('InCond: %5.1f[℃]/ %5.1f[bar]'%(InCond.T-273.15, InCond.p/1.0e5))
        print('OutCond: %5.1f[℃]/ %5.1f[bar]'%(OutCond.T-273.15, OutCond.p/1.0e5))
        print('Cond mdot: %5.3f[kg/s]'%(InCond.m))
        print('Evap fluid: %s'%(InEvap.Y))
        print('InEvap: %5.1f[℃]/ %5.1f[bar]'%(InEvap.T-273.15, InEvap.p/1.0e5))
        print('OutEvap: %5.1f[℃]/ %5.1f[bar]'%(OutEvap.T-273.15, OutEvap.p/1.0e5))
        print('Evap mdot: %5.3f[kg/s]'%(InEvap.m))
if __name__ == '__main__':
    # 대체 냉매 실험장치
    #input_ref = 'REFPROP::R32[0.76]&R125[0.07731]&CF3I[0.16269]'
    input_ref = 'REFPROP::R466A.mix'
    cycle_inputs = Cycle_Inputs()
    comp_inputs = Comp_Inputs()
    cond_inputs = PHX_Inputs()
    evap_inputs = PHX_Inputs()
    Tcond_in = 32.1+273.15
    Pcond = 101300.0
    Tevap_in = 12.1+273.15
    Pevap = 101300.0
    mcond = 2.63*PropsSI("D","T",Tcond_in,"P",Pcond,"Water")/3600
    mevap = 2.35*PropsSI("D","T",Tevap_in,"P",Pevap,"Water")/3600
    InCond = Fluid_flow(Y='water',m=mcond,T=Tcond_in, p = Pcond)
    OutCond = Fluid_flow(Y='water',m=mcond,T=0.0, p = Pcond)
    InEvap = Fluid_flow(Y='water',m=mevap, T=Tevap_in, p = Pevap)
    OutEvap = Fluid_flow(Y='water',m=mevap, T=0.0, p = Pevap)
    InCond_REF = Fluid_flow(Y=input_ref)
    OutCond_REF = Fluid_flow(Y=input_ref)
    InEvap_REF = Fluid_flow(Y=input_ref)
    OutEvap_REF = Fluid_flow(Y=input_ref)
    
    outputs = Outputs()

    cycle_inputs.layout = 'bas'
    cycle_inputs.DSH = 7.9 - PropsSI("T","P",(8.5+1.01300)*1.0e5, "Q", 1.0, "REFPROP::R466A.mix") -273.15
    dsc = PropsSI("T","P",(22.6+1.01300)*1.0e5,"Q",0.0,"REFPROP::R466A.mix")-31.9-273.15
    cycle_inputs.DSC = dsc
    cycle_inputs.cond_x = 0.01
    cycle_inputs.M_ref = 8.40615
    cycle_inputs.V_rec = 0.006
    
    comp_inputs.n_poly = 0.0
    comp_inputs.V_dis = 43.0e-6
    comp_inputs.frequency = 60
    comp_inputs.C_gap = 0.01
    comp_inputs.extra_work = 0.15

    cond_inputs.N_element = 100
    cond_inputs.N_plate = 40
    cond_inputs.thk_plate = 0.0005
    cond_inputs.thk_tot = 0.0916
    cond_inputs.L_vert = 0.455
    cond_inputs.L_width = 0.115
    cond_inputs.beta = 60
    cond_inputs.A_flow = 2.39
    cond_inputs.type = 'phx'
    cond_inputs.cor = True
    cond_inputs.mult_pri = 1.0
    cond_inputs.mult_sec = 1.0
    cond_inputs.mult_A = 1.0

    evap_inputs.N_element = 100
    evap_inputs.N_plate = 50
    evap_inputs.thk_plate = 0.0005
    evap_inputs.thk_tot = 0.120
    evap_inputs.L_vert = 0.442
    evap_inputs.L_width = 0.110
    evap_inputs.beta = 60
    evap_inputs.A_flow = 2.8
    evap_inputs.type = 'phx'
    evap_inputs.cor = True
    evap_inputs.mult_pri = 1.0
    evap_inputs.mult_sec = 1.0
    evap_inputs.mult_A = 0.4
    
    bas_off = VCHP_off('h','dsc')
    (InCond, OutCond, InEvap, OutEvap, noCond, noEvap) = bas_off.Input_Processing(InCond, OutCond, InEvap, OutEvap)
    (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = bas_off.OffDesign_Solver(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, cycle_inputs, comp_inputs, cond_inputs, evap_inputs, outputs, noCond, noEvap)
    bas_off.Result_summary(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
    # 스마트플랫폼
    '''
    input_ref = 'R410A'
    cycle_inputs = Cycle_Inputs()
    comp_inputs = Comp_Inputs()
    cond_inputs = PHX_Inputs()
    evap_inputs = PHX_Inputs()
    InEvap = Fluid_flow(Y='water',m=0.0, T=285.15, p = 101300.0)
    OutEvap = Fluid_flow(Y='water',m=0.0, T=280.15, p = 101300.0)
    InCond = Fluid_flow(Y='water',m=0.0,T=305.15, p = 101300.0)
    OutCond = Fluid_flow(Y='water',m=0.0,T=310.15, p = 101300.0)
    InEvap_REF = Fluid_flow(Y=input_ref)
    OutEvap_REF = Fluid_flow(Y=input_ref)
    InCond_REF = Fluid_flow(Y=input_ref)
    OutCond_REF = Fluid_flow(Y=input_ref)
    outputs = Outputs()

    cycle_inputs.layout = 'bas'
    cycle_inputs.DSH = 5.0
    cycle_inputs.DSC = 0.5
    cycle_inputs.cond_x = 0.01
    cycle_inputs.M_ref = 1.0
    cycle_inputs.V_rec = 0.005
    cycle_inputs.cond_Q = 4700.0
    
    comp_inputs.n_poly = 2.35
    comp_inputs.V_dis = 14.1e-6
    comp_inputs.frequency = 60*2900/3600
    comp_inputs.C_gap = 0.09
    comp_inputs.extra_work = 0.05

    cond_inputs.N_element = 150
    cond_inputs.N_plate = 40
    cond_inputs.thk_plate = 0.0005
    cond_inputs.thk_tot = 0.0644
    cond_inputs.L_vert = 0.244
    cond_inputs.L_width = 0.064
    cond_inputs.beta = 60
    cond_inputs.A_flow = 1.06
    cond_inputs.type = 'phx'
    cond_inputs.cor = True
    cond_inputs.mult_pri = 1.0
    cond_inputs.mult_sec = 1.0
    cond_inputs.mult_A = 0.85

    evap_inputs.N_element = 150
    evap_inputs.N_plate = 40
    evap_inputs.thk_plate = 0.0005
    evap_inputs.thk_tot = 0.0644
    evap_inputs.L_vert = 0.244
    evap_inputs.L_width = 0.064
    evap_inputs.beta = 60
    evap_inputs.A_flow = 1.06
    evap_inputs.type = 'phx'
    evap_inputs.cor = True
    evap_inputs.mult_pri = 1.0
    evap_inputs.mult_sec = 1.0
    evap_inputs.mult_A = 0.6
    
    bas_off = VCHP_off('h','dsc')
    (InCond, OutCond, InEvap, OutEvap, noCond, noEvap) = bas_off.Input_Processing(InCond, OutCond, InEvap, OutEvap)
    (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = bas_off.OffDesign_Solver(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, cycle_inputs, comp_inputs, cond_inputs, evap_inputs, outputs, noCond, noEvap)
    '''