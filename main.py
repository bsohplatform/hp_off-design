from HP_dataclass import *
from COMP_module import *
from HX_module import *
from CoolProp.CoolProp import PropsSI
from copy import *

class VCHP_off:
    def __init__(self, input_ref):
        self.pcr = PropsSI('PCRIT','',0,'',0,input_ref)
        self.input_ref = input_ref
    def Basic_property(self, Fluid_flow):
        Fluid_flow.c = PropsSI("C","T",Fluid_flow.T,"P",Fluid_flow.p,Fluid_flow.Y)
        Fluid_flow.v = PropsSI("V","T",Fluid_flow.T,"P",Fluid_flow.p,Fluid_flow.Y)
        Fluid_flow.pr = PropsSI("Prandtl","T",Fluid_flow.T,"P",Fluid_flow.p,Fluid_flow.Y)
        Fluid_flow.l = PropsSI("L","T",Fluid_flow.T,"P",Fluid_flow.p,Fluid_flow.Y)
        Fluid_flow.d = PropsSI("D","T",Fluid_flow.T,"P",Fluid_flow.p,Fluid_flow.Y)
            
        return (Fluid_flow)
    
    def OffDesign_outer_solver(self, InCond, OutCond, InEvap, OutEvap, Cycle_Inputs, Comp_Inputs, Cond_Inputs, Evap_Inputs):
        OutComp_REF = Fluid_flow(pcr=self.pcr, Y=self.input_ref)
        InCond_REF = Fluid_flow(pcr=self.pcr, Y=self.input_ref)
        OutCond_REF = Fluid_flow(pcr=self.pcr, Y=self.input_ref)
        InEvap_REF = Fluid_flow(pcr=self.pcr, Y=self.input_ref)
        OutEvap_REF = Fluid_flow(pcr=self.pcr, Y=self.input_ref)
        outputs = Outputs()
        
        comp = COMP_module(Comp_Inputs.mode)
        cond = HX_module(hx_type=Cond_Inputs.htype, cor = Cond_Inputs.cor, Inputs=Cond_Inputs)
        evap = HX_module(hx_type=Evap_Inputs.htype, cor = Evap_Inputs.cor, Inputs=Evap_Inputs)
        
        InCond_n = deepcopy(InCond)
        InEvap_n = deepcopy(InEvap)
        
        OutCond.T = InCond.T+5
        OutCond.p = InCond.p
        
        OutEvap.T = InEvap.T-5
        OutEvap.p = InEvap.p
        
        InCond = self.Basic_property(InCond)
        OutCond = self.Basic_property(OutCond)
        InEvap = self.Basic_property(InEvap)
        OutEvap = self.Basic_property(OutEvap)
        
        a_cond = 1
        n_cond = 0
        while a_cond:
            n_cond = n_cond+1
            a_evap = 1
            n_evap = 0
            while a_evap:
                n_evap = n_evap+1
                (InCond_n, OutCond, InEvap_n, OutEvap, OutComp_REF, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = self.OffDesign_inner_solver(InCond_n, OutCond, InEvap_n, OutEvap, OutComp_REF, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, comp, cond, evap, Cycle_Inputs, Comp_Inputs, Cond_Inputs, Evap_Inputs, outputs)
                err_evap = InEvap_n.T - InEvap.T
                OutEvap.T = OutEvap.T - err_evap
                if abs(err_evap) < 0.1:
                    a_evap = 0
                elif n_evap > 3:
                    a_evap = 0
            
            err_cond = InCond_n.T - InCond.T
            OutCond.T = OutCond.T - err_cond
            if abs(err_cond) < 0.1:
                a_cond = 0
            elif n_cond > 3:
                a_cond = 0
                
                    
        return (InCond, OutCond, InEvap, OutEvap, OutComp_REF, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
        
    def OffDesign_inner_solver(self, InCond, OutCond, InEvap, OutEvap, OutComp_REF, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, comp, cond, evap, Cycle_Inputs, Comp_Inputs, Cond_Inputs, Evap_Inputs, Outputs):
        comp_p_lb = PropsSI('P','T',InCond.T, 'Q', 1.0, InCond_REF.Y)    
        comp_p_ub = self.pcr
        
        a_m = 1
        while a_m:
            OutComp_P = 0.5*(comp_p_ub + comp_p_lb)
            
            evap_p_ub = PropsSI('P','T',InEvap.T, 'Q', 1.0, InEvap_REF.Y)    
            evap_p_lb = max(101300.0, PropsSI("PTRIPLE",InEvap_REF.Y))
            a_dsh = 1
            
            while a_dsh:
                OutEvap_P = 0.5*(evap_p_ub + evap_p_lb) 
                
                OutEvap_REF.p = OutEvap_P
                OutEvap_REF.Ts = PropsSI('T','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                OutEvap_REF.hl = PropsSI('H','P',OutEvap_REF.p,'Q',0.0, OutEvap_REF.Y )
                OutEvap_REF.hg = PropsSI('H','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                
                OutComp_REF.p = OutComp_P
                OutComp_REF.Ts = PropsSI('T','P',OutComp_REF.p,'Q',1.0, OutComp_REF.Y )
                OutComp_REF.hl = PropsSI('H','P',OutComp_REF.p,'Q',0.0, OutComp_REF.Y )
                OutComp_REF.hg = PropsSI('H','P',OutComp_REF.p,'Q',1.0, OutComp_REF.Y )  
                        
                (OutEvap_REF, OutComp_REF, Outputs.comp_W, Outputs.comp_eff_isen, Outputs.DSH, a_dsh) = comp.Off(OutEvap_REF, OutComp_REF, Comp_Inputs, Cycle_Inputs.DSH)
                
                InCond_REF.m = OutComp_REF.m
                OutCond_REF.m = InCond_REF.m
                InEvap_REF.m = OutEvap_REF.m
                
                InCond_REF = deepcopy(OutComp_REF)
                InCond_REF.p = OutComp_REF.p*(1-Cycle_Inputs.flowmeter_dp)
                InCond_REF.T = PropsSI("T","P",InCond_REF.p,"H",InCond_REF.h,InCond_REF.Y)
                        
                if Cond_Inputs.htype == 'phx':                    
                    (InCond_REF, OutCond_REF, InCond, OutCond, cond_Q, cond_rho, Outputs.cond_T_pp, err_p_cond)=cond.PHX('cond',InCond_REF, OutCond_REF, InCond, OutCond)
                elif Cond_Inputs.htypeb == 'fthx':
                    (InCond_REF, OutCond_REF, InCond, OutCond, cond_Q, cond_rho, Outputs.cond_T_pp, err_p_cond)=cond.FTHX('cond',InCond_REF, OutCond_REF, InCond, OutCond)   

                if err_p_cond == 1:
                    comp_p_lb = min(comp_p_lb*1.01,PropsSI("P","T",OutEvap.T,"Q",0.0,InCond_REF.Y))
                    print("!!응축기 온도 역전 발생!! %.2f[℃](냉매입구측), %.2f[℃](공정출구측)" %(InCond_REF.T-273.15, OutCond.T-273.15))
                    print("!!냉매 고압 압력!! %.3f[bar]" %(InCond_REF.p/1.0e5))
                    print("")
                    a_dsh = 0
                
                InEvap_REF.p = OutEvap_P/(1-Evap_Inputs.dp)
                InEvap_REF.h = OutCond_REF.h
                InEvap_REF.T = PropsSI("T","H",InEvap_REF.h,"P",InEvap_REF.p,InEvap_REF.Y)
                InEvap_REF.hl = PropsSI("H","P",InEvap_REF.p,"Q",0.0,InEvap_REF.Y)
                InEvap_REF.hg = PropsSI("H","P",InEvap_REF.p,"Q",1.0,InEvap_REF.Y)
                
                if Evap_Inputs.htype == 'phx':
                    (InEvap_REF, OutEvap_REF, InEvap, OutEvap, evap_Q, evap_rho, Outputs.evap_T_pp, err_p_evap)=evap.PHX('evap',InEvap_REF, OutEvap_REF, InEvap, OutEvap)
                elif Evap_Inputs.htype == 'fthx':
                    (InEvap_REF, OutEvap_REF, InEvap, OutEvap, evap_Q, evap_rho, Outputs.evap_T_pp, err_p_evap)=evap.FTHX('evap',InEvap_REF, OutEvap_REF, InEvap, OutEvap)
                
                if err_p_evap == 1:
                    #comp_p_ub = max(comp_p_ub*1.01, PropsSI("P","T",OutCond.T,"Q",1.0,InEvap_REF.Y))
                    comp_p_lb = 0.5*(comp_p_ub+comp_p_lb)
                    print("!!증발기 온도 역전 발생!! %.2f[℃](냉매입구측), %.2f[℃](공정출구측)" %(InEvap_REF.T-273.15, OutEvap.T-273.15))
                    print("!!냉매 고압 압력!! %.3f[bar]" %(InCond_REF.p/1.0e5))
                    print("")
                    a_dsh = 0
                
                dsh = OutEvap_REF.T - PropsSI("T","P",OutEvap_REF.p,"Q",1.0,OutEvap_REF.Y)
                err_dsh = OutEvap_REF.h - PropsSI("H","T",OutEvap_REF.Ts+Outputs.DSH,"P",OutEvap_REF.p,OutEvap_REF.Y)
                
                if err_dsh > 0:
                    evap_p_lb = 0.5*(evap_p_lb + evap_p_ub)
                else:
                    evap_p_ub = 0.5*(evap_p_lb + evap_p_ub)
                
                if abs(err_dsh/1000) < Cycle_Inputs.tol:
                    a_dsh = 0
                elif (evap_p_ub - evap_p_lb)/101300 < Cycle_Inputs.tol:
                    a_dsh = 0
                                    
                M_comp2cond = Cycle_Inputs.V_comp2cond*PropsSI("D","H",InCond_REF.h,"P",InCond_REF.p, InCond_REF.Y)
                M_cond2tev = Cycle_Inputs.V_cond2tev*PropsSI("D","H",OutCond_REF.h,"P",OutCond_REF.p, OutCond_REF.Y)
                M_tev2evap = Cycle_Inputs.V_tev2evap*PropsSI("D","H",InEvap_REF.h, "P",InEvap_REF.p,InEvap_REF.Y)
                M_evap2comp = Cycle_Inputs.V_evap2comp*PropsSI("D","H",OutEvap_REF.h, "P",OutEvap_REF.p,OutEvap_REF.Y)
                
                Outputs.M_cond = cond.V_p*cond_rho
                Outputs.M_evap = evap.V_p*evap_rho
                
                if Comp_Inputs.type == 'low':
                    d_free_volume = PropsSI("D","T",OutEvap_REF.T,"P",OutEvap_REF.p,OutEvap_REF.Y)
                else:
                    d_free_volume = PropsSI("D","T",InCond_REF.T,"P",InCond_REF.p,InCond_REF.Y)
                    
                Outputs.M_comp = Comp_Inputs.V_free*d_free_volume
                Outputs.M_oil = Comp_Inputs.V_oil*Comp_Inputs.d_oil*Comp_Inputs.frac_ref_in_oil
                Outputs.M_pipe = M_comp2cond+M_cond2tev+M_tev2evap+M_evap2comp
                Outputs.M_ref = Outputs.M_cond+Outputs.M_evap+Outputs.M_comp+Outputs.M_oil+Outputs.M_pipe
                
            if err_p_cond == 0 and err_p_evap == 0:
                if Cycle_Inputs.BC == 'm':
                    target = Cycle_Inputs.M_ref
                    param = Outputs.M_ref
                    err_m = (param - target)/target
                    
                elif Cycle_Inputs.BC == 'dsc':
                    target = PropsSI("H","T",OutCond_REF.Ts-Cycle_Inputs.DSC,"P",OutCond_REF.p,OutCond_REF.Y)
                    param = OutCond_REF.h
                    err_m = (target - param)/1000
                
                if err_m > 0:
                    comp_p_ub = 0.5*(comp_p_ub + comp_p_lb)
                else:
                    comp_p_lb = 0.5*(comp_p_ub + comp_p_lb)
                    
                if abs(err_m) < Cycle_Inputs.tol:
                    a_m = 0
                elif (comp_p_ub - comp_p_lb)/101300 < Cycle_Inputs.tol:
                    a_m = 0

        Outputs.cond_x = (OutCond_REF.h - OutCond_REF.hl)/(OutCond_REF.hg - OutCond_REF.hl)
        Outputs.DSC = OutCond_REF.Ts - OutCond_REF.T
        Outputs.evap_Q = evap_Q
        Outputs.cond_Q = cond_Q
        Outputs.COP_c = evap_Q/Outputs.comp_W
        Outputs.COP_h = cond_Q/Outputs.comp_W
        
        return (InCond, OutCond, InEvap, OutEvap, OutComp_REF, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Outputs)
    
    def Result_summary(self, InCond, OutCond, InEvap, OutEvap, OutComp_REF, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Outputs):
        print('--------------------- Primary Results ---------------------')
        print('Refrigerant: %s'%(InCond_REF.Y))
        print('COP: %5.3f'%(Outputs.COP_h))
        print('mass flow rate: %5.3f[kg/s]'%(InCond_REF.m)) 
        print('OutComp_REF: %5.3f[℃]/ %5.3f[bar]'%(OutComp_REF.T-273.15, OutComp_REF.p/1.0e5))
        print('InCond_REF: %5.3f[℃]/ %5.3f[bar]'%(InCond_REF.T-273.15, InCond_REF.p/1.0e5))
        print('OutCond_REF: %5.3f[℃]/ %5.3f[bar]'%(OutCond_REF.T-273.15, OutCond_REF.p/1.0e5))
        print('InEvap_REF: %5.3f[℃]/ %5.3f[bar]'%(InEvap_REF.T-273.15, InEvap_REF.p/1.0e5))
        print('OutEvap_REF: %5.3f[℃]/ %5.3f[bar]'%(OutEvap_REF.T-273.15, OutEvap_REF.p/1.0e5))
        print('Total Inventory: %5.5f[kg]'%(Outputs.M_ref))
        print('Condenser outlet quality: %5.3f'%(Outputs.cond_x))
        print('--------------------- Components Results ---------------------')
        print('condQ: %5.3f [kW]'%(Outputs.cond_Q/1.0e3))
        print('evapQ: %5.3f [kW]'%(Outputs.evap_Q/1.0e3))
        print('compW: %5.3f [kW]'%(Outputs.comp_W/1.0e3))
        print('cond_Tpp: %5.3f [℃]'%(Outputs.cond_T_pp))
        print('evap_Tpp: %5.3f [℃]'%(Outputs.evap_T_pp))
        print('cond_dp: %5.3f [-]'%((InCond_REF.p-OutCond_REF.p)/InCond_REF.p*100.0))
        print('evap_dp: %5.3f [-]'%((InEvap_REF.p-OutEvap_REF.p)/OutEvap_REF.p*100.0))
        print('comp_eff: %5.3f [%%]'%(Outputs.comp_eff_isen*100))
        print('DSH: %5.3f [℃]'%(OutEvap_REF.T - PropsSI("T","P",OutEvap_REF.p,"Q",1.0,OutEvap_REF.Y)))
        print('DSC: %5.3f [℃]'%(PropsSI("T","P",OutCond_REF.p,"Q",1.0,OutCond_REF.Y)-OutCond_REF.T))
        print('--------------------- Secondary Results ---------------------')
        print('Cond fluid: %s'%(InCond.Y))
        print('InCond: %5.3f[℃]/ %5.3f[bar]'%(InCond.T-273.15, InCond.p/1.0e5))
        print('OutCond: %5.3f[℃]/ %5.3f[bar]'%(OutCond.T-273.15, OutCond.p/1.0e5))
        print('Cond mdot: %5.3f[kg/s]'%(InCond.m))
        print('Evap fluid: %s'%(InEvap.Y))
        print('InEvap: %5.3f[℃]/ %5.3f[bar]'%(InEvap.T-273.15, InEvap.p/1.0e5))
        print('OutEvap: %5.3f[℃]/ %5.3f[bar]'%(OutEvap.T-273.15, OutEvap.p/1.0e5))
        print('Evap mdot: %5.3f[kg/s]'%(InEvap.m))
        
if __name__ == '__main__':
    cycle_inputs = Cycle_Inputs()
    comp_inputs = Comp_Inputs()
    cond_inputs = PHX_Inputs()
    evap_inputs = PHX_Inputs()
    
    input_ref = 'R290'
    Tcond_in = 40+273.15
    Tcond_out = 0
    Pcond = 101300.0

    Tevap_in = 10+273.15
    Tevap_out = 0
    Pevap = 101300.0

    cond_fluid = 'water'
    evap_fluid = 'water'

    mevap = 0.3/3600*PropsSI("D","T",Tevap_in,"P",Pevap,evap_fluid)
    mcond = 0.35/3600*PropsSI("D","T",Tcond_in,"P",Pcond,cond_fluid)

    InCond = Fluid_flow(Y=cond_fluid,m=mcond,T=Tcond_in, p = Pcond)
    OutCond = Fluid_flow(Y=cond_fluid,m=mcond,T=Tcond_out, p = Pcond)
    InEvap = Fluid_flow(Y=evap_fluid,m=mevap, T=Tevap_in, p = Pevap)
    OutEvap = Fluid_flow(Y=evap_fluid,m=mevap, T=Tevap_out, p = Pevap)
    
    cycle_inputs.BC = 'm'
    cycle_inputs.DSH = 5
    cycle_inputs.layout = 'bas'
    cycle_inputs.M_ref = 0.25
    cycle_inputs.tol = 1.0e-2
    cycle_inputs.V_comp2cond = 0.0
    cycle_inputs.V_cond2tev = 0.0
    cycle_inputs.V_tev2evap = 0.0
    cycle_inputs.V_evap2comp = 0.0

    comp_inputs.mode = 'poly'
    comp_inputs.n_poly = 1.0
    comp_inputs.V_dis = 22.2e-6
    comp_inputs.V_free = 0.0
    comp_inputs.frequency = 50.0
    comp_inputs.C_gap = 0.05
    comp_inputs.extra_work = 0.42
    comp_inputs.V_oil = 0.0
    comp_inputs.frac_ref_in_oil = 0.2

    evap_inputs.N_element = 100
    evap_inputs.N_plate = 18
    evap_inputs.thk_plate = 0.00014
    evap_inputs.pitch_p = 0.00234
    evap_inputs.pitch_s = 0.00234
    evap_inputs.crg_pitch_p = 3*evap_inputs.pitch_p
    evap_inputs.crg_pitch_s = 3*evap_inputs.pitch_s
    evap_inputs.enlargement_p = 1.13
    evap_inputs.enlargement_s = 1.13
    evap_inputs.Nch_p = 8
    evap_inputs.Nch_s = 9
    evap_inputs.L_vert = 0.479
    evap_inputs.L_width = 0.117
    evap_inputs.beta = 60
    evap_inputs.A_flow = 1.01
    evap_inputs.type = 'phx'
    evap_inputs.cor = True
    evap_inputs.mult_pri = 0.5
    evap_inputs.mult_sec = 0.5
    evap_inputs.mult_A = 0.9
    evap_inputs.k = 237.0

    cond_inputs.N_element = 100
    cond_inputs.N_plate = 20
    cond_inputs.thk_plate = 0.00014
    cond_inputs.pitch_p = 0.00234
    cond_inputs.pitch_s = 0.00234
    cond_inputs.crg_pitch_p = 3*cond_inputs.pitch_p
    cond_inputs.crg_pitch_s = 3*cond_inputs.pitch_s
    cond_inputs.enlargement_p = 1.17
    cond_inputs.enlargement_s = 1.17
    cond_inputs.Nch_p = 9
    cond_inputs.Nch_s = 10
    cond_inputs.L_vert = 0.479
    cond_inputs.L_width = 0.117
    cond_inputs.beta = 60
    cond_inputs.A_flow = 1.13
    cond_inputs.type = 'phx'
    cond_inputs.cor = True
    cond_inputs.mult_pri = 0.7
    cond_inputs.mult_sec = 0.7
    cond_inputs.mult_A = 0.9
    cond_inputs.k = 237.0
    
    outputs = Outputs()
    bas_off = VCHP_off(input_ref)
    (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = bas_off.OffDesign_outer_solver(InCond, OutCond, InEvap, OutEvap, cycle_inputs, comp_inputs, cond_inputs, evap_inputs)
    bas_off.Result_summary(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)

    print('충진량 분포----전체: %.2f, 응축기: %.2f [g], 증발기: %.2f [g], 압축기: %.2f [g], 오일: %.2f [g], 배관: %.2f [g]' %(outputs.M_ref*1000, outputs.M_cond*1000, outputs.M_evap*1000, outputs.M_comp*1000, outputs.M_oil*1000, outputs.M_pipe*1000))
    print('충진량: %.2f [g/kW]' %(outputs.M_ref*1.0e6/outputs.cond_Q))