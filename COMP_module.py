from CoolProp.CoolProp import PropsSI
from HP_dataclass import *

class COMP_module:
    def __init__(self, mode='poly'):
        self.mode = mode
        
    def Off(self, primary_in, primary_out, Inputs, DSH=5.0):
        if self.mode == 'poly':
            primary_in.T =  DSH + primary_in.Ts
            primary_in.d = Aux_fn.PropCal(primary_in, 'D', 'T', 'P')
            primary_in.h = Aux_fn.PropCal(primary_in, 'H', 'T', 'P')
            Z = Aux_fn.PropCal(primary_in, 'Z', 'T', 'P')
            s = Aux_fn.PropCal(primary_in, 'S', 'T', 'P')
            
            n_comp = Aux_fn.PropCal(primary_in,'Cpmass','T','P')/Aux_fn.PropCal(primary_in,'Cvmass','T','P') if Inputs.n_poly == 0 else Inputs.n_poly
            n_comp = max(n_comp, 1.0)
            V_comp = Inputs.V_dis
            f_comp = Inputs.frequency
            C_comp = Inputs.C_gap
            
            eff_vol = 1 + C_comp - C_comp*pow(primary_out.p/primary_in.p,1/n_comp)
            primary_in.m = eff_vol*primary_in.d*V_comp*f_comp
            primary_out.m = primary_in.m
            
            w_comp = (n_comp/(n_comp-1))*Z*(primary_in.p/primary_in.d)*(pow(primary_out.p/primary_in.p,(n_comp-1)/n_comp) - 1); 
            w_comp = w_comp*(1+Inputs.extra_work)
            comp_W = primary_in.m*w_comp/Inputs.eff_mech
            
            primary_out.h = primary_in.h + w_comp
            primary_out.T = Aux_fn.PropCal(primary_out, 'T', 'H', 'P')
            h_comp_out_ideal = PropsSI('H','P',primary_out.p,'S',s,primary_out.Y)
            comp_eff_isen = (h_comp_out_ideal - primary_in.h)/(primary_out.h - primary_in.h)
            
            try:
                primary_out.Ts = PropsSI('T','P',primary_out.p, 'Q', 1.0, primary_out.Y)
            except:
                primary_out.Ts = PropsSI('TCRIT','',0,'',0,primary_out.Y)
            if abs(primary_out.T - primary_out.Ts) < 0.1:
                DSH = DSH + 1.0
                a = 0
            else:
                a = 1
            
            return (primary_in, primary_out, comp_W, comp_eff_isen, DSH, a)
            

if __name__ == '__main__':
    comp_in = Fluid_flow(Y='R410A', T=290.0, p = 1.0e6)
    comp_in.Ts = PropsSI('T','P',1.0e6,'Q',1.0,comp_in.Y)
    comp_out = Fluid_flow(Y='R410A',p = 3.0e6)
    inputs = Comp_Inputs()
    inputs.comp_n_poly = 2.3
    inputs.comp_V_dis = 9.3e-6
    inputs.comp_frequency = 60
    inputs.comp_C_gap = 0.05
    comp = COMP_module()
    (comp_in, comp_out, comp_W, comp_eff_isen, DSH, cond_a) = comp.Off(comp_in, comp_out, inputs, DSH = 5.0)
    
    print(comp_W)