from CoolProp.CoolProp import PropsSI
from HP_dataclass import *
import joblib
from math import *

class COMP_module:
    def __init__(self, mode):
        self.mode = mode
        
    def Off(self, primary_in, primary_out, Inputs, DSH):
        
        primary_in.T =  DSH + primary_in.Ts
        primary_in.d = PropsSI("D","T",primary_in.T,"P",primary_in.p,primary_in.Y)
        primary_in.h = PropsSI("H","T",primary_in.T,"P",primary_in.p,primary_in.Y)
        Z = PropsSI("Z","T",primary_in.T,"P",primary_in.p,primary_in.Y)
        s = PropsSI("S","T",primary_in.T,"P",primary_in.p,primary_in.Y)
        
        n_comp = PropsSI("Cpmass","T",primary_in.T,"P",primary_in.p,primary_in.Y)/PropsSI("Cvmass","T",primary_in.T,"P",primary_in.p,primary_in.Y)
        
        V_comp = Inputs.V_dis
        f_comp = Inputs.frequency
        C_comp = Inputs.C_gap
        
        h_comp_out_ideal = PropsSI('H','P',primary_out.p,'S',s,primary_out.Y)
        
        if self.mode == 'map':
            pr_now = primary_out.p/primary_in.p
            if len(Inputs.rpm_list) == 1:
                comp_eff_isen = self.Interpol_map(Inputs.pr_mat, Inputs.eff_mat, pr_now)
                primary_in.m = self.Interpol_map(Inputs.pr_mat, Inputs.mdot_mat, pr_now)
            else:
                # 이건 나중에 RPM 여러개일 때 개선
                a = 1
                
            primary_out.h = (h_comp_out_ideal - primary_in.h)/comp_eff_isen + primary_in.h
            primary_out.T = PropsSI("T","H",primary_out.h, "P",primary_out.p,primary_out.Y)
            w_comp = (primary_out.h - primary_in.h)
            primary_out.m = primary_in.m
            comp_W = primary_in.m*w_comp/Inputs.eff_mech
        else:
            if self.mode == 'poly':
                n_comp = n_comp*Inputs.n_poly
                w_comp = (n_comp/(n_comp-1))*Z*(primary_in.p/primary_in.d)*(pow(primary_out.p/primary_in.p,(n_comp-1)/n_comp) - 1)
                w_comp = w_comp*(1+Inputs.extra_work)
                
                primary_out.h = primary_in.h + w_comp
                primary_out.T = PropsSI("T","H",primary_out.h,"P",primary_out.p,primary_out.Y)
                comp_eff_isen = (h_comp_out_ideal - primary_in.h)/(primary_out.h - primary_in.h)
                
            elif self.mode == 'pkl':
                comp_eff_isen = joblib.load(Inputs.comp_pkl).predict([[DSH, primary_in.p, primary_out.p]])
                comp_eff_isen = float(comp_eff_isen)
                
                primary_out.h = (h_comp_out_ideal - primary_in.h)/comp_eff_isen + primary_in.h
                primary_out.T = PropsSI("T","H",primary_out.h, "P",primary_out.p,primary_out.Y)
                
            eff_vol = 1 + C_comp - C_comp*pow(primary_out.p/primary_in.p,1/n_comp)
            primary_in.m = eff_vol*primary_in.d*V_comp*f_comp
            primary_out.m = primary_in.m
            comp_W = primary_in.m*w_comp/Inputs.eff_mech
            
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
    
    def Interpol_map(self, x_list, y_list, x_now):
        if x_now < x_list[0]:
            y_now = y_list[0]
            #print('범위를 벗어났습니다.(Inner)')
        elif x_now > x_list[-1]:
            y_now = y_list[-1]
            #print('범위를 벗어났습니다.(Outer)')
        else:
            x_idx = [i for i in range(len(x_list)-1) if x_list[i] <= x_now < x_list[i+1]]
            y_now = (y_list[x_idx[0]+1]-y_list[x_idx[0]])/(x_list[x_idx[0]+1]-x_list[x_idx[0]])*(x_now-x_list[x_idx[0]])+y_list[x_idx[0]]
        return y_now
        
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