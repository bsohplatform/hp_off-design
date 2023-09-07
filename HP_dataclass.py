from dataclasses import dataclass, field
from CoolProp.CoolProp import PropsSI

@dataclass
class Fluid_flow:
    Y: str
    m: float=0.0
    T: float=0.0
    Ts: float=0.0
    p: float=0.0
    h: float=0.0
    hl: float=0.0
    hg: float=0.0
    Tcr: float=0.0
    pcr: float=0.0
    c: float=0.0
    v: float=0.0
    pr: float=0.0
    l: float=0.0
    d: float=0.0
    
@dataclass
class Cycle_Inputs:
    layout: str = 'bas'
    DSH: float = 0.0
    DSC: float = 0.0
    cond_x: float = 0.0
    M_ref: float = 0.0
    cond_Q: float = 0.0
    evap_Q: float = 0.0
    V_rec: float = 0.0
@dataclass
class Comp_Inputs:
    V_dis: float = 0.0 # Displacement volume
    C_gap: float = 0.0 # Clearance factor (clearance/displacement)
    n_poly: float = 0.0 # Polytropic number
    eff_mech: float = 1.0 # mechanical_efficiency
    frequency: float = 60.0 # Compressor frequency [Hz]
    extra_work: float = 0.0

@dataclass
class PHX_Inputs:
    N_element: int = 30
    N_plate: int = 0 # Number of Plates
    thk_plate: float = 0.0 # Single plate thickness
    thk_tot: float = 0.0 # Total thickness of PHE
    L_vert: float = 0.0 # Vertical length of PHE (Center to center of inlet and outlet ports)
    L_width: float = 0.0 # total horizontal length of PHE
    beta: float = 0.0 # chevron angle
    crg_pitch: float = 0.0 # corrugation pitch
    crg_depth: float = 0.0 # corrugation depth
    A_flow: float = 0.0
    UA: float = 0.0
    dp: float = 0.0
    mdot_nominal:float = 0.0
    cor: bool = False
    htype: str = 'phx'
    qloss_frac: float = 0.0
    mult_pri: float = 1.0
    mult_sec: float = 1.0
    mult_A: float = 1.0
@dataclass
class Outputs:
    comp_W: float = 0.0
    evap_Q: float = 0.0
    cond_Q: float = 0.0
    DSC: float = 0.0
    cond_x: float = 0.0
    M_ref: float = 0.0
    COP_c: float = 0.0
    COP_h: float = 0.0
    cond_T_pp: float = 0.0
    evap_T_pp: float = 0.0
    
    comp_eff_isen: float = 0.0


class Aux_fn:
    @staticmethod
    def PropCal(Fluid_flow, outpar: str, par1: str, par2: str):
        prop = {'T':Fluid_flow.T, 'P':Fluid_flow.p, 'H':Fluid_flow.h, 
                         'D':Fluid_flow.d, 'C': Fluid_flow.c, 'L': Fluid_flow.l}
                
        outval = PropsSI(outpar,par1,prop[par1],par2,prop[par2],Fluid_flow.Y)
        
        return outval
    

if __name__ == '__main__':
    cond = Fluid_flow(Y='Water', T=300, p = 101300)
    
    cond.d = Aux_fn.PropCal(cond, 'D', 'T', 'P')
    print(cond.d)