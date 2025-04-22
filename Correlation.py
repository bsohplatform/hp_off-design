from math import pi, sin, log10
from CoolProp.CoolProp import PropsSI, get_aliases

class PHX:
    @staticmethod
    def TP_state_prop(h, p, fluid):
        T = PropsSI('T','P',p,'H', h, fluid)
        hl = PropsSI('H','P',p,'Q',0.0, fluid)
        hg = PropsSI('H','P',p,'Q',1.0, fluid)
        x = max(min((h-hl)/(hg-hl),1),0)
        
        return (T, hl, hg, x)
    
    @staticmethod
    def TP_transport_prop(G_ref, x, p, fluid):
        dl = PropsSI('D','P',p,'Q',0.0, fluid)
        dg = PropsSI('D','P',p,'Q',1.0, fluid)
        d = 1/(x/dg+(1-x)/dl)
        cl = PropsSI('C','P',p,'Q',0.0, fluid)
        cg = PropsSI('C','P',p,'Q',1.0, fluid)
        vl = PropsSI('V','P',p,'Q',0.0, fluid)
        vg = PropsSI('V','P',p,'Q',1.0, fluid)
        ll = PropsSI('L','P',p,'Q',0.0, fluid)
        lg = PropsSI('L','P',p,'Q',1.0, fluid)
        i = PropsSI('I','Q',x,'P',p,fluid)
        
        # Lockhart-Martinelli correction mass flux
        G = G_ref*((1-x)+x*(dl/dg)**0.5)
        
        return (G, d, dl, dg, cl, cg, vl, vg, ll, lg, i)
    
    @staticmethod
    def SP_transport_prop(T, p, fluid):
        d = PropsSI('D','T',T,'P',p,fluid)
        c = PropsSI('C','T',T,'P',p,fluid)
        v = PropsSI('V','T',T,'P',p,fluid)
        l = PropsSI('L','T',T,'P',p,fluid)
        
        return (d, c, v, l)
        
    @staticmethod
    def cor_2ph(purpose, Tp, Tw, p, pcr, qq, x, G, Geq, hl, hg, cl, cg, dl, dg, dm, ll, lg, vl, vg, i, Dh, beta, L, enlargement, fluid, cor):
        # Define parameters for correlations
        grav = 9.80665
        Re = Geq*Dh/vl
        Rel = Geq*Dh/vl*(1-x)
        Prl = vl*cl/ll
        Xtt = ((1-x)/x)**0.9*(dg/dl)**0.5*(vl/vg)**0.1
        phil = (1+12/Xtt+1/Xtt**2)**0.5
        Ga = dl*(dl-dg)*grav*Dh**3/vl**2
        P_red = p/pcr
        w = -log10(P_red)
        Bo = qq/Geq/(hg-hl)
        Fr = Geq**2/(dm**2*grav*Dh)
        Co = (dg/dl)**0.5*((1-x)/x)**0.8
        Bd = (dl-dg)*grav*Dh**2/i
        We = Geq**2*Dh/dm/i
        beta_max = 70.0
        
        if purpose == 'cond':
            f = 94.75*P_red**0.8*Bo**0.5*(G*Dh/vl)**-0.4*Re**-0.0467
            
            if cor == 'yan':
                # Yan (1999), Condensation heat transfer and pressure drop of refrigerant R-134a in a plate heat exchanger
                htc = 4.118*(ll/Dh)*Re**0.4*Prl**0.333333
            elif 'palmer' in cor:
                if cor == 'palmer_r32': # palmer (2000), Evaporation and Condensation Heat Transfer Performance of Flammable Refrigerants   in a Brazed Plate Heat Exchanger
                    a = 1
                    b = 0.298
                    c = 0
                    d = 0.346
                    e = 0
                    f = 1.5
                    g = 1.5
                else:
                    a = 1
                    b = 0.387
                    c = 0.0824
                    d = 0.346
                    e = 0
                    f = 1.5
                    g = 1.5
                Nul = 0.16*Rel**0.89*Prl**0.3 # cooling
                htc = (ll/Dh)*a*Nul**b*phil**c*Ga**d*Xtt**e*P_red**f*w**g
            else: # longo (2015b), A new computational procedure for refrigerant condensation inside herringbone-type Brazed Plate Heat Exchangers
                if Re < 1600:
                    htc = 0.943*((ll**3*dl**2*grav*(hg-hl))/(vl*(Tp-Tw)*L))**0.25
                else:
                    htph = 1.875*enlargement*(ll/Dh)*Re**0.445*Prl**0.333333
                    hsph = 0.2267*(lg/Dh)*(G*Dh/vg)**0.631*(vg*cg/lg)**0.333333
                    Tsat = PropsSI("T","P",p,"Q",0.0,fluid)
                    htc = htph+(Tp-Tsat)/(Tsat-Tw)*(hsph+cg*(qq*(1-x))/(hg-hl))
                    
        else:
            C = 2.125*(beta/beta_max)**9.993+0.955
            f = 4*C*15.698*We**-0.475*Bd**0.255*(dl/dg)**-0.571
                
            if cor == 'yan':
                htc = (ll/Dh)*1.926*Prl**0.333333*Bo**0.3*Re**0.5
                if Re < 6000:
                    f = 6.947e5*Re**-1.109*(G*Dh/vl)**-0.5
                else:
                    f = 31.21*Re**0.04557*(G*Dh/vl)**-0.5
            elif 'palmer' in cor:
                if cor == 'palmer_r32':
                    a = 1.0
                    b = 0.42
                    c = 0
                    d = 0.088
                    e = 1.5
                    f = 1.5
                    g = 0
                    h = 1.5
                    Nul = 0.16*Rel**0.89*Prl**0.4 # heating
                    htc = (ll/Dh)*a*Nul**b*Bo**c*Fr**d*w**e*Co**f*Xtt**g*(PropsSI("M",fluid)*1000)**h
                else:
                    a = 2.7
                    b = 0.55
                    c = 0.5
                    htc = (ll/Dh)*a*Rel**b*Prl**c
            elif cor == 'amalfi':
                # Almalfi (2016), Flow boiling and frictional pressure gradients in plate heat exchangers. Part 2: Comparison of literature methods to database and new prediction methods    
                
                if Bd < 4:
                    htc = 982.0*(ll/Dh)*(beta/beta_max)**1.101*We**0.315*Bo**0.32*(dl/dg)**-0.224
                else:
                    htc = 18.495*(ll/Dh)*(beta/beta_max)**0.248*(x*Geq*Dh/vg)**0.135*(Geq*Dh/vl)**0.351*Bd**0.235*Bo**0.198*(dl/dg)**-0.223
    
                
            else: # longo (2015), A new model for refrigerant boiling inside Brazed Plate Heat Exchangers (BPHEs)
                ref_list = ["Methane","Ethene","Ethane","Propene","Propane","n-Butane","i-Butane","n-Pentane","i-Pentane","n-Hexane","Cyclohexane","n-Heptane","Benzene","Toluene","Diphenyl","Methanol","Ethanol","n-Propanol","i-Propanol","n-Butanol","i-Butanol","2-Butanol","Ethanal","Acetone","CO2","R23","R32","R125","R134a","R143a","R152a","R1234yf","R227ea","R507","RC318","R14","R123","R11","R12","R13","R13B1","R22","R113","R114","R115","R502","CH3Cl","CCl4","SF6","H2O","NH3","O2","N2","Ar","Ne","H2","He"]
                h0_list = [7.13e3,5.29e3,4.94e3,4.25e3,4.12e3,3.7e3,3.7e3,3.3e3,3.33e3,3.11e3,2.97e3,2.89e3,3.03e3,2.96e3,2.47e3,498e3,4.44e3,3.52e3,3.97e3,3.35e3,3.41e3,3.33e3,3.71e3,3.48e3,5.35e3,5.03e3,4.7e3,4.43e3,4.26e3,4.5e3,4.21e3,387e3,3.93e3,4.46e3,3.72e3,5.98e3,3.53e3,3.45e3,4e3,4.77e3,4.17e3,4.37e3,3.2e3,3.61e3,4.06e3,4.1e3,4.18e3,2.97e3,4.53e3,4.13e3,4.9e3,8.37e3,9.12e3,8.54e3,21.9e3,21.5e3,75.8e3]
                h0_dict = dict(zip(ref_list,h0_list))
                if fluid.startswith('REFPROP::'):
                    key = fluid.split("::")[-1]
                    if key in h0_dict:
                        h0 = h0_dict[key]
                    else:
                        print(f"'{key}' not found in h0 list")
                        return None
                
                aliases = get_aliases(fluid)
                for name in aliases:
                    if name in h0_dict:
                        h0 = h0_dict[name]
                        break
                        
                    else:
                        print(f"'{fluid}' not found in h0 list")
                        return None

                hcb = 0.122*enlargement*(ll/Dh)*Re**0.8*Prl**0.333333
                C_ra = 1 # (Ra/0.4)**0.1333
                Fp = 1.2*P_red**0.27+(2.5+1/(1-P_red))*P_red
                hnb = 0.58*enlargement*h0*C_ra*Fp*(qq/20000)**0.467
                htc = max(hcb,hnb)
        
        return (htc, f)
    
    @staticmethod
    def cor_1ph(G, c, d, l, v, v_ratio, Dh, beta, depth, pitch, enlargement, cor):
        Re = G*Dh/v
        Pr = v*c/l
        
        if cor == 'yan':
            # Yan (1999), Evaporation Heat Transfer and Pressure Drop of Refrigerant R-134a in a Plate Heat Exchanger
            # Water-to-Water
            htc = 0.2121*(l/Dh)*Re**0.78*Pr**0.333333*v_ratio**0.14    
        elif cor == 'muley':
            # Muley (1999), Experimental Study of Turbulent Flow Heat Transfer and Pressure Drop in a Plate Heat Exchanger with Chevron plates
            # Water-to-Water
            C1 = 0.2668-0.006967*beta+7.244e-5*beta**2
            D1 = 20.78-50.94*enlargement+41.16*enlargement**2-10.51*enlargement**3
            P1 = 0.728+0.0543*sin(pi*beta/45.0+3.7)
            C2 = 2.917-0.1277*beta+2.016e-3*beta**2
            D2 = 5.474-19.02*enlargement+18.93*enlargement**2-5.341*enlargement**3
            P2 = 0.2+0.0577*sin(pi*beta/45+2.1)
            htc = (l/Dh)*C1*D1*Re**P1*Pr**0.333333*v_ratio**0.14
            f = C2*D2*Re**-P2
        elif cor == 'wan':
            # Wanniarachchi (1995), Approximate Correlations for Chevron-Type Plate Heat Exchangers
            beta1 = max(min(62, beta),20)
            m = 0.646+0.0011*beta1
            Nut = 12.6*beta1**-1.142*enlargement**(1-m)*Re**m
            Nu1 = 3.65*beta1**-0.455*enlargement**0.661*Re**0.339
            htc = (l/Dh)*(Nu1**3+Nut**3)**0.333333*Pr**0.333333*v_ratio**0.17
            
            p = 0.00423*beta1+0.0000223*beta1**2
            f1 = 1774*beta1**-1.026*enlargement**2/Re
            ft = 46.6*beta1**-1.08*enlargement**(1+p)*Re**-p
            f = (f1**3+ft**3)**0.333333
        else: # Yang (2017), Heat transfer correlations for single-phase flow in plate heat exchangers based on experimental data
            C2 = 2.917-0.1277*beta+2.016e-3*beta**2
            D2 = 5.474-19.02*enlargement+18.93*enlargement**2-5.341*enlargement**3
            P2 = 0.2+0.0577*sin(pi*beta/45+2.1)
            C3 = -1.343e-4*beta**2+1.808e-2*beta-0.0075
            P3 = -7.956e-5*beta**2+9.687e-3*beta+0.3155
            gamma = 2*depth/pitch
            htc = (l/Dh)*C3*Re**P3*Re**(enlargement/beta)*Re**(gamma/beta)*Pr**0.333333*v_ratio**0.14
            f = C2*D2*Re**-P2
            
        return (htc, f)
        
    
class FTHX:
    @staticmethod
    def htc_dittus(Re, pr, Dh, l, purpose):
        if purpose == 'cond':
            n = 0.3
        else:
            n = 0.4
        
        htc = (l/Dh)*0.023*Re**0.8*pr*n
        
        return htc