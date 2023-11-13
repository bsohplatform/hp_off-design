from math import pi, sin

class PHX:
    @staticmethod
    def cor_cond_2ph(p, pcr, qq, G, Re, pr, hl, hg, ll, vl, Dh):
        # Yan (1998)
        Bo = qq/G/(hg-hl)
        htc = 4.118*(ll/Dh)*pow(Re,0.4)*pow(pr, 0.333333)
        f = 94.75*pow(p/pcr,0.8)*pow(Bo,0.5)*pow(G*Dh/vl,-0.4)*pow(Re, -0.0467)
        '''
        if G < 35: # Han (2003), Experiments on the characteristics of evaporation of R410A in brazed plate heat exchangers with different geometric configurations
            Ge1 = 11.22*pow(crg_pitch/Dh, -2.83)*pow(pi/2-pi*beta/180, -4.5)
            Ge2 = 0.35*pow(crg_pitch/Dh, 0.23)*pow(pi/2-pi*beta/180, 1.48)
            Ge3 = 3521.1*pow(crg_pitch/Dh, 4.17)*pow(pi/2-pi*beta/180,-7.75)
            Ge4 = -1.024*pow(crg_pitch/Dh,0.0925)*pow(pi/2-pi*beta/180,-1.3)
            htc = Ge1*(ll/Dh)*pow(Re,Ge2)*pow(pr, 0.333333)
            f = Ge3*pow(Re,Ge4)
        elif G > 60: # Yan (1999), Condensation heat transfer and pressure drop of refrigerant R-134a in a plate heat exchanger
            Bo = qq/G/(hg-hl)
            htc = 4.118*(ll/Dh)*pow(Re,0.4)*pow(pr, 0.333333)
            f = 94.75*pow(p/pcr,0.8)*pow(Bo,0.5)*pow(G*Dh/vl,-0.4)*pow(Re, -0.0467)
        else: # Interpolation
            Ge1 = 11.22*pow(crg_pitch/Dh, -2.83)*pow(pi/2-pi*beta/180, -4.5)
            Ge2 = 0.35*pow(crg_pitch/Dh, 0.23)*pow(pi/2-pi*beta/180, 1.48)
            Ge3 = 3521.1*pow(crg_pitch/Dh, 4.17)*pow(pi/2-pi*beta/180,-7.75)
            Ge4 = -1.024*pow(crg_pitch/Dh,0.0925)*pow(pi/2-pi*beta/180,-1.3)
            Bo = qq/G/(hg-hl)
            htc = (G-35)/25*4.118*(ll/Dh)*pow(Re,0.4)*pow(pr, 0.333333) + (60-G)/25*Ge1*(ll/Dh)*pow(Re,Ge2)*pow(pr, 0.333333)
            f = (G-35)/25*94.75*pow(p/pcr,0.8)*pow(Bo,0.5)*pow(G*Dh/vl,-0.4)*pow(Re, -0.0467) + (60-G)/25*Ge3*pow(Re,Ge4)
        ''' 
        return (htc, f)
    
    @staticmethod
    def cor_evap_2ph(qq, x, G, Geq, hl, hg, dl, dg, ll, vl, vg, i, Dh, beta):
        g = 9.80665
        # Almalfi (2016), Flow boiling and frictional pressure gradients in plate heat exchangers. Part 2: Comparison of literature methods to database and new prediction methods
        Bd = (dl-dg)*g*Dh**2/i
        Bo = qq/G/(hg-hl)
        beta_max = 70.0
        dm = 1/(x/dg + (1-x)/dl);
        if Bd < 4:
            htc = 982.0*(ll/Dh)*pow(beta/beta_max,1.101)*pow(Geq**2*Dh/dm/i,0.315)*pow(dl/dg,-0.224)*pow(Bo,0.32)
        else:
            htc = 18.495*pow(beta/beta_max,0.248)*pow(x*G*Dh/vg,0.135)*pow(G*Dh/vl,0.351)*pow(dl/dg,0.223)*pow(Bd,0.235)*pow(Bo,0.198)
            htc = htc*(ll/Dh)
        C = 2.125*pow(beta/beta_max,9.993)+0.955
        f = 4*C*15.698*pow(Geq**2*Dh/dm/i,-0.475)*pow((dl-dg)*g*Dh**2/i,0.255)*pow(dl/dg,-0.571)
        
        return (htc, f)
    
    @staticmethod
    def cor_1ph(Re, pr, l, v_ratio, Dh, beta, enlargement):
        # Wanniarachchi
        beta1 = max(min(62, beta),20)
        p = 0.00423*beta1+0.0000223*beta1**2
        ft = 46.6*pow(beta1,-1.08)*pow(enlargement,1+p)*pow(Re,-p)
        f1 = 1774*pow(beta1,-1.026)*enlargement**2/Re
        m = 0.646+0.0011*beta1
        Nut = 12.6*pow(beta1,-1.142)*pow(enlargement,1-m)*pow(Re,m)
        Nu1 = 3.65*pow(beta1,-0.455)*pow(enlargement,0.661)*pow(Re,0.339)
        htc = (l/Dh)*pow(pow(Nu1,3)+pow(Nut,3),1/3)*pow(pr,1/3)*pow(v_ratio,0.17)
        
            
        f = pow(pow(f1,3)+pow(ft,3),1/3)
        '''
        if Re > 1.0e3: # Muley (1999), Experimental Study of Turbulent Flow Heat Transfer and Pressure Drop in a Plate Heat Exchanger With Chevron Plates
            C1 = 0.2668-0.006967*beta+7.244e-5*beta**2
            D1 = 20.78-50.94*enlargement+41.16*enlargement**2-10.51*pow(enlargement,3)
            P1 = 0.728+0.0543*sin(pi*beta/45.0+3.7)
            C2 = 2.917-0.1277*beta+2.016e-3*pow(beta,2)
            D2 = 5.474-19.02*enlargement+18.93*enlargement**2-5.341*pow(enlargement,3)
            P2 = 0.2+0.0577*sin(pi*beta/45+2.1)
            htc = (l/Dh)*C1*D1*pow(Re,P1)*pow(pr,0.333333)*pow(v_ratio, 0.14)
            f = C2*D2*pow(Re, P2)
        elif Re < 400: # Muley (1999), Enhanced Heat Transfer Characteristics of Viscous Liquid Flows in a Chevron Plate Heat Exchanger
            htc = (l/Dh)*1.6774*pow(Dh/L, 0.333333)*pow(beta/30,0.38)*pow(Re, 0.5)*pow(pr, 0.333333)*pow(v_ratio, 0.14)
            f = pow(pow(30.20/Re,5)+pow(6.28/pow(Re,0.5),5),0.2)*pow(beta/30,0.83)
        else:
            C1 = 0.2668-0.006967*beta+7.244e-5*beta**2
            D1 = 20.78-50.94*enlargement+41.16*enlargement**2-10.51*pow(enlargement,3)
            P1 = 0.728+0.0543*sin(pi*beta/45.0+3.7)
            C2 = 2.917-0.1277*beta+2.016e-3*pow(beta,2)
            D2 = 5.474-19.02*enlargement+18.93*enlargement**2-5.341*pow(enlargement,3)
            P2 = 0.2+0.0577*sin(pi*beta/45+2.1)
            htc = (Re - 400)/600*(l/Dh)*1.6774*pow(Dh/L, 0.333333)*pow(beta/30,0.38)*pow(Re, 0.5)*pow(pr, 0.333333)*pow(v_ratio, 0.14)+\
            (1000 - Re)/600*(l/Dh)*1.6774*pow(Dh/L, 0.333333)*pow(beta/30,0.38)*pow(Re, 0.5)*pow(pr, 0.333333)*pow(v_ratio, 0.14)
            
            f = (Re - 400)/600*C2*D2*pow(Re, P2)+(1000 - Re)/600*pow(pow(30.20/Re,5)+pow(6.28/pow(Re,0.5),5),0.2)*pow(beta/30,0.83)'''
            
        return (htc, f)
    
    @staticmethod
    def htc_water(Re, pr, beta):
        # Chisholm and Wanniarachchi
        htc = 0.724*pow(beta/30,0.646)*pow(Re, 0.583)*pow(pr, 1/3)
        
        return htc

class FTHX:
    @staticmethod
    def htc_dittus(Re, pr, Dh, l, purpose):
        if purpose == 'cond':
            n = 0.3
        else:
            n = 0.4
        
        htc = (l/Dh)*0.023*Re**0.8*pr*n
        
        return htc