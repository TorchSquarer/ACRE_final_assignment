import numpy as np

class model():
    
    def __init__(self)-> list[np.matrix, np.array]:
        # --------------------------------------------------------
        # Parameters
        # --------------------------------------------------------
        # Flow rates [m3/s]
        F_aq  = 1e-4
        F_org = 5e-5
        F_R   = F_aq + F_org

        # Phase volumes [m3]
        V_aq  = 1e-2
        V_org = 2*V_aq
        V_R   = V_aq + V_org

        # Feed concentration [mol/m3]
        C_Glu = 500.0

        # Rate constants [1/s], pseudo-first-order
        k1 = 1e-3    # Glucose -> Fructose
        k2 = 5e-4    # Fructose -> Glucose (reverse)
        k3 = 2e-3    # Fructose -> HMF
        k4 = 3e-4    # HMF(aq) -> side products (rehydration)
        k5 = 1e-4    # HMF(aq) -> humins
        k6 = 5e-3    # HMF(org) -> FDCA

        # Mass transfer
        db_LL = 1e-4        # Liquid bubble diameter
        db_GL = 1e-3        # Gas bubble diamer
        db_LS = 1e-3        # Solid particle diameter
        k_aq2org = 2e-2     # MT constant (aq) -> (org)
        a_aq2org = 1/db_LL  # 
        m_aq2org = 1.0      # 
        k_g2org  = 5e-2     # MT constant (g)  -> (org)
        a_g2org  = 1/db_GL  # Hebben wij deze nodig?!?!?!?!?!?!?!?!?!?
        m_g2org  = 1.0
        k_org2p  = 1.0e-2
        a_org2p  = 1/db_LS
        C_O2sat  = 8.0      # O2 saturation in organic [mol/m3]
        
        ### Creating a Matrix
        R1 = k1*V_aq
        R2 = k2*V_aq
        R3 = k3*V_aq
        R4 = k4*V_aq
        R5 = k5*V_aq
        R6 = k6*V_org
        MT_aq2org = k_aq2org * a_aq2org
        MT_g2org  = k_g2org * a_g2org # ??
        MT_org2p  = k_org2p * a_org2p


        self.comps = {0: "Glu(aq)", 1: "Fru(aq)", 2: "HMF(aq)", 3: "HMF(org)", 
                      4: "HMF(p)", 5: "O2(g)", 6: "O2(org)", 7: "O2 (aq)", 8: "O2(p)", 
                      9: "FDCA(p)", 10:"FDCA(org)", 11: "FDCA(aq)", 
                      12:"Formic and Glaric Acids", 13: "Humins"}

        NOfComps = len(self.comps)
        self.balanceMatrix = np.zeros((NOfComps,NOfComps))
        self.balanceMatrix[0, 0:2] = [-(F_aq + R1), R2]                                 # Glucose (aq)
        self.balanceMatrix[1, 0:2] = [R1, -(R2 + R3)]                                   # Fructose (aq)
        self.balanceMatrix[2, 1:4] = [R3, -(R4 + R5 + MT_aq2org), MT_aq2org*m_aq2org]   # HMF (aq)
        self.balanceMatrix[3, 2:4] = [MT_aq2org, -(MT_aq2org*m_aq2org + MT_org2p)]      # HMF (org)
        self.balanceMatrix[4, 3:5] = [MT_org2p, -R6]                                    # HMF (p)
        self.balanceMatrix[5, 0] # O2(g)
        self.balanceMatrix[6, ]



        self.ddt = np.array([
            - F_aq*C_Glu, # Glucose (aq)
            0, # Fructose (aq)
            0,
            0,
            0,
            0
        ])

        print(self.balanceMatrix)
        print(self.ddt)
        """
        self.balanceMatrix = np.matrix(
            [ -k1*V_aq, k2*V_aq, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0], # Glucose (aq)
            [ k1*V_aq, -(k2*V_aq + k3*V_aq), , , , , , , , , , , , ], # Fructose (aq)
            [ , , , , , , , , , , , , , ], # HMF (aq)
            [ , , , , , , , , , , , , , ], # HMF (org)
            [ , , , , , , , , , , , , , ], # Hummins (aq)
            [ , , , , , , , , , , , , , ], # Formic and Glaric Acid (aq)
            [ , , , , , , , , , , , , , ], # O2 (g)
            [ , , , , , , , , , , , , , ], # O2 (org)
            [ , , , , , , , , , , , , , ], # O2 (cat)
            [ , , , , , , , , , , , , , ], # HMF (cat)
            [ , , , , , , , , , , , , , ], # FDCA (cat)
            [ , , , , , , , , , , , , , ], # FDCA (aq)
            [ , , , , , , , , , , , , , ]  # FDCA (org)
            )
        """       

    def components():
        """
        *Creates an array of components.*
        Assumptions:
        * Concentration of Water (l) and Organic Phase molecule (org) is constant
        * Place 
        """
        conc_aq = [
            """
            List containing components in Aqueous Solution
            #  -> Component reacts
            ## -> Component does NOT react
            """
            # Glucose (l)
            # Fructose (l)
            # HMF (l)
            ## Levulinic Acid (l)
            ## Formic Acid (l)
            ## Humins (l)
        ]

        conc_org = [
            # 
        ]

    def matrix():
        """
        Requirements:
        * Matrix needs to have as many columns as the sum of components above
        """
        np.matrix() = [
            [], 
            [],
            [],
            [],
        ]




if __name__ == '__main__':
    pass