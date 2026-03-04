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
        V_org = 5e-3
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
        kLa_LL  = 2e-2  # (aq) -> (org)
        P_HMF   = 3.0   # partition coefficient HMF: C_org/C_aq [-]
        kLa_GL  = 5e-2  # (g)  -> (org)
        C_O2sat = 8.0   # O2 saturation in organic [mol/m3]
        
        ### Creating a Matrix
        self.balanceMatrix = np.zeros((13,13))
        self.balanceMatrix[0, 0:1] = [-k1*V_aq, k2*V_aq]                # Glucose (aq)
        self.balanceMatrix[1, 0:1] = [k1*V_aq, -(k2*V_aq + k3*V_aq)]    # Fructose (aq)
        self.balanceMatrix[2, 1:2] = [k2*V_aq]



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