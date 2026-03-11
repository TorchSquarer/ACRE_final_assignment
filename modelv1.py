import numpy as np
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

class model():
    def __init__(self)-> None:
        # --------------------------------------------------------
        # Parameters
        # --------------------------------------------------------
        # Flow rates [m3/s]
        F_aq  = 2
        F_org = 5
        F_R   = F_aq + F_org

        # Phase volumes [m3]
        V_aq  = 0.1
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
        db_LS = 1e-3        # Solid particle diameter
        k_aq2org = 2e2     # MT constant (aq) -> (org)
        a_aq2org = 1/db_LL  # 
        m_aq2org = 1.0      # 
        m_org2p  = 1.0      # 
        k_org2p  = 1.0e2
        a_org2p  = 1/db_LS 
        C_O2sat  = 8.0      # O2 saturation in organic [mol/m3]
        
        # ----------------------------------
        # Matrix setup Steady State)
        # ----------------------------------

        ### Creating a Matrix elemental shorthands
        R1 = k1*V_aq
        R2 = k2*V_aq
        R3 = k3*V_aq
        R4 = k4*V_aq
        R5 = k5*V_aq
        R6 = k6*V_org
        MT_aq2org = k_aq2org * a_aq2org
        MT_org2p  = k_org2p * a_org2p

        # Calculation of O2 beforehand, assuming effective reaction Glucose <=> Fructose is slow
        self.print_state("Setup of oxygen balances")
        self.O2_Balances = np.array([
            [-MT_org2p, MT_org2p*m_org2p],            # O2 (org)
            [MT_org2p,  -(R6 + MT_org2p*m_org2p)]     # O2 (p)
        ])

        self.O2_b = np.array([F_org*C_O2sat, 0])
        self.print_state("O2 Balances:", self.O2_Balances, "O2 b:", self.O2_b, "Solving matrix")

        self.O2_Results = np.linalg.solve(self.O2_Balances, self.O2_b)
        self.O2_cat = self.O2_Results[1]
        self.print_state("O2 concentration on catalyst:", self.O2_cat)

        R6_prime = R6*self.O2_cat
        ### Listing components/variables
        #self.comps = {0: "Glu(aq)", 1: "Fru(aq)", 2: "HMF(aq)", 3: "HMF(org)", 
        #              4: "HMF(p)", 5: "O2(org)", 6: "O2(p)", 7: "FDCA(p)", 8:"FDCA(org)", 
        #              9:"Formic and Glaric Acids", 10: "Humins"}
        self.comps = {0: "Glu(aq)", 1: "Fru(aq)", 2: "HMF(aq)", 3: "HMF(org)", 
                      4: "HMF(p)", 5: "FDCA(p)", 6:"FDCA(org)", 7: "FDCA(aq)",
                      8:"Formic and Glaric Acids", 9: "Humins"}



        N_Of_Comps = len(self.comps)
        self.ReactantBalances = np.zeros((N_Of_Comps,N_Of_Comps))
        self.ReactantBalances[0, 0:2] = [-(F_aq + R1), R2]                    # Glucose (aq)
        self.ReactantBalances[1, 0:2] = [R1,         -(F_aq+ R2 + R3)]        # Fructose (aq)
        self.ReactantBalances[2, 1:4] = [R3,         -(R4 + R5 + MT_aq2org), 
                                         MT_aq2org*m_aq2org]                  # HMF (aq)
        self.ReactantBalances[3, 2:4] = [MT_aq2org,  -(MT_aq2org*m_aq2org + MT_org2p)]  # HMF (org)
        self.ReactantBalances[4, 3:5] = [MT_org2p ,  -(R6_prime + MT_org2p*m_org2p)]    # HMF (p)
        self.ReactantBalances[5, 4:6] = [R6_prime,    -MT_org2p]              # FDCA (p) (no m_org2p and TM one direction)
        self.ReactantBalances[6, 5:8] = [MT_org2p, -(F_org + MT_aq2org), MT_aq2org]      # FDCA (org) (no m_org2p or m_aq2org)
        self.ReactantBalances[7, 6:8] = [MT_aq2org, -(F_aq + MT_aq2org)]      # FDCA (aq) (no m_org2p)
        self.ReactantBalances[8, 3]   = R4                                    # Acids
        self.ReactantBalances[8, 8]   = F_org
        self.ReactantBalances[9, 3]   = R5                                    # Humins
        self.ReactantBalances[9, 9]   = F_org
        
        ### Matrix check
        self.print_state("Reactant balances initialized as:", self.ReactantBalances)

        self.ddt = np.zeros((N_Of_Comps, 1))
        self.ddt[0] = - F_aq*C_Glu
        self.print_state("Printing b:", self.ddt)

    def print_state(self, *args)->None:
        msg_border = "***" + 64*"-" + "***"
        print(msg_border)
        for item in args:
            print(item)

        #print(msg_border)
    
    def solve_balances(self)->None:
        self.conc_results = np.linalg.solve(self.ReactantBalances, self.ddt)
        self.print_state("Found concentrations:", self.conc_results)


if __name__ == '__main__':
    CSTR = model()
    CSTR.solve_balances()