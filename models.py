import numpy as np

np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

class ReactorModel:
    def __init__(self):
        # 1. Store parameters as a dictionary for easy sensitivity analysis
        self.params = {
            'F_aq': 2.0, 'F_org': 5.0, 'V_aq': 0.1, 'C_Glu': 500.0,
            'k1': 1e-3, 'k2': 5e-4, 'k3': 2e-3, 'k4': 3e-4, 'k5': 1e-4, 'k6': 5e-3,
            'k_aq2org': 2e2, 'k_org2p': 1.0e2, 'a_aq2org': 1e4, 'a_org2p': 1e3,
            'm_aq2org': 1.0, 'm_org2p': 1.0, 'C_O2sat': 8.0
        }
        self.V_org = 2 * self.params['V_aq']

        self.print_state(self.params, "V_org:", self.V_org)

    def get_balances(self):
        p = self.params
        V_org = self.V_org
        
        # Intermediate Calculations
        R1, R2, R3 = p['k1']*p['V_aq'], p['k2']*p['V_aq'], p['k3']*p['V_aq']
        R4, R5, R6 = p['k4']*p['V_aq'], p['k5']*p['V_aq'], p['k6']*V_org
        MT_aq2org = p['k_aq2org'] * p['a_aq2org']
        MT_org2p  = p['k_org2p'] * p['a_org2p']

        # O2 Pre-calculation
        self.print_state("calculating o2")
        O2_A = np.array([[-MT_org2p, MT_org2p*p['m_org2p']], 
                         [MT_org2p,  -(R6 + MT_org2p*p['m_org2p'])]])
        O2_b = np.array([p['F_org']*p['C_O2sat'], 0])
        O2_results = np.linalg.solve(O2_A, O2_b)
        O2_cat = O2_results[1]
        
        self.print_state("Found O2 on catalyst:", O2_cat)

        R6_prime = R6 * O2_cat

        # Construct Main Matrix A
        self.print_state("Setting up Matrix balances")
        N = 10
        A = np.zeros((N, N))
        A[0, 0:2] = [-(p['F_aq'] + R1), R2]
        A[1, 0:2] = [R1, -(p['F_aq'] + R2 + R3)]
        A[2, 1:4] = [R3, -(R4 + R5 + MT_aq2org), MT_aq2org*p['m_aq2org']]
        A[3, 2:4] = [MT_aq2org, -(MT_aq2org*p['m_aq2org'] + MT_org2p)]
        A[4, 3:5] = [MT_org2p, -(R6_prime + MT_org2p*p['m_org2p'])]
        A[5, 4:6] = [R6_prime, -MT_org2p]
        A[6, 5:8] = [MT_org2p, -(p['F_org'] + MT_aq2org), MT_aq2org]
        A[7, 6:8] = [MT_aq2org, -(p['F_aq'] + MT_aq2org)]
        A[8, 3], A[8, 8] = R4, -p['F_org']
        A[9, 3], A[9, 9] = R5, -p['F_org']
        self.print_state("Constructed matrix:", A)

        # Construct vector b
        b = np.zeros(N)
        b[0] = -p['F_aq'] * p['C_Glu']
        self.print_state("Constructed b:", b)
        
        return A, b

    def solve_steady_state(self):
        A, b = self.get_balances()
        return np.linalg.solve(A, b)

    def ode_system(self, t, y):
        # For dynamics: dy/dt = A*y - b (rearranged from 0 = Ay - b)
        # Note: Since b was defined for the LHS, check signs based on your balance
        A, b = self.get_balances()
        return A @ y - b

    def print_state(self, *args)->None:
        msg_border = "***" + 64*"-" + "***"
        print(msg_border)
        for item in args:
            print(item)
    

if __name__ == '__main__':
    CSTR = ReactorModel()
    print(CSTR.solve_steady_state())