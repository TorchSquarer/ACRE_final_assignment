import numpy as np
from scipy.integrate import solve_ivp # New import for ODEs
import matplotlib.pyplot as plt

np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

class ReactorModel:
    def __init__(self):
        self.params = {
            'F_aq': 2.0, 'F_org': 5.0, 'V_aq': 0.1, 'C_Glu': 500.0,
            'k1': 1e-3, 'k2': 5e-4, 'k3': 2e-3, 'k4': 3e-4, 'k5': 1e-4, 'k6': 5e-3,
            'k_aq2org': 2e2, 'k_org2p': 1.0e2, 'a_aq2org': 1e4, 'a_org2p': 1e3,
            'm_aq2org': 1.0, 'm_org2p': 1.0, 'C_O2sat': 8.0
        }
        self.V_org = 2 * self.params['V_aq']

    def get_balances(self):
        p = self.params
        V_org = self.V_org
        
        R1, R2, R3 = p['k1']*p['V_aq'], p['k2']*p['V_aq'], p['k3']*p['V_aq']
        R4, R5, R6 = p['k4']*p['V_aq'], p['k5']*p['V_aq'], p['k6']*V_org
        MT_aq2org = p['k_aq2org'] * p['a_aq2org']
        MT_org2p  = p['k_org2p'] * p['a_org2p']

        O2_A = np.array([[-MT_org2p, MT_org2p*p['m_org2p']], 
                         [MT_org2p,  -(R6 + MT_org2p*p['m_org2p'])]])
        O2_b = np.array([p['F_org']*p['C_O2sat'], 0])
        O2_results = np.linalg.solve(O2_A, O2_b)
        O2_cat = O2_results[1]
        
        R6_prime = R6 * O2_cat

        N = 10
        A = np.zeros((N, N))
        A[0, 0:2] = [-(p['F_aq'] + R1), R2]
        A[1, 0:2] = [R1, -(p['F_aq'] + R2 + R3)]
        A[2, 1:4] = [R3, -(R4 + R5 + MT_aq2org), MT_aq2org*p['m_aq2org']]
        A[3, 2:4] = [MT_aq2org, -(MT_aq2org*p['m_aq2org'] + MT_org2p)]
        # Corrected index 4 row (HMF_p)
        A[4, 3:5] = [MT_org2p, -(R6_prime + MT_org2p*p['m_org2p'])]
        A[5, 4:6] = [R6_prime, -MT_org2p]
        A[6, 5:8] = [MT_org2p, -(p['F_org'] + MT_aq2org), MT_aq2org]
        A[7, 6:8] = [MT_aq2org, -(p['F_aq'] + MT_aq2org)]
        A[8, 3], A[8, 8] = R4, -p['F_org']
        A[9, 3], A[9, 9] = R5, -p['F_org']

        b = np.zeros(N)
        b[0] = -p['F_aq'] * p['C_Glu']
        
        return A, b

    def solve_steady_state(self):
        A, b = self.get_balances()
        return np.linalg.solve(A, b)

    def ode_system(self, t, y):
        """Standard form for solve_ivp: dy/dt = f(t, y)"""
        A, b = self.get_balances()
        return A @ y - b

    def solve_dynamics(self, t_end=3600, n_points=1000):
        """Solves the ODE system from t=0 to t_end."""
        # Initial concentrations: all zero
        y0 = np.zeros(10) 
        t_span = (0, t_end)
        t_eval = np.linspace(0, t_end, n_points)

        # Solving the IVP
        sol = solve_ivp(self.ode_system, t_span, y0, t_eval=t_eval, method='BDF')
        
        return sol

    def print_state(self, *args):
        msg_border = "***" + 64*"-" + "***"
        print(msg_border)
        for item in args:
            print(item)

if __name__ == '__main__':
    CSTR = ReactorModel()
    
    # Solve Dynamics
    solution = CSTR.solve_dynamics(t_end=500)
    
    # Plotting example
    plt.figure(figsize=(10, 6))
    plt.plot(solution.t, solution.y[0], label="Glucose (aq)")
    plt.plot(solution.t, solution.y[2], label="HMF (aq)")
    plt.plot(solution.t, solution.y[7], label="FDCA (aq)")
    plt.xlabel("Time [s]")
    plt.ylabel("Concentration [mol/m3]")
    plt.legend()
    plt.grid(True)
    plt.show()