'''
------------------------------------------------------------------------
Solves the dynamic programming problem of the firm with:
- Quadratic adjustment costs
- TFP/Profit shocks
- No taxes
- No external finance
- Partial equilibrium

This py-file calls the following other file(s):

This py-file creates the following other file(s):
    (make sure that an OUTPUT folder exists)
            graphs/
------------------------------------------------------------------------
'''

# Import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import numba
import ar1_approx as ar1

# Create directory if images directory does not already exist
cur_path = os.path.split(os.path.abspath(__file__))[0]
output_fldr = 'images'
output_dir = os.path.join(cur_path, output_fldr)
if not os.access(output_dir, os.F_OK):
    os.makedirs(output_dir)

'''
------------------------------------------------------------------------
Define functions
------------------------------------------------------------------------
'''


@numba.jit
def VFI_loop(EV, e, betafirm, Pi, sizez, sizek, Vmat):
    '''
    ------------------------------------------------------------------------
    This function loops over the state and control variables, operating on the
    value function to update with the last iteration's value function
    ------------------------------------------------------------------------
    INPUTS:
    EV       = (sizez, sizek) matrix, expected value function (expectations
               over z')
    e        = (sizek, sizek) matrix, cash flow values for each possible
               combination of capital stock today (state) and choice of capital
               stock tomorrow (control)
    betafirm = scalar in [0, 1], the discount factor of the firm
    Pi       = (sizez, sizez) matrix, transition probabilities between points
               in the productivity state space
    sizez    = integer, number of grid points for firm productivity shocks
               state space
    sizek    = integer, the number of grid points in capital space
    Vmat     = (sizek, sizek) matrix, matrix with values of firm at each
               combination of state (k) and control (k')

    OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None

    OBJECTS CREATED WITHIN FUNCTION: None

    FILES CREATED BY THIS FUNCTION: None

    RETURNS: Vmat
    ------------------------------------------------------------------------
    '''
    for i in range(sizez):  # loop over z
        for j in range(sizek):  # loop over k
            for m in range(sizek):  # loop over k'
                Vmat[i, j, m] = e[i, j, m] + betafirm * EV[i, m]

    return Vmat


'''
------------------------------------------------------------------------
Specify Parameters
------------------------------------------------------------------------
beta      = scalar in (0, 1), rate of time preference
alpha_k   = scalar in [0, 1], exponent on capital in firm production function
alpha_l   = scalar in [0, 1], exponent on labor in firm production function
delta     = scalar in [0, 1], depreciation rate on capital
psi       = scalar, coefficient in quadratic adjustment costs for capital
w         = scalar, exogenous wage rate
r         = scalar, risk free interest rate, in eqm, r = (1 / beta) - 1
beta_firm = scalar in [0, 1], the discount factor of the firm
sigma_eps = scalar > 0, standard deviation of profit/productivity shocks to
            AR(1) process for firm productivity
mu        = scalar, unconditional mean of productivity shocks
rho       = scalar in [0, 1], persistence of productivity shocks
sizez     = integer, number of grid points for firm productivity shocks state
            space
------------------------------------------------------------------------
'''
beta = 0.96
alpha_k = 0.29715
alpha_l = 0.65
delta = 0.154
psi = 1.08
w = 1.3
r = ((1 / beta) - 1)
betafirm = (1 / (1 + r))
sigma_eps = 0.213
mu = 0
rho = 0.7605
sizez = 9


'''
-------------------------------------------------------------------------
Discretizing state space for productivity shocks
-------------------------------------------------------------------------
sigma_z   = scalar, standard deviation of ln(z)
num_sigma = scalar, number of standard deviations around mean to include in
            grid space for z
step      = scalar, distance between grid points in the productivity state
            space
Pi        = (sizez, sizez) matrix, transition probabilities between points in
            the productivity state space
z         = (sizez,) vector, grid points in the productivity state space
-------------------------------------------------------------------------
'''
# We will use the Rouwenhorst (1995) method to approximate a continuous
# distribution of shocks to the AR1 process with a Markov process.
sigma_z = sigma_eps / ((1 - rho ** 2) ** (1 / 2))
num_sigma = 3
step = (num_sigma * sigma_z) / (sizez / 2)
Pi, z = ar1.rouwen(rho, mu, step, sizez)
Pi = np.transpose(Pi)  # make so rows are where start, columns where go
z = np.exp(z)  # because the AR(1) process was for the log of productivity

'''
-------------------------------------------------------------------------
Discretizing state space for capital
-------------------------------------------------------------------------
dens   = integer, density of the grid: number of grid points between k and
         (1 - delta) * k
kstar  = scalar, capital stock choose w/o adjustment costs and mean
         productivity shock
kbar   = scalar, maximum capital stock the firm would ever accumulate
ub_k   = scalar, upper bound of capital stock space
lb_k   = scalar, lower bound of capital stock space
krat   = scalar, difference between upper and lower bound in log points
numb   = integer, the number of steps between the upper and lower bounds for
         the capital stock. The number of grid points is dens*numb.
K      = (sizek,) vector, grid points in the capital state space, from high
         to low
kvec  = (sizek,) vector, capital grid points
sizek = integer, the number of grid points in capital space
-------------------------------------------------------------------------
'''
dens = 5
# put in bounds here for the capital stock space
kstar = ((((1 / betafirm - 1 + delta) * ((w / alpha_l) **
                                         (alpha_l / (1 - alpha_l)))) /
         (alpha_k * (z[(sizez - 1) // 2] ** (1 / (1 - alpha_l))))) **
         ((1 - alpha_l) / (alpha_k + alpha_l - 1)))
kbar = kstar * 500
lb_k = 0.001
ub_k = kbar
krat = np.log(lb_k / ub_k)
numb = np.ceil(krat / np.log(1 - delta))
K = np.empty(int(numb * dens))
for j in range(int(numb * dens)):
    K[j] = ub_k * (1 - delta) ** (j / dens)
kvec = K[::-1]
sizek = kvec.shape[0]

'''
-------------------------------------------------------------------------
Generating possible cash flow levels
-------------------------------------------------------------------------
op = (sizez, sizek) matrix, operating profits for each point in capital stock
     and productivity shock grid spaces
e  = (sizez, sizek, sizek) array, cash flow values for each possible
     combination of current productivity shock (state), capital stock today
     (state), and choice of capital stock tomorrow (control)
-------------------------------------------------------------------------
'''
op = np.empty((sizez, sizek))
for i in range(sizez):
    for j in range(sizek):
        op[i, j] = ((1 - alpha_l) * ((alpha_l / w) **
                                     (alpha_l / (1 - alpha_l))) *
                    ((z[i] * kvec[j] ** alpha_k) **
                     (1 / (1 - alpha_l))))

e = np.empty((sizez, sizek, sizek))
for i in range(sizez):
    for j in range(sizek):
        for m in range(sizek):
            e[i, j, m] = (op[i, j] - kvec[m] + ((1 - delta) * kvec[j]) -
                          ((psi / 2) *
                           (((kvec[m] - ((1 - delta) * kvec[j])) ** 2)
                            / kvec[j])))


'''
------------------------------------------------------------------------
Value Function Iteration
------------------------------------------------------------------------
VFtol     = scalar, tolerance required for value function to converge
VFdist    = scalar, distance between last two value functions
VFmaxiter = integer, maximum number of iterations for value function
VFiter    = integer, current iteration number
Vmat      = (sizez, sizek, sizek) array, array with values of firm at each
            combination of state (z, k) and control (k')
Vstore    = (sizez, sizek, VFmaxiter) array, value function at each iteration
            of VFI
V & TV    = (sizez, sizek) matrix, store the value function at each iteration
            (V being the most current value and TV the one prior)
EV        = (sizez, sizek) matrix, expected value function (expectations over
            z')
PF        = (sizez, sizek) matrix, indicies of choices (k') for all states
            (z, k)
VF        = (sizez, sizek) matrix, matrix of value functions for each possible
            value of the state variables (k)
------------------------------------------------------------------------
'''
VFtol = 1e-6
VFdist = 7.0
VFmaxiter = 3000
V = np.empty((sizez, sizek))  # initial guess at value function
Vmat = np.empty((sizez, sizek, sizek))  # initialize Vmat matrix
Vstore = np.empty((sizez, sizek, VFmaxiter))  # initialize Vstore array
VFiter = 1
start_time = time.clock()
while VFdist > VFtol and VFiter < VFmaxiter:
    TV = V
    EV = np.dot(Pi, V)  # expected VF (expectation over z')
    Vmat = VFI_loop(EV, e, betafirm, Pi, sizez, sizek, Vmat)

    Vstore[:, :, VFiter] = V.reshape(sizez, sizek)  # store value function at
    # each iteration for graphing later
    V = Vmat.max(axis=2)  # apply max operator to Vmat (to get V(k))
    PF = np.argmax(Vmat, axis=2)
    VFdist = (np.absolute(V - TV)).max()  # check distance between value
    # function for this iteration and value function from past iteration
    VFiter += 1

VFI_time = time.clock() - start_time
if VFiter < VFmaxiter:
    print('Value function converged after this many iterations:', VFiter)
else:
    print('Value function did not converge')
print('VFI took ', VFI_time, ' seconds to solve')


VF = V  # solution to the functional equation

'''
------------------------------------------------------------------------
Find optimal capital and investment policy functions
------------------------------------------------------------------------
optK = (sizez, sizek) vector, optimal choice of k' for each (z, k)
optI = (sizez, sizek) vector, optimal choice of investment for each (z, k)
------------------------------------------------------------------------
'''
optK = kvec[PF]
optI = optK - (1 - delta) * kvec


'''
------------------------------------------------------------------------
Graph Value Function and Decision Rules
------------------------------------------------------------------------
'''

# Plot value function
# plt.figure()
fig, ax = plt.subplots()
ax.plot(kvec, VF[0, :], 'k--', label='z = ' + str(z[0]))
ax.plot(kvec, VF[(sizez - 1) // 2, :], 'k:', label='z = ' + str(z[(sizez - 1)
                                                                  // 2]))
ax.plot(kvec, VF[-1, :], 'k', label='z = ' + str(z[-1]))
# Now add the legend with some customizations.
legend = ax.legend(loc='lower right', shadow=True)
# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')
# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')
for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
plt.xlabel('Size of Capital Stock')
plt.ylabel('Value Function')
plt.title('Value Function - stochastic firm w/ adjustment costs')
output_path = os.path.join(output_dir, 'V_firm4')
plt.savefig(output_path, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()


# Plot value function at several iterations
# plt.figure()
fig, ax = plt.subplots()
ax.plot(kvec, Vstore[(sizez - 1) // 2, :, 0], 'k--', label='1st iter')
ax.plot(kvec, Vstore[(sizez - 1) // 2, :, 1], 'k:', label='2nd iter')
ax.plot(kvec, Vstore[(sizez - 1) // 2, :, 2], 'k:', label='3rd iter')
ax.plot(kvec, Vstore[(sizez - 1) // 2, :, 3], 'k:', label='4th iter')
ax.plot(kvec, Vstore[(sizez - 1) // 2, :, VFiter-1], 'k', label='Last iter')
# Now add the legend with some customizations.
legend = ax.legend(loc='lower right', shadow=True)
# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')
# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')
for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
plt.xlabel('Size of Capital Stock')
plt.ylabel('Value Function, z = ' + str(z[(sizez - 1) // 2]))
plt.title('Value Function - stochastic firm w/ adjustment costs')
output_path = os.path.join(output_dir, 'V_firm4_iter')
plt.savefig(output_path, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()


# Plot optimal capital stock rule as a function of firm size
# plt.figure()
fig, ax = plt.subplots()
ax.plot(kvec, optK[0, :], 'k--', label='z = ' + str(z[0]))
ax.plot(kvec, optK[(sizez - 1) // 2, :], 'k:', label='z = ' + str(z[(sizez - 1)
                                                                    // 2]))
ax.plot(kvec, optK[-1, :], 'k', label='z = ' + str(z[-1]))
ax.plot(kvec, kvec, 'k:', label='45 degree line')
# Now add the legend with some customizations.
legend = ax.legend(loc='upper left', shadow=True)
# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')
# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')
for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
plt.xlabel('Size of Capital Stock')
plt.ylabel('Optimal Choice of Capital Next Period')
plt.title('Policy Function, Next Period Capital - stochastic firm w/ ' +
          'adjustment costs')
output_path = os.path.join(output_dir, 'Kprime_firm4')
plt.savefig(output_path, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()


# Plot investment rule as a function of firm size
# plt.figure()
fig, ax = plt.subplots()
ax.plot(kvec, optI[(sizez - 1) // 2, :], 'k--', label='Investment')
ax.plot(kvec, kvec, 'k:', label='45 degree line')
# Now add the legend with some customizations.
legend = ax.legend(loc='upper left', shadow=True)
# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')
# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')
for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
plt.xlabel('Size of Capital Stock')
plt.ylabel('Optimal Investment')
plt.title('Policy Function, Investment - stochastic firm w/ adjustment ' +
          'costs')
output_path = os.path.join(output_dir, 'invest_firm4')
plt.savefig(output_path, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()
