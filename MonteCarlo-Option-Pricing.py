
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st


# Parameters of the option
K = 50 # Strike price
S0 = 50 # Current underlying price
B = 45 # Barrier level
sigma = 0.2 # Volatility of the underlying
r = .1 # Risk free rate
T = .5 # Half a year to expiry


def simulateUnderlying(W, t, r=r, sigma=sigma,  S0=S0):
	'''
	Calculates the price of the underlying using Geometric Brownian motion and risk free rate r
	'''
	X = (r-0.5*sigma**2)*t + sigma*W
	return S0*np.exp(X)

def calculateOptionPrice(S):
	'''
	Returns the option price taking an array of the evolution of the underlying as input
	'''
	# Value of the option: Discounted expected return of the option
	price = (S[-1]-K)*np.exp(-r*T)
	# If we do not hit the barrier set option price to 0
	if S[-1] < K:
		price = 0
	# If we are knocked out by the lower barrier the option becomes worthless
	if min(S) <= B:
		price = 0
	return(price)

def addCumulativeAverage(array, element):
	'''
	Appends the cumulative average of the option price for each additional simulation to the input array.
	This array is only needed to plot the evolution of the option as a function of the number of simulations.
	'''
	length=len(array)
	if length == 0:
		array.append(element)
	else:
		array.append((array[length-1]*length+element)/(length+1))


def calculatePrices_MC(K=K, S0=S0, B=B, sigma=sigma, r=r, T=T):
	'''
	Calculates the price of the option and the deltas using MC simulations and finite difference method.
	'''

	# Result arrays
	C = [] # Price of the option, as simulated in each run
	C_avg = [] # Price of the option, averaged after each run
	delta = []
	delta_avg = []
	gamma = []
	gamma_avg = []
	vega = []
	vega_avg = []
	theta = []
	theta_avg = []
	rho = []
	rho_avg = []
	
	# Simulation parameters
	dt = 0.0001 # Time interval
	N = round(T/dt) # Number of steps in the Brownian motion
	t = np.linspace(0, T, N)
	Nsim = 1000000 # Number of simulations

	# Simulation parameters for calculating the greeks via FDM
	dS0 = S0*0.01
	dsigma = sigma*0.01
	dr = r*0.01

	for i in range(Nsim):
		# Standard Brownian motion 
		W = np.random.standard_normal(size = N) 
		W = np.cumsum(W)*np.sqrt(dt)

		# The evolution of the underlying, as needed for calculating the price of the option (using risk-free rate r)
		S_std = simulateUnderlying(W, t)
		optionPrice_std = calculateOptionPrice(S_std)
		C.append(optionPrice_std)
		addCumulativeAverage(C_avg, optionPrice_std)

		# Calculation of delta
		S_delta = simulateUnderlying(W, t, S0 = S0+dS0)
		optionPrice_delta = calculateOptionPrice(S_delta)
		delta_sim = (optionPrice_delta - optionPrice_std) / dS0
		delta.append(delta_sim)
		addCumulativeAverage(delta_avg, delta_sim)

		# Calculation of gamma
		S_gamma = simulateUnderlying(W, t, S0 = S0-dS0)
		optionPrice_gamma = calculateOptionPrice(S_gamma)
		gamma_sim = (optionPrice_delta - 2*optionPrice_std + optionPrice_gamma) / dS0**2
		gamma.append(gamma_sim)
		addCumulativeAverage(gamma_avg, gamma_sim)

		# Calculation of vega
		S_vega = simulateUnderlying(W, t, sigma = sigma+dsigma)
		optionPrice_vega = calculateOptionPrice(S_vega)
		vega_sim = (optionPrice_vega - optionPrice_std) / dsigma
		vega.append(vega_sim)
		addCumulativeAverage(vega_avg, vega_sim)

		# Calculation of theta
		Ndt = 10 # Number of steps of dt being used for the finite difference of T
		W_theta = W[:-Ndt]
		t_theta = np.linspace(0,T-Ndt*dt,N-Ndt)
		S_theta = simulateUnderlying(W_theta, t = t_theta)
		optionPrice_theta = calculateOptionPrice(S_theta)*np.exp(r*Ndt*dt)
		theta_sim = (optionPrice_theta - optionPrice_std) / (Ndt*dt)
		theta.append(theta_sim)
		addCumulativeAverage(theta_avg, theta_sim)

		# Calculation of rho
		S_rho = simulateUnderlying(W, t, r = r+dr)
		optionPrice_rho = calculateOptionPrice(S_rho)*np.exp(-dr*T)
		rho_sim = (optionPrice_rho - optionPrice_std) / dr
		rho.append(rho_sim)
		addCumulativeAverage(rho_avg, rho_sim)

	return C, C_avg, delta, delta_avg, gamma, gamma_avg, vega, vega_avg, theta, theta_avg, rho, rho_avg


def calculate_dj(S0, r, sigma, T):
	d1 = (np.log(S0/K) + (r + 0.5*sigma**2)*T)/(sigma*np.sqrt(T))
	d2 = d1 - sigma*np.sqrt(T)
	return d1, d2

def calculateVanillaPrice(S0, r, sigma, T):
	d1, d2 = calculate_dj(S0, r, sigma, T)
	return S0*st.norm.cdf(d1) - K*np.exp(-r*T)*st.norm.cdf(d2)

def calculateBarrierPrice(S0=S0, r=r, sigma=sigma, T=T):
	alpha = 0.5*(1-r/(0.5*sigma**2))
	return calculateVanillaPrice(S0, r, sigma, T) - (S0/B)**(2*alpha)*calculateVanillaPrice(B**2/S0, r, sigma, T) 


def calculatePrices_Analytic(K=K, S0=S0, B=B, sigma=sigma, r=r, T=T):
	'''
	Calculates the price of the option analytically.
	The deltas are calculated using the finite difference method.
	'''

	# W can use smaller increments for the parameters here, in order to achieve a more precise value of the greeks
	dS0_analytic = S0*0.0000001
	dsigma_analytic = sigma*0.0000001
	dr_analytic = r*0.0000001
	dt_analytic = T*0.0000001

	C = calculateBarrierPrice(S0, r, sigma, T)
	delta = (calculateBarrierPrice(S0+dS0_analytic, r, sigma, T) - optionPrice_analytic) / dS0_analytic
	gamma = (calculateBarrierPrice(S0+dS0_analytic, r, sigma, T) - 2*optionPrice_analytic + calculateBarrierPrice(S0-dS0_analytic, r, sigma, T)) / dS0_analytic**2
	vega = (calculateBarrierPrice(S0, r, sigma+dsigma_analytic, T) - optionPrice_analytic) / dsigma_analytic
	theta = (calculateBarrierPrice(S0, r, sigma, T-dt_analytic) - optionPrice_analytic) / dt_analytic
	rho = (calculateBarrierPrice(S0, r+dr_analytic, sigma, T) - optionPrice_analytic) / dr_analytic

	return C, delta, gamma, vega, theta, rho





C, C_avg, delta, delta_avg, gamma, gamma_avg, vega, vega_avg, theta, theta_avg, rho, rho_avg = calculatePrices_MC()

print('Option price MC: {}'.format(np.average(C)))
print('Delta MC: {}'.format(np.average(delta)))
print('Gamma MC: {}'.format(np.average(gamma)))
print('Vega MC: {}'.format(np.average(vega)))
print('Theta MC: {}'.format(np.average(theta)))
print('Rho MC: {}'.format(np.average(rho)))

C_analytic, delta_analytic, gamma_analytic, vega_analytic, theta_analytic, rho_analytic = calculatePrices_Analytic()

print('Option price analytic: {}'.format(C_analytic))
print('Delta analytic: {}'.format(delta_analytic))
print('Gamma analytic: {}'.format(gamma_analytic))
print('Vega analytic: {}'.format(vega_analytic))
print('Theta analytic: {}'.format(theta_analytic))
print('Rho analytic: {}'.format(rho_analytic))


#Plot the errors vs the of number of simulations
Nlower = 100000 # lower range of simulations from where to start the plots
Nsim_array = range(1, Nsim+1)
fig1, ax1 = plt.subplots()
ax1.plot(Nsim_array[Nlower-1:], (C_avg[Nlower-1:]/optionPrice_analytic - 1)*100)
ax1.set_title('Percentage error of the option price')
fig2, ax2 = plt.subplots()
ax2.plot(Nsim_array[Nlower-1:], (delta_avg[Nlower-1:]/delta_analytic - 1)*100)
ax2.set_title('Percentage error of delta')
fig3, ax3 = plt.subplots()
ax3.plot(Nsim_array[Nlower-1:], (gamma_avg[Nlower-1:]/gamma_analytic - 1)*100)
ax3.set_title('Percentage error of gamma')
fig4, ax4 = plt.subplots()
ax4.plot(Nsim_array[Nlower-1:], (vega_avg[Nlower-1:]/vega_analytic - 1)*100)
ax4.set_title('Percentage error of vega')
fig5, ax5 = plt.subplots()
ax5.plot(Nsim_array[Nlower-1:], (theta_avg[Nlower-1:]/theta_analytic - 1)*100)
ax5.set_title('Percentage error of theta')
fig6, ax6 = plt.subplots()
ax6.plot(Nsim_array[Nlower-1:], (rho_avg[Nlower-1:]/rho_analytic - 1)*100)
ax6.set_title('Percentage error of rho')

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
	ax.set_xlabel('Number of simulations')
	ax.set_ylabel('Error in percent')

plt.show()


