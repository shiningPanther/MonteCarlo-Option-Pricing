# MonteCarlo-Option-Pricing


#### Overview
In this script I calculated the price and greeks of a European Down-and-Out barrier option using Monte Carlo simulations. The greeks are obtained by finited difference method.
I compared the results to the analytic calculations of the price and greeks.

In this example I used a strike price of K = 50, barrier B = 45, underlying asset S<sub>0</sub> = 50, volatility of the underlying &sigma; = 0.2, risk free rate r = 0.1 and time to expiry T = 0.5. 

#### Discussion
In order to obtain a good estimate of the exact price and in particular the greeks, one needs to run 10<sup>6</sup> simulations (Note: This number is much smaller for European Vanilla options). Such a large number of simulations takes an unsatisfying long time to run. Solutions to overcome this issue that lead to a faster conversion and that I will implement in the future include:
  
  1. Use of variance reduction techniques, like antithetic variates.
  2. Use of higher order approximations for the greeks.
  3. Use of quasi-random numbers instead of pseudorandom numbers.
  
#### Examples
<img src="/example_price.png" width="500"/>
<img src="/example_delta.png" width="500"/>
<img src="/example_gamma.png" width="500"/>
<img src="/example_theta.png" width="500"/>
<img src="/example_vega.png" width="500"/>
<img src="/example_rho.png" width="500"/>
