# This file was created by Giovanni De Luca and Andrea Montanino  
# Downloading the Litecoin data from Yahoo Finance for a specific period
litecoin = getSymbols("LTC-USD", from= "2018-01-01", to="2025-04-30", src="yahoo", auto.assign = FALSE)

# Define the variable RETURNS
Prlitecoin = ts(litecoin$`LTC-USD.Adjusted`)
returns = diff(log(Prlitecoin)*100)
Acf(returns)
Pacf(returns)
Acf(returns^2)
Pacf(returns^2)
plot(returns, type="l")

# Estimation of the GARCH model
garch_rlitecoin = ugarchspec(mean.model=list(armaOrder=c(0,0)), 
                             variance.model= list(garchOrder=c(1,1), model="sGARCH"), 
                             distribution.model="std")
garchfit_rlitecoin = ugarchfit(garch_rlitecoin, data=returns)
garchfit_rlitecoin
residui_returns=residuals(garchfit_rlitecoin, standardize=T)
LB_m1 = Box.test(residui_returns, lag=1, type="Ljung-Box")
LB_m1 
LB_m7 = Box.test(residui_returns, lag=7, type="Ljung-Box")
LB_m7 
LB_m14 = Box.test(residui_returns, lag=14, type="Ljung-Box")
LB_m14 
LB2_m1 = Box.test(residui_returns^2, lag=1, type="Ljung-Box")
LB2_m1
LB2_m7 = Box.test(residui_returns^2, lag=7, type="Ljung-Box")
LB2_m7
LB2_m14 = Box.test(residui_returns^2, lag=14, type="Ljung-Box")
LB2_m14

residui_core_returns = coredata(residui_returns)
length(residui_core_returns)

# Define the variable VOLUMES
volumes = ts(log(litecoin$`LTC-USD.Volume`), frequency=7)
length(volumes)
volumes_d = diff(volumes)
Acf(volumes_d)
Pacf(volumes_d)
volumes_dd = diff(volumes_d, lag=7)
Acf(volumes_dd)
Pacf(volumes_dd)
Acf(volumes_dd^2)
Pacf(volumes_dd^2)

# Estimation of the SARIMA model
sarima = arima(volumes_dd, order=c(1,0,1), seasona=list(order=c(0,0,1)))
sarima
coeff <- coef(sarima)
std <- sqrt(diag(vcov(sarima))) 
t_values <- coeff/std
p_values <- 2 * (1 - pnorm(abs(t_values)))
results <- data.frame(Coefficients = coeff, Std_Errors = std, T_values = t_values, P_values = p_values)
print(round(results,3))
resid_sarima = sarima$residuals
Acf(resid_sarima)
Pacf(resid_sarima)
Acf(resid_sarima^2)
Pacf(resid_sarima^2)

LB_acf1 = Box.test(resid_sarima, lag = 1, type="Ljung-Box")
LB_acf1
LB_acf7 = Box.test(resid_sarima, lag = 7, type="Ljung-Box")
LB_acf7
LB_acf14 = Box.test(resid_sarima, lag = 14, type="Ljung-Box")
LB_acf14

# Estimation of the GARCH model
garch_volumes = ugarchspec(mean.model=list(armaOrder=c(0,0)), 
                          variance.model= list(garchOrder=c(1,2), model="eGARCH"),
                          distribution.model="sstd")
garchfit_volumes = ugarchfit(garch_volumes, data=resid_sarima)
garchfit_volumes
#plot(garchfit_volumes)
residui_volumes = residuals(garchfit_volumes, standardize=T)
LB_acf1 = Box.test(residui_volumes, lag=1, type="Ljung-Box")
LB_acf1 
LB_acf7 = Box.test(residui_volumes, lag=7, type="Ljung-Box")
LB_acf7 
LB_acf14 = Box.test(residui_volumes, lag=14, type="Ljung-Box")
LB_acf14
LB2_acf1 = Box.test(residui_volumes^2, lag=1, type="Ljung-Box")
LB2_acf1
LB2_acf7 = Box.test(residui_volumes^2, lag=7, type="Ljung-Box")
LB2_acf7
LB2_acf14 = Box.test(residui_volumes^2, lag=14, type="Ljung-Box")
LB2_acf14
residui_volumes_core = coredata(residui_volumes)
residui_volumes_true= as.numeric(residui_volumes_core)
length(residui_volumes_true)

residui_returns_true = residui_core_returns[8:2676]
length(residui_returns_true)


# Computation of the empirical distrbution functions
pobs_returns = pobs(residui_returns_true)
pobs_volumes = pobs(residui_volumes_true)
plot(pobs_returns, pobs_volumes)

# Application of the copula model to the residuals
data_pobs = cbind(pobs_returns, pobs_volumes)
j=joeCopula()
j90 = rotCopula(joeCopula(), flip=c(T,F))

# Definition of the mixture
Mixt_cop_pobs <- list(j, j90)
Mix_cop_pobs <- mixCopula(Mixt_cop_pobs, w=c())
fit_Mix_cop_pobs <- fitCopula(Mix_cop_pobs, data_pobs, method="ml")
print(round(fit_Mix_cop_pobs@estimate,3))
standard_error <- sqrt(diag(fit_Mix_cop_pobs@var.est))
print(round(standard_error,3))

j_est = joeCopula(fit_Mix_cop_pobs@estimate[1])
j90_est = rotCopula(joeCopula(fit_Mix_cop_pobs@estimate[2]), flip=c(T,F))

# Definiton of the estimated copula
Mixture_cop_est <- list(j_est, j90_est)
Mix_cop_est <- mixCopula(Mixture_cop_est, w=c(fit_Mix_cop_pobs@estimate[3], fit_Mix_cop_pobs@estimate[4]))

# Goodness-of-fit measure 
testGOF_mixture_not_parallel <- gofCopula(Mix_cop_est, data_pobs, N = 1000, method = "Sn",
                                          estim.method =  "mpl", simulation = "mult", test.method = "family")
testGOF_mixture_not_parallel

# Creation of the contour plot
wireframe2(Mix_cop_est, dCopula, n = 30, drawlabels = FALSE,
           scales = list(arrows = FALSE, 
                         y = list(rot = 90, distance = 2)),  # Aumenta la distanza dell'asse Y
           col = "blue", 
           xlab = "Returns", 
           ylab = "Volumes", 
           zlab = "",
           screen = list(z = 15, x = -55)) 

# Tail dependence coefficients
lambda_uu = fit_Mix_cop_pobs@estimate[3]*BiCopPar2TailDep(6,fit_Mix_cop_pobs@estimate[1])$upper
round(lambda_uu,3)
lambda_lu = fit_Mix_cop_pobs@estimate[4]*BiCopPar2TailDep(6,fit_Mix_cop_pobs@estimate[2])$upper
round(lambda_lu,3)

