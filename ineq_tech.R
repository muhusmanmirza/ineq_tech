#Library----
library(easypackages)
libraries('phaseR', 'deSolve', 'rootSolve', 'ggplot2')
load("C:/MinGW/msys/1.0/home/mirza009/ineq_tech/ineq_tech.RData")

#Model abiotic----
n       <- 100
alpha   <- 0.2
c       <- 1
delta   <- 0.2
gamma   <- 0.2
k       <- 1
lambda1 <- 20
m       <- 1000
mu      <- 0.1
phi     <- 1
rho     <- 0.5
#s       <- abs(rnorm(100, 0.15, 0.01))
s       <- 0.25
theta1  <- 50
#iv      <- c(a = runif(100, 0, 100), X = 500)
iv      <- c(a = 1, X = 500)
t       <- seq(0, 1000, 1)

X_f <- function(X, c, m, n) {c*(m - X) - n*X}
plot(0:20, X_f(0:20, 2, 15, 0), type = 'l'); abline(h = 0, col = 'grey'); abline(v = 0, col = 'grey')

model_abio_2d <- function(t, y, parameters) {
  a <- y[1]
  X <- y[2]
  ext <- k*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)*phi*a^(alpha + delta)*X^gamma 
  da <- a^(alpha+alpha+delta)*s*k^2*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)^2*phi*X^gamma - a*mu
  dX <- c*(m - X) - n*ext
  list(c(da, dX))
}

out_2d <- ode(iv, t, model_abio_2d, NULL)
matplot(out_2d[,1], out_2d[,-1], type = 'l')

ff <- flowField(model_abio_2d, x.lim = c(0, 1000), y.lim = c(0, 1000), points = 20, 
                system = 'two.dim', add = F, xlab = 'a', ylab = 'X')
null <- nullclines(model_abio_2d, x.lim = c(0, 500), y.lim = c(0, 500), 
                   system = "two.dim", add = F)
traj <- trajectory(logistic, y0 = seq(1, 200, 2), t.start = 0, t.end = 100, t.step = 1, 
                   system = "one.dim", lty = 1:10 , colour = 1:10, add = T)
pp <- phasePortrait(logistic, y.lim = c(-10, 120))


model_abio <- function(t, y, parameters) {
  a <- y[1:100]
  X <- y[101]
  ext <- k*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1-rho)*phi*a^(alpha+delta)*X^gamma 
  da <- a^(alpha+alpha+delta)*s*k^2*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)^2*phi*X^gamma - a*mu
  dX <- c*(1 - X/m) - sum(ext)
  list(c(da, dX))
}

out <- ode(iv, t, model_abio, NULL)
matplot(out[,1], out[,2:101], type = 'l')
matplot(out[,1], out[,102], type = 'l')

#Model biotic----
n       <- 100
alpha   <- 0.2
c       <- 30
delta   <- 0.2
gamma   <- 0.3
k       <- 1
lambda1 <- 20
m       <- 100
mu      <- 0.1
n       <- 100
phi     <- 1
rho     <- 0.5
s       <- 0.25
theta1  <- 50
tau     <- 10
chi     <- 1e-4
iv <- c(a = 100, X = 100)
t = seq(0, 100, 1)

model_bio_2d <- function(t, y, parameters) {
  a <- y[1]
  X <- y[2]
  ext <- k*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)*a^(alpha + delta)*X^gamma 
  da <- a^(alpha + alpha + delta)*s*k^2*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)^2*X^gamma - a*mu
  dX <- tau*(X - b)*(X - c)*(1 - X/m) - n*ext
  list(c(da, dX))
}

out <- ode(iv, t, model_bio_2d, NULL)
matplot(out[,1], out[,-1], type = 'l')

ff <- flowField(model_bio_2d, x.lim = c(0, 1000), y.lim = c(0, 1000), points = 20, 
                system = 'two.dim', add = F, xlab = 'a', ylab = 'X')
null <- nullclines(model_bio_2d, x.lim = c(0, 2000), y.lim = c(0, 600), 
                   system = "two.dim", add = F, xlab = 'Wealth (a)', ylab = 'Resource (X)')
lines(out[,2], out[,3], lwd = 2)
traj <- trajectory(model_bio_2d, y0 = iv, t.start = 0, t.end = 1, t.step = 0.01, 
                   system = "two.dim", add = T)
pp <- phasePortrait(logistic, y.lim = c(-10, 120))


model_bio <- function(t, y, parameters) {
  a <- y[1:100]
  X <- y[101]
  ext <- k*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)*phi*a^(alpha + delta)*X^gamma 
  da <- a^(alpha + alpha + delta)*s*k^2*((rho*a^lambda1)/(a^lambda1 + theta1^lambda1) + 1 - rho)^2*phi*X^gamma - a*mu
  dX <- tau*(X - b)*(X - c)*(1 - X/m) - sum(ext)
  list(c(da, dX))
}

out <- ode(iv, t, model_bio, NULL)
matplot(out[,1], out[,2:101], type = 'l')
matplot(out[,1], out[,102], type = 'l')


#Stability----
iv <- c(a = 4e4, X = 100)
points(iv[1], iv[2], pch = 20)
roots <- steady(y = iv, func = model_bio_2d_total)
points(roots$y[1], roots$y[2], pch = 20, col = 'red')
c_point <- steady(y = iv, func = model_bio_2d_rate, method = "runsteady", time = c(0, 1e2))
points(roots$y[1], roots$y[2], pch = 20, col = 'blue')
eigen(jacobian.full(y = c(roots$y[1], roots$y[2]), func =  model_bio_2d_total))

par <- seq(0.25, 0.7, 0.001)
stab <- matrix(NA, length(par), 5); colnames(stab) <- c('Par', 'a', 'X', '1st eigen value', '2nd eigen value')
iv <- c(a = 29738.5030, X = 93.6035)
for (i in 1:length(par)) {
  s <- par[i]
  stab[i, 1] <- par[i]
  roots <- steady(y = iv, func = model_bio_2d_total, positive = T)
  iv <- c(a = roots$y[1], X = roots$y[2])
  stab[i, c(2, 3)] <- c(roots$y[1], roots$y[2])
  stab[i, c(4, 5)] <- eigen(jacobian.full(y = c(roots$y[1], roots$y[2]), func =  model_bio_2d_total))$values
}
null <- nullclines(model_bio_2d, x.lim = c(0, 3000), y.lim = c(0, 200), system = "two.dim", 
                   add = F, xlab = 'Wealth (a)', ylab = 'Resource (X)', panel.first = grid())
points(stab[,c(2, 3)], pch = 20, col = 'red')

#Bifurcations----
data <- read.table('b.t_1db', fill = T, header = T, skip = 18)
data <- data[data[,1] == 1, ]
data <- as.data.frame(lapply(data, function(f) {as.numeric(as.character(f))}))
data <- data[1:5412,]

  #for X wrt s
ggplot() + geom_path(data = data[1:which(data$TY == 2)[1],], aes(x = s, y = x), size = 1) +
  geom_path(data = data[which(data$TY == 9)[2]:nrow(data),], aes(x = s, y = x), size = 1) +
  geom_path(data = data[which(data$TY == 2)[1]:which(data$TY == 2)[2],], aes(x = s, y = x), colour = 'red', linetype = 2, size = 1) + 
  geom_path(data = data[which(data$TY == 2)[2]:which(data$TY == 2)[3],], aes(x = s, y = x), size = 1) + 
  geom_path(data = data[which(data$TY == 2)[3]:which(data$TY == 2)[4],], aes(x = s, y = x), colour = 'red', linetype = 2, size = 1) +
  geom_path(data = data[which(data$TY == 2)[4]:which(data$TY == 2)[5],], aes(x = s, y = x), colour = 'red', linetype = 2, size = 1) + 
  geom_path(data = data[which(data$TY == 2)[5]:which(data$TY == 2)[6],], aes(x = s, y = x), colour = 'red', linetype = 2, size = 1) + 
  geom_path(data = data[which(data$TY == 2)[6]:which(data$TY == -9),], aes(x = s, y = x), colour = 'red', linetype = 2, size = 1) + 
  geom_point(data = data[which(data$TY == 2),], aes(x = s, y = x), colour = 'red', size = 2) + 
  theme_light()

  #for a wrt s
ggplot() + geom_path(data = data[1:which(data$TY == 2)[1],], aes(x = s, y = a), size = 1) +
  geom_path(data = data[which(data$TY == 9)[2]:nrow(data),], aes(x = s, y = a), size = 1) +
  geom_path(data = data[which(data$TY == 2)[1]:which(data$TY == 2)[2],], aes(x = s, y = a), colour = 'red', linetype = 2, size = 1) + 
  geom_path(data = data[which(data$TY == 2)[2]:which(data$TY == 2)[3],], aes(x = s, y = a), size = 1) + 
  geom_path(data = data[which(data$TY == 2)[3]:which(data$TY == 2)[4],], aes(x = s, y = a), colour = 'red', linetype = 2, size = 1) +
  geom_path(data = data[which(data$TY == 2)[4]:which(data$TY == 2)[5],], aes(x = s, y = a), colour = 'red', linetype = 2, size = 1) + 
  geom_path(data = data[which(data$TY == 2)[5]:which(data$TY == 2)[6],], aes(x = s, y = a), colour = 'red', linetype = 2, size = 1) + 
  geom_path(data = data[which(data$TY == 2)[6]:which(data$TY == -9),], aes(x = s, y = a), colour = 'red', linetype = 2, size = 1) + 
  geom_point(data = data[which(data$TY == 2),], aes(x = s, y = a), colour = 'red', size = 2) + 
  geom_hline(yintercept = 0, size = 1) +
  theme_light()

data_2d <- read.table('b.t_2db', fill = T, header = T, skip = 19)
data_cusp <- data_2d[data_2d[,1] == 2, ]
data_cusp <- as.data.frame(lapply(data_cusp, function(x) {as.numeric(as.character(x))}))
data_lp <- data_2d[data_2d[,1] == 8, ]
data_lp <- as.data.frame(lapply(data_lp, function(x) {as.numeric(as.character(x))}))
poly_1 <- data.frame( x = c(0, 0, 0.1356484, 0.1261166, 0), y = c(0, 1, 1, 0, 0))

ggplot() + geom_path(data = data_cusp[1:which(data_cusp$TY == -24)[1],], aes(x = s, y = rho), size = 1) + 
  geom_path(data = data_cusp[which(data_cusp$TY == -24)[1]:which(data_cusp$TY == -24)[2],], aes(x = s, y = rho), size = 1) +
  geom_vline(xintercept =  data_lp$s[1], size = 1) +
  geom_point(data = data_cusp[which(data_cusp$TY == 22),], aes(x = s, y = rho), colour = 'red', size = 2) +
  theme_light() + coord_cartesian(c(0, 1), c(0, 1))

#Rate induced tipping----
n       <- 100
alpha   <- 0.2
c       <- 30 #100
delta   <- 0.2
gamma   <- 0.3
k       <- 1
lambda1 <- 20
m       <- 100 #300
mu      <- 0.1
n       <- 100
rho     <- 0.5
theta1  <- 50 #100
tau     <- 10   #2
s       <- 0.25
sm <- 0.9 #HB point is 0.654 and SN point is 0.657
s0 <- 0.25
r1 <- 5e-4
r2 <- 5e-4
iv <- c(a = 29738.5030, X = 93.6035, s = s0)
t = seq(0, 3e3, 1)

model_bio_2d_rate <- function(t, y, parameters) {
  a <- y[1]
  X <- y[2]
  s <- y[3]
  if (X < 0 | is.na(X)) {X <- 0}
  if (a < 0 | is.na(a)) {a <- 0}
  t1 <- (sm - s0)/r1
  t2 <- (sm - s0)/r2 + t1
  if (t >= 0 & t < t1) {r <- r1}
  else if (t >= t1 & t < t2) {r <- -r2}
  else {r <- 0}
  ext <- k*((rho*(a/n)^lambda1)/((a/n)^lambda1 + theta1^lambda1) + 1 - rho)*(a/n)^(alpha + delta)*X^gamma
  ds <- r
  da <- n*(a/n)^(alpha + alpha + delta)*s*k^2*((rho*(a/n)^lambda1)/((a/n)^lambda1 + theta1^lambda1) + 1 - rho)^2*X^gamma - a*mu
  dX <- tau*(X)*(X - c)*(1 - X/m) - n*ext
  list(c(da, dX, ds))
}

out_r <- ode(iv, t, model_bio_2d_rate, NULL)
matplot(out_r[,'time'], out_r[,c('s')], type = 'l', xlab = 'Time', ylab = 'Savings', lwd = 2); grid()
matplot(out_r[,'time'], out_r[,c('a')], type = 'l', xlab = 'Time', ylab = 'Wealth', lwd = 2); grid()
matplot(out_r[,'time'], out_r[,c('X')], type = 'l', xlab = 'Time', ylab = 'Resource', lwd = 2); grid()

rates <- expand.grid(r1 = seq(5e-4, 5e-2, length.out = 50), r2 = seq(5e-4, 5e-2, length.out = 50))
rate_tip <- matrix(NA, nrow(rates), 3, dimnames = list(NULL, c('a','X', 's')))
for (i in 1:nrow(rates)) {
  cat("Iteration", i, "\n")
  r1 <- rates[i, 1]
  r2 <- rates[i, 2]
  out <- ode(iv, t, model_bio_2d_rate, NULL)
  rate_tip[i,] <- out[nrow(out),2:4]
}
rate_tip <- cbind(rates, rate_tip) 
rate_tip[is.na(rate_tip)] <- 0
ggplot(data = rate_tip, aes(x = r1, y = r2, fill = a)) + geom_tile() + 
  scale_fill_gradient(low = 'lightblue', high = 'darkblue', breaks = c(1, 2.9e4), labels = c('Tipping to collapse', 'Return to safety'), name = '') + 
  theme_minimal()

ggplot(data = rate_tip, aes(x = r1, y = r2, z = a)) + geom_contour(bins = 1.5) + 
  scale_fill_gradient(low = 'lightblue', high = 'darkblue', breaks = c(1, 98), labels = c('Tipping to collapse', 'Return to safety'), name = '') + 
  theme_minimal() 

sm1 <- seq(0.655, 0.66, length.out = 6)
r2_c1 <- c(1.66265e-05, 5.704450e-05, 0.0001077848326, 0.0001718082, 0.0002532070, 0.0003590450) 
plot(sm1, r2_c1, pch = 19, col = 'red')
lines(sm1, r2_c1, col = 'red')

#Static equilibrum----
model_bio_2d_total <- function(t, y, parameters) {
  a <- y[1]
  X <- y[2]
  ext <- k*((rho*(a/n)^lambda1)/((a/n)^lambda1 + theta1^lambda1) + 1 - rho)*(a/n)^(alpha + delta)*X^gamma
  da <- n*(a/n)^(alpha + alpha + delta)*s*k^2*((rho*(a/n)^lambda1)/((a/n)^lambda1 + theta1^lambda1) + 1 - rho)^2*X^gamma - a*mu
  dX <- tau*(X)*(X - c)*(1 - X/m) - n*ext
  list(c(da, dX))
}
t = seq(0, 1e3, 1)
iv <- c(a = 6e4, X = 200)
out_t <- ode(iv, t, model_bio_2d_total, NULL)
plot(out_t)
null <- nullclines(model_bio_2d_total, x.lim = c(0, 1e5), y.lim = c(0, 600), 
                   system = "two.dim", add = F, xlab = 'Wealth (a)', ylab = 'Resource (X)')
iv <- c(a = 4e4, X = 100)
points(iv[1], iv[2], pch = 20)
roots <- steady(y = iv, func = model_bio_2d_total)
points(roots$y[1], roots$y[2], pch = 20, col = 'red')
lines(out_t[,'a'], out_t[,'X'], lwd = 2)

#Tracking the dynamic equilibrium---- 
par <- out_r[,c('s')]; par[par > 0.654] <- 0.654
ss <- matrix(NA, length(par), 2, dimnames = list(NULL, c('a','X')))
iv <- c(a = 29738.5030, X = 93.6035)
roots <- NA
for (i in 1:length(par)) {
  s <- par[i]
  roots <- steady(y = iv, func = model_bio_2d_total, positive = T)
  ss[i,] <- c(roots$y[1], roots$y[2])
  iv <- c(roots$y[1], roots$y[2])
  #roots <- ode(iv, t, model_bio_2d_total, NULL)[1001, 2:3]
  #ss[i,] <- c(roots[1], roots[2])
}
ss <- as.data.frame(cbind(par, ss, out_r)); ss[ss[,1] == 0.654, 2:3] <- NA
names(ss) <- c("par",  "a_st",    "X_st",    "time", "a_dy",    "X_dy",    "s_dy"  )

ggplot(data = ss, aes(x = time)) + geom_line(aes(y = X_st, color = 'Static State'), size = 1) + 
  geom_line(aes(y = X_dy, color = 'Dynamic State'), size = 1) +
  scale_color_discrete(name = '') + coord_cartesian(c(0, 2500), c(60, 100)) +
  theme_light()

ggplot(data = ss, aes(x = time)) + geom_line(aes(y = a_st, color = 'Static State'), size = 1) + 
  geom_line(aes(y = a_dy, color = 'Dynamic State'), size = 1) + 
  scale_color_discrete(name = '') + coord_cartesian(c(0, 2500)) +
  theme_light()


