library(latex2exp)
library(ggplot2)
library(patchwork)

################
### Figure 1 ###
################
x <- c(2, 4, 6, 8)
N <- c(1, 2, 5, 10)
y <- lapply(1:4, function(i) x[i] * rnorm(N[i]))
# load('data/y.RData')
y <- lapply(1:4, function(i) {
  sort(y[[i]])
}) # sort observed values
M <- plcm(N) # least common multiple of N_i
yM <- t(sapply(1:4, function(i) {
  rep(y[[i]], each = M / N[i])
  })) # n by M

df <- data.frame(y = rep(x, N),
                 x = unlist(y), 
                 subj = factor(rep(1:4, N)))
# 1 observations
g1 <- ggplot() +
  geom_point(data = df, aes(x = x, y = y, color = subj, shape = subj), size = 2) +
  geom_line(data = data.frame(y = c(sapply(x, function(xi) xi + 9 * dnorm(seq(-20, 20, 0.1), mean = 0, sd = xi))), 
                              x = rep(seq(-20, 20, 0.1), 4), 
                              subj = factor(rep(1:4, each = 401))), aes(x = x, y = y, color = subj, group = subj)) +
  coord_flip() +
  labs(x = TeX("$Y_{ij}$"), y = TeX('$Z_i$')) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=20), axis.title.y=element_text(angle = 0, vjust = 0.5))

# 2 empirical quantile functions
g2 <- ggplot(df, aes(x, color = subj)) +
  stat_ecdf(geom = "step") +
  # geom_hline(yintercept = 1:10/10, linetype = 2, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = '', limits = c(-20, 20)) +
  scale_y_continuous(name = 'Probability') +
  theme_bw() +
  theme(legend.position = 'none', text = element_text(size=20))

# barycenter and REM
res1 <- grem(y, x, 5)
res2 <- grem(y, x, 3)
comp <- c('z=5', 'z=3')
g3 <- ggplot(data.frame(x = c(res1$qp, res2$qp), subj = factor(rep(comp, each = 10), levels = comp[2:1])), aes(x, linetype = subj)) +
  stat_ecdf(geom = "step") +
  # geom_hline(yintercept = 1:10/10, linetype = 2, alpha = 0.5) +
  coord_flip() +
  labs(x = '') +
  scale_y_continuous(name = 'Probability') +
  scale_linetype_discrete(labels = c(TeX(r"($\hat{m}_G(3)$)"), TeX(r"($\hat{m}_G(5)$)"))) +
  theme_bw() +
  theme(legend.position = 'top', legend.title = element_blank(), text = element_text(size=20))

# density
dFit <- list(density(res1$qp, n = 1000, bw = 3.5), 
             density(res2$qp, n = 1000, bw = 2.5))
df <- data.frame(d = c(sapply(dFit, function(dFiti) dFiti$y)), 
                 dsup = c(sapply(dFit, function(dFiti) dFiti$x)), 
                 comp = factor(rep(comp, each = 1000), levels = comp[2:1]))

g4 <- ggplot(data = df) + 
  geom_line(aes(x = dsup, y = d, linetype = comp)) +# , color = '#F8766D'
  labs(x = '') + 
  scale_y_continuous(name = 'Density', limits = c(0, 0.15)) +
  scale_linetype_discrete(labels = c(TeX(r"($\hat{m}_G(3)$)"), TeX(r"($\hat{m}_G(5)$)"))) +
  theme_bw() +
  theme(legend.position = 'top', legend.title = element_blank(), text = element_text(size=20))

pdf(file = "latex/Biometrika/img/diagram.pdf", width = 14, height = 8)
(g1 | g2) / (g3 | g4)
dev.off()

################
### Figure 2 ###
################
Q <- 100# number of simulations
nVec <- c(50, 100, 200, 500, 1000)# number of subjects

load('data/iseg.RData')
dfg <- data.frame(err = c(iseg0, isegp, iseg, iseg0T, isegpT, isegT), 
                  method = factor(rep(rep(c('Fully observed', 'Petersen and Müller (2019)', 'REM'), each = Q * length(nVec)), 2), levels = c('Fully observed', 'Petersen and Müller (2019)', 'REM')), 
                  groupn = rep(rep(nVec, each = Q), 6), 
                  gl = rep(c('Setting I', 'Setting III'), each = Q * length(nVec) * 3))
load('data/isel.RData')
dfl <- data.frame(err = c(isel0, iselp, isel, isel0T, iselpT, iselT), 
                  method = factor(rep(rep(c('Fully observed', 'Petersen and Müller (2019)', 'REM'), each = Q * length(nVec)), 2), levels = c('Fully observed', 'Petersen and Müller (2019)', 'REM')), 
                  groupn = rep(rep(nVec, each = Q), 6), 
                  gl = rep(c('Setting II', 'Setting IV'), each = Q * length(nVec) * 3))
df <- rbind(dfg, dfl)
pdf(file = 'latex/img/ise.pdf', width = 12, height = 10)
ggplot(data = df, aes(x = factor(groupn), y = err, fill = method, linetype = method)) +
  facet_wrap(vars(gl), dir = 'v', scales = 'free') +
  geom_boxplot(notch = TRUE, outlier.alpha = 0.2) +
  scale_linetype_manual(values = c('dashed', 'longdash', 'solid')) +
  labs(x = 'n', y = 'ISE') +
  theme_bw() +
  theme(legend.position = 'top', legend.title = element_blank(), text = element_text(size=20))
dev.off()

################
### Figure 3 ###
################
load('data/isegb.RData')
dfg <- data.frame(err = c(iseg0, isegp, iseg, iseg0T, isegpT, isegT), 
                  method = factor(rep(rep(c('Fully observed', 'Petersen and Müller (2019)', 'REM'), each = Q * length(nVec)), 2), levels = c('Fully observed', 'Petersen and Müller (2019)', 'REM')), 
                  groupn = rep(rep(nVec, each = Q), 6), 
                  gl = rep(c('Setting I', 'Setting III'), each = Q * length(nVec) * 3))
load('data/iselb.RData')
dfl <- data.frame(err = c(isel0, iselp, isel, isel0T, iselpT, iselT), 
                  method = factor(rep(rep(c('Fully observed', 'Petersen and Müller (2019)', 'REM'), each = Q * length(nVec)), 2), levels = c('Fully observed', 'Petersen and Müller (2019)', 'REM')), 
                  groupn = rep(rep(nVec, each = Q), 6), 
                  gl = rep(c('Setting II', 'Setting IV'), each = Q * length(nVec) * 3))
df <- rbind(dfg, dfl)
pdf(file = 'latex/img/iseb.pdf', width = 12, height = 10)
ggplot(data = df, aes(x = factor(groupn), y = err, fill = method, linetype = method)) +
  facet_wrap(vars(gl), dir = 'v', scales = 'free') +
  geom_boxplot(notch = TRUE, outlier.alpha = 0.2) +
  scale_linetype_manual(values = c('dashed', 'longdash', 'solid')) +
  labs(x = 'n', y = 'ISE') +
  theme_bw() +
  theme(legend.position = 'top', legend.title = element_blank(), text = element_text(size=20))
dev.off()
