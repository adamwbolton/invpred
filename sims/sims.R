# single sceneario
library(invpred)
library(doParallel)
library(foreach)
library(knitr)
library(kableExtra)

arsenic <- stats::lm(measured ~ actual, data = investr::arsenic)

(cal <- investr::calibrate(arsenic, y0 = 3, interval = "inversion"))
(cal <- investr::calibrate(arsenic, y0 = 3, interval = "Wald"))

sd1 <- sd2 <- sigma(arsenic)
y0 = c(6,6.1)
invpred(arsenic, y0, interval = "t")
y0 = c(3,3.1)
invpred(arsenic, y0, interval = "welch")

y0 = rnorm(1e3, mean = 3, sd = sd2)
invpred(arsenic, y0, interval = "t")
invpred(arsenic, y0, interval = "bf")
invpred(arsenic, y0, interval = "welch")

### simlulate data ####
out <- do_sim_multiple(object = arsenic, x0 = seq(1,6, by = 1), m = 3, sd1 = 1e-1, sd2 = 1.1e-1,sims = 1e2)
res <- out$results


sum <- res %>% group_by(interval, x0) %>%
  summarise(cov = mean((lower < x0) * (upper > x0)), lwr = mean(lower), upr = mean(upper))


sum %>%
  ggplot(aes(x0, cov, colour = interval)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  lims(y = c(0,1)) +
  geom_hline(yintercept = 0.95, linetype = 2)

# coverage
sum %>%
  mutate(sd1 = out$sd1, sd2 = out$sd2, m = out$m) %>%
  select(sd1, sd2, interval, cov, x0,m) %>%
  mutate(val = paste(format(cov,digits = 3))) %>%
  select(sd1,sd2, interval, val, x0) %>%
  arrange(sd1,sd2,interval) %>%
  pivot_wider(names_from = x0, values_from = val) %>%
  kable(format = "latex", digits = 4, longtable = T)


sum %>%
  ggplot(aes(x0, y = upr - lwr, colour = interval)) +
  geom_point() +
  geom_line() +
  theme_bw()





# multiple scenarios


sd1s <- c(1e-2,1e-3)
sd2s <- c(1e-2,1e-3)
ms <- c(2,5, 10, 20)
params <- expand.grid(sd1 =sd1s,sd2 =sd2s, m = ms)
x0s <- seq(1,6, by = 0.5)
iters <- 1e2

out <- foreach(i = 1:nrow(params), .combine = rbind) %do%
  (do_sim_multiple(object = arsenic, x0 = x0s, m = params[i,'m'], sd1 = params[i,'sd1'], sd2 = params[i,'sd2'],
                  sims = 1e2)$results %>% mutate(sd1 = params[i,'sd1'], sd2 = params[i,'sd2'], m = params[i,'m']))


# print results
print_x0s <- seq(1,6)
sum <- out %>% group_by(interval, x0, sd1, sd2, m) %>%
  summarise(cov = mean((lower < x0) * (upper > x0)), lwr = mean(lower), upr = mean(upper), wid = mean(upper - lower), sdwid = sd(upper -lower)) %>%
  filter(x0 %in% print_x0s)


covtb <- sum %>%
  mutate(cov = paste(format(cov,digits = 4))) %>%
  select(sd1,sd2,m,interval, cov, x0) %>%
  arrange(sd1,sd2,m,interval) %>%
  pivot_wider(names_from = x0, values_from = cov)
colnames(covtb) <- c("$\\sigma_1$", "$\\sigma_2$", "m", "Method", print_x0s)

kable(covtb, align = "c", "latex",escape = FALSE, booktabs = T, longtable = T) %>%
  kable_styling(latex_options = c("repeat_header"), font_size = 7) %>%
  kable_paper(full_width = F) %>%
  kable_paper() %>%
  column_spec(1:4, bold = T) %>%
  collapse_rows(columns = 1:3, valign = "middle",latex_hline = "custom", custom_latex_hline = c(1,2,3)) %>%
  write(file = "~/Desktop/cov.tex")




# kable(covtb, align = "c", "latex",escape = FALSE, booktabs = T) %>%
#   kable_paper(full_width = F) %>%
#   column_spec(1:4, bold = T) %>%
#   collapse_rows(columns = 1:3, valign = "top",custom_latex_hline = c(1,2,3,4))


widtb <- sum %>%
  ungroup() %>%
  mutate(val = paste0(format(wid,digits = 3), " (", format(sdwid,digits =3), ")")) %>%
  select(sd1,sd2, m,interval,x0,val) %>%
  arrange(sd1,sd2,m,interval) %>%
  pivot_wider(names_from = x0, values_from = c(val))
# widtb[,1:5]

colnames(widtb) <- c("$\\sigma_1$", "$\\sigma_2$", "m", "Method", print_x0s)
kable(widtb, align = "c", "latex",escape = FALSE, booktabs = T, longtable = T) %>%
  kable_styling(latex_options = c("repeat_header"), font_size = 7) %>%
  kable_paper() %>%
  column_spec(1:4, bold = T) %>%
  collapse_rows(columns = 1:3, valign = "middle",latex_hline = "custom", custom_latex_hline = c(1,2,3)) %>%
  write(file = "~/Desktop/wid.tex")





# graphs
sum %>%
  mutate(m = factor(m)) %>%
  ggplot(aes(x0, cov, colour = m, linetype = interval)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_grid(rows = vars(sd1), cols = vars(sd2)) +
  lims(y = c(0,1)) +
  geom_hline(yintercept = 0.95, linetype = 2)


sum %>%
  mutate(m = factor(m)) %>%
  ggplot(aes(x0, upr - lwr, colour = m, linetype = interval)) +
  geom_line() +
  theme_bw() +
  facet_grid(rows = vars(sd1), cols = vars(sd2))


out %>%
  mutate(m = factor(m)) %>%
  group_by(sd1, sd2, m, interval, x0) %>%
  summarise(wid = mean(upper - lower), sdwid = sd(upper - lower)) %>%
  ggplot(aes(x0, wid, colour = m, linetype = interval)) +
  geom_errorbar(aes(ymin=wid-sdwid, ymax=wid+sdwid), width=.1) +
  theme_bw() +
  facet_grid(rows = vars(sd1), cols = vars(sd2))
