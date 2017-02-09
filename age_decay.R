# Run lines 161-176 and 184-186 of resistance_regressions.R first
require(ggplot2)

env_lm_coef <- summary(lm(lasso_y ~ filtered_X$age_d[-outliers] + filtered_X$carried[-outliers]))$coefficients
x = seq(0, 365*10, 30)
y = exp(env_lm_coef[1,1])*exp(env_lm_coef[2,1]*x)

dat = cbind(x,y)

ggplot(as.data.frame(dat), aes(x=x, y=y)) + 
  geom_smooth(method = "glm", method.args = list(family = gaussian(lin="log"))) + 
  xlab("Age (days)") + ylab("Mean duration of carriage episode (days)") + 
  xlim(0,365*3) + ylim(0,70) + theme_bw(base_size = 16)
