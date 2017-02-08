# Run lines 161-176 and 184-186 of resistance_regressions.R first

summary(lm(lasso_y ~ filtered_X$age_d[-outliers] + filtered_X$carried[-outliers]))
env_lm_coef <- summary(lm(lasso_y ~ filtered_X$age_d[-outliers] + filtered_X$carried[-outliers]))$coefficients
env_lm_coef
x = seq(0, 365*10, 30)
x
y = exp(env_lm_coef[0,0])*exp(env_lm_coef[1,0]*x)
plot(x,y)
y
exp(env_lm_coef[0,0])
env_lm_coef
env_lm_coef[0]
env_lm_coef[,0]
env_lm_coef[,1]
env_lm_coef[0,1]
env_lm_coef[1,1]
y = exp(env_lm_coef[1,1])*exp(env_lm_coef[2,1]*x)
y
plot(x,y)
env_lm <- lm(lasso_y ~ filtered_X$age_d[-outliers])
x
predict(env_lm)
predict(env_lm, data=seq(0, 365*10, 30))
?predict
predict(env_lm, seq(0, 365*10, 30))
env_lm
require(ggplot2)
dat = cbind(x,y)
dat
ggplot(dat, aes(x=x, y=y)) + geom_point()
ggplot(as.data.frame(dat), aes(x=x, y=y)) + geom_point()
ggplot(as.data.frame(dat), aes(x=x, y=y)) + geom_smooth(method = "glm", family = gaussian(lin="log"))
ggplot(as.data.frame(dat), aes(x=x, y=y)) + geom_smooth(method = "glm", method.args = list(family = gaussian(lin="log")))
ggplot(as.data.frame(dat), aes(x=x, y=y)) + geom_smooth(method = "glm", method.args = list(family = gaussian(lin="log"))) + xlab("Age in days") + ylab("Mean duration of carriage episode")
ggplot(as.data.frame(dat), aes(x=x, y=y)) + geom_smooth(method = "glm", method.args = list(family = gaussian(lin="log"))) + xlab("Age in days") + ylab("Mean duration of carriage episode") + xlim(0,365*3)
ggplot(as.data.frame(dat), aes(x=x, y=y)) + geom_smooth(method = "glm", method.args = list(family = gaussian(lin="log"))) + xlab("Age in days") + ylab("Mean duration of carriage episode") + xlim(0,365*3) + theme_bw(base_size = 16)
