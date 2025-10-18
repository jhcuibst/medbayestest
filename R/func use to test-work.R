usethis::create_package("d:/13-Git/mediationBayes")

library(medbayestest)
getNamespaceExports("medbayestest")
find.package("medbayestest")
devtools::document()
devtools::install(clean = TRUE, upgrade = "never")
medbayestest::medbayes()

devtools::install_github("jhcuibst/medbayestest")

# remove.packages("medbayestest")
# find.package("medbayestest")
library(mediationBayes)
library(brms)
data(example_data)
prior.m = c(
  set_prior("student_t(3, 0, 2.5)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "hu") )

prior.y = set_prior("student_t(3, 0, 2.5)", class="b")


hnb.m <- brm(bf(mnb~ x + age, hu ~ x + age), data = example_data, prior = prior.m,
             family = hurdle_negbinomial(), chains = 4, iter = 2000)

hnb.y <- brm(y ~  x + mnb + age + im, data = example_data, prior = prior.y,
             family = bernoulli(), chains = 4, iter = 2000)


outmed <- medbayestest::medbayes(model.m = hnb.m, model.y = hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",
                   control.value = 0, treat.value = 1)

outmed2 <- medbayes_zim(model.m = hnb.m, model.y = hnb.y, mediator="mnb", treat="x", outcome = "y", ind_mediator = "im",
             control.value = 0, treat.value = 1)

outmed2$effects.rr




