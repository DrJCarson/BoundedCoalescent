#Test simulation and likelihood functions behave as expected
context("Test simulation and likelihood functions behave as expected")

test_that("Simulation and likelihood functions agree on likelihood value.", {
  set.seed(0)
  t=c(0.4, 0.8, 1.2, 1.6, 2)
  l=rep(2, 5)
  ne=1.1
  b=-0.1
  expect_silent(sim<-bounded_times_sample(t=t,l=l , ne = ne, b = b, nsam = 1))
  expect_is(sim,'list')
  expect_silent(l<-bounded_times_likelihood(t=t,l=l,ne=ne,b=b,c=sim$times))
  expect_is(l,'numeric')
  expect_equal(sim$likelihood,l)
})

test_that("Likelihood function gives expected result on small homochronous unbounded example.", {
  t=c(2000,2000,2000)
  l=c(1,1,1)
  c=c(1999,1998)
  ne=2
  b=0
  a=bounded_times_likelihood(t=t,l=l,c=c,ne=ne,b=b)
  b=dexp(1,3/2)*dexp(1,1/2)
  expect_equal(a,b)
})
