#Test simulation and likelihood functions behave as expected
context("Test simulation and likelihood functions behave as expected")

test_that("Simulation and likelihood functions agree on likelihood value, using no phylo object and no topology.", {
  set.seed(0)
  t=c(2000.4, 2000.8, 2001.2, 2001.6, 2002)
  l=rep(2, 5)
  ne=1.1
  b=2000
  expect_silent(sim<-bounded_sample_times(t=t,l=l , ne = ne, b = b))
  expect_is(sim,'list')
  expect_silent(lik<-bounded_likelihood(t=t,l=l,c=sim$times,ne=ne,b=b,topology = F))
  expect_is(lik,'numeric')
  expect_equal(sim$likelihood,lik)
})

test_that("Likelihood function gives expected result on small unbounded example.", {
  t=c(2000,2001)
  l=c(2,1)
  c=c(1999,1998)
  ne=2
  b=-Inf
  a=bounded_likelihood(t=t,l=l,c=c,ne=ne,b=b,topology = F)
  b=dexp(1,3/2)*dexp(1,1/2)
  expect_equal(a,b)
})

test_that("Simulation and likelihood functions agree on likelihood value, using phylo object and topology.", {
  set.seed(1)
  t=c(2000.4, 2000.8, 2001.2, 2001.6, 2002)
  l=rep(3, 5)
  ne=1.1
  b=2000
  expect_silent(sim<-bounded_sample_phylo(t=t,l=l , ne = ne, b = b))
  expect_is(sim,'list')
  expect_silent(lik<-bounded_likelihood_phylo(phy=sim$phylo,ne=ne,b=b,topology = T))
  expect_is(lik,'numeric')
  expect_equal(sim$likelihood,lik)
})

test_that("Likelihood function works on a rtree.",{
  library(ape)
  set.seed(2)
  t=rtree(100)
  t$root.time=2000
  expect_silent(lik<-bounded_likelihood_phylo(phy=t,ne=0.9,b=1999,topology = T))
  expect_is(lik,'numeric')
})
