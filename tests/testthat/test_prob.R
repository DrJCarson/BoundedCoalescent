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

test_that("Simulation can be done using either direct method or rejection sampling.", {
  set.seed(3)
  t=c(2000.4, 2000.8, 2001.2, 2001.6, 2002)
  l=rep(3, 5)
  ne=1.1
  b=2000
  expect_silent(sim<-bounded_sample_phylo(t=t,l=l , ne = ne, b = b,method='direct'))
  expect_is(sim,'list')
  expect_gt(bounded_likelihood_phylo(phy=sim$phylo,ne=0.9,b=1999,topology = T),0)
  expect_silent(sim<-bounded_sample_phylo(t=t,l=l , ne = ne, b = b,method='rejection'))
  expect_is(sim,'list')
  expect_gt(bounded_likelihood_phylo(phy=sim$phylo,ne=0.9,b=1999,topology = T),0)
})

test_that("Likelihood function works for a multiPhylo.",{
  library(ape)
  set.seed(2)
  t=c(rtree(10),rtree(20))
  t[[1]]$root.time=t[[2]]$root.time=2000
  expect_silent(lik<-bounded_likelihood_phylo(phy=t,ne=0.9,b=1999,topology = T))
  expect_is(lik,'numeric')
  expect_equal(length(lik),2)
})

test_that("Lineage number remains constant if dt=0.",{
  expect_equal(coalescent_probability(10,10,0,2),1)
})


test_that("Forward probabilities add up to 1.",{
  expect_silent(b<-bounded_forward_algorithm(2010:2019,rep(1,10),1.2,2000))
  expect_equal(max(abs(colSums(b$probabilities)-1)),0)
})


test_that("Can simulate with rejection.",{
  set.seed(0)
  t=2000:2020
  l=rep(2, length(t))
  ne=1.3
  b=1999
  expect_silent(sim<-bounded_sample_times(t=t,l=l , ne = ne, b = b,nsam=5,method='rejection'))
  expect_is(sim,'list')
  expect_equal(nrow(sim$times),5)
  expect_silent(sim<-bounded_sample_phylo(t=t,l=l , ne = ne, b = b,nsam=5,method='rejection'))
  expect_is(sim,'list')
  expect_is(sim$phylo,'multiPhylo')
  expect_equal(length(sim$phylo),5)
})
