#include "gtest/gtest.h"
#include "TestVariables.hpp"

class Algebra : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    gen = std::mt19937(rd());

  }

  std::random_device rd;
  std::mt19937 gen;

  prob::distribution<double,A> pA, qA;
  prob::distribution<double,B,C> pBC, qBC;
  prob::distribution<double,A,B,C> pABC, qABC;
  prob::distribution<double,A, prob::given, B, C> pAgBC, qAgBC;
  prob::distribution<double,B, C, prob::given, A> pBCgA;
  prob::distribution<double,A, B, C, prob::given, D> pABCgD;
  prob::distribution<double, A, prob::given, B, C,D> pAgBCD, qAgBCD;
  prob::distribution<double, B, C, prob::given, D> pBCgD;
};

TEST_F(Algebra, JoinMarginalize)
{
  prob::init::random(pA, gen);
  prob::init::random(pBC, gen);

  pABC = prob::join(pA,pBC);
  qA = pABC.marginalize<0>();
  qBC = pABC.marginalize<1,2>();

  EXPECT_LT((pA-qA).array().abs().sum(), 1e-10);
  EXPECT_LT((pBC-qBC).array().abs().sum(),1e-10);
}

TEST_F(Algebra, UnconditionCondition)
{
    prob::init::random(pAgBC, gen);
    prob::init::random(pBC, gen);

    pABC = prob::uncondition(pAgBC, pBC);
    prob::condition(pABC, pBC, qAgBC);

    EXPECT_LT((qAgBC-pAgBC).array().abs().sum(), 1e-10);
}

TEST_F(Algebra, PartialUnconditionCondition)
{
  prob::init::random(pAgBCD, gen);
  prob::init::random(pBCgD, gen);

  pABCgD = prob::partial_uncondition(pAgBCD, pBCgD);
  prob::condition_conditionals(pABCgD, pBCgD, qAgBCD);

  EXPECT_LT((qAgBCD-pAgBCD).array().abs().sum(), 1e-10);
}

TEST_F(Algebra, Bayes)
{
  prob::init::random(pAgBC, gen);
  prob::init::random(pBC, gen);

  auto pA = prob::uncondition(pAgBC, pBC).marginalize<0>();

  pBCgA = prob::bayes(pAgBC, pA, pBC);
  qAgBC = prob::bayes(pBCgA, pBC, pA);

  EXPECT_LT((qAgBC-pAgBC).array().abs().sum(), 1e-10);
}

