#include "gtest/gtest.h"
#include "prob"

RVAR_STATIC(X,4)
RVAR_STATIC(Y,4)
RVAR_STATIC(Z,2)


class Information : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    gen = std::mt19937(rd());

    pXY(X(0),Y(0)) = 1/8.0;
    pXY(X(0),Y(1)) = 1/16.0;
    pXY(X(0),Y(2)) = 1/16.0;
    pXY(X(0),Y(3)) = 1/4.0;

    pXY(X(1),Y(0)) = 1/16.0;
    pXY(X(1),Y(1)) = 1/8.0;
    pXY(X(1),Y(2)) = 1/16.0;
    pXY(X(1),Y(3)) = 0;

    pXY(X(2),Y(0)) = 1/32.0;
    pXY(X(2),Y(1)) = 1/32.0;
    pXY(X(2),Y(2)) = 1/16.0;
    pXY(X(2),Y(3)) = 0;

    pXY(X(3),Y(0)) = 1/32.0;
    pXY(X(3),Y(1)) = 1/32.0;
    pXY(X(3),Y(2)) = 1/16.0;
    pXY(X(3),Y(3)) = 0;

    pX = pXY.marginalize<0>();
    pY = pXY.marginalize<1>();

    prob::condition(pXY, pY, pXgY);
    prob::condition(pXY.marginalize<1,0>(), pX, pYgX);

    p(Z(0)) = 1/2.0;
    p(Z(1)) = 1/2.0;

    q(Z(0)) = 3/4.0;
    q(Z(1)) = 1/4.0;
  }

  std::random_device rd;
  std::mt19937 gen;

  prob::distribution<double,X,Y> pXY;
  prob::distribution<double,X, prob::given, Y> pXgY;
  prob::distribution<double,X> pX;
  prob::distribution<double,Y, prob::given, X> pYgX;
  prob::distribution<double,Y> pY;

  prob::distribution<double,Z> p;
  prob::distribution<double,Z> q;
};

TEST_F(Information, Entropy)
{
  EXPECT_LT(abs(prob::it::entropy(pXY) - 3.375), 1e-10);
  EXPECT_LT(abs(prob::it::entropy(pX) - 1.75), 1e-10);
  EXPECT_LT(abs(prob::it::entropy(pY) - 2.0), 1e-10);

  EXPECT_LT(abs(prob::it::conditional_entropy(pYgX, pX) - 1.625), 1e-10);
  EXPECT_LT(abs(prob::it::conditional_entropy(pXgY, pY) - 1.625), 1e-10);
}

TEST_F(Information, MutualInformation)
{
  EXPECT_LT(abs(prob::it::mutual_information(pYgX, pX) - 0.375), 1e-10);
  EXPECT_LT(abs(prob::it::mutual_information(pXgY, pY) - 0.375), 1e-10);
}

TEST_F(Information, KullbackLeiblerDivergence)
{
  EXPECT_LT(abs(prob::it::kl_divergence(p, q) - 0.2075), 1e-10);
  EXPECT_LT(abs(prob::it::kl_divergence(q, p) - 0.1887), 1e-10);
}

TEST_F(Information, JensenShannonDivergence)
{
  EXPECT_LT(abs(prob::it::js_divergence(p, q, 0.5) - 0.048795), 1e-10);
  EXPECT_LT(abs(prob::it::js_divergence(q, p, 0.5) - 0.048795), 1e-10);
}

