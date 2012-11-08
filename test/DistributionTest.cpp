#include "gtest/gtest.h"
#include "TestVariables.hpp"

TEST(Distribution, DistributionProperties)
{
  prob::distribution<double,A,B> pAB;
  prob::distribution<double,A,B, prob::given, C> pABgC;

  prob::distribution<double,X> pX(X(10));
  prob::distribution<double,X,prob::given,Y,Z> pXgYZ(X(10)|Y(5),Z(7));

  EXPECT_FALSE((prob::distribution<double,A,B>::conditional_distribution()));
  EXPECT_FALSE((prob::distribution<double,X,Y>::conditional_distribution()));

  EXPECT_TRUE((prob::distribution<double,A,B, prob::given, C>::conditional_distribution()));
  EXPECT_TRUE((prob::distribution<double,X, prob::given, Y, Z>::conditional_distribution()));

  auto req_ab_extent = std::make_tuple(A(A::extent()), B(B::extent()));
  auto req_c_extent = std::make_tuple(C(C::extent()));

  auto req_x_extent = std::make_tuple(X(10));
  auto req_yz_extent = std::make_tuple(Y(5), Z(7));

  EXPECT_EQ(pAB.posterior_extents(), req_ab_extent);
  EXPECT_EQ(pABgC.posterior_extents(), req_ab_extent);
  EXPECT_EQ(pABgC.conditional_extents(), req_c_extent);
  EXPECT_EQ(pAB.col_extents(), req_ab_extent);
  EXPECT_EQ(pABgC.col_extents(), req_ab_extent);
  EXPECT_EQ(pABgC.row_extents(), req_c_extent);

  EXPECT_EQ(pX.posterior_extents(), req_x_extent);
  EXPECT_EQ(pXgYZ.posterior_extents(), req_x_extent);
  EXPECT_EQ(pXgYZ.conditional_extents(), req_yz_extent);
  EXPECT_EQ(pX.col_extents(), req_x_extent);
  EXPECT_EQ(pXgYZ.col_extents(), req_x_extent);
  EXPECT_EQ(pXgYZ.row_extents(), req_yz_extent);
}

// Distribution Assignment
TEST(Distribution, DistributionAssignment)
{
  Eigen::MatrixXd mAB;
  Eigen::MatrixXd mABgC;

  prob::distribution<double,A,B> pAB;
  prob::distribution<double,A,B, prob::given, C> pABgC;

  prob::distribution<double,A,B> qAB;
  prob::distribution<double,A,B, prob::given, C> qABgC;

  std::random_device rd;
  std::mt19937 gen(rd());

  prob::init::random(pAB,gen);
  prob::init::random(pABgC,gen);

  mAB = pAB;
  mABgC = pABgC;

  qAB = pAB;
  qABgC = pABgC;

  EXPECT_EQ(qAB, pAB);
  EXPECT_EQ(qABgC, pABgC);

  EXPECT_EQ(qAB.posterior_extents(), pAB.posterior_extents());
  EXPECT_EQ(qAB.conditional_extents(), pAB.conditional_extents());
  EXPECT_EQ(qABgC.posterior_extents(), pABgC.posterior_extents());
  EXPECT_EQ(qABgC.conditional_extents(), pABgC.conditional_extents());

  prob::distribution<double,A,B> rAB(mAB, pAB.row_extents(), pAB.col_extents());
  prob::distribution<double,A,B, prob::given, C> rABgC(mABgC, pABgC.row_extents(), pABgC.col_extents());

  EXPECT_EQ(pAB, rAB);
  EXPECT_EQ(pABgC, rABgC);

  qAB.each_index([&rAB, &qAB] (const A& a, const B& b) { EXPECT_EQ(qAB(a,b),rAB(a,b)); });
  qABgC.each_index([&rABgC, &qABgC] (const A& a, const B& b, prob::given g, const C& c) { EXPECT_EQ(qABgC(a,b,g,c),rABgC(a,b,g,c)); });
}

// Reshape
TEST(Distribution, DistributionReshape)
{
  prob::distribution<double,X> pX(X(10));
  prob::distribution<double,X,prob::given,Y,Z> pXgYZ(X(10)|Y(5),Z(7));

  pX.reshape_dimensions(std::make_tuple(),std::make_tuple(X(5)));
  pXgYZ.reshape_dimensions(std::make_tuple(Y(2),Z(4)), std::make_tuple(X(5)));

  EXPECT_EQ(pX.posterior_extents(), std::make_tuple(X(5)));
  EXPECT_EQ(pXgYZ.posterior_extents(), std::make_tuple(X(5)));
  EXPECT_EQ(pXgYZ.conditional_extents(), std::make_tuple(Y(2),Z(4)));

  std::random_device rd;
  std::mt19937 gen(rd());
  prob::init::random(pX,gen);
  prob::init::random(pXgYZ,gen);

  prob::distribution<double,X> qX(pX);
  prob::distribution<double,X,prob::given,Y,Z> qXgYZ(pXgYZ);

  qX.reshape(X(10));
  qXgYZ.reshape(X(10)|Y(5),Z(2));

  pX.each_index([&pX,&qX] (const X& x) { EXPECT_LT(abs(pX(x)-qX(x)), 1e-10); } );

}

// Value access

// Posterior access

// Index iteration

// Normalization

// Sum

// Sum by conditional

// Map

// Map copy

// Map by conditional

// Map copy by conditional

// Grouped map sum

// Grouped sum

// Output / Input
