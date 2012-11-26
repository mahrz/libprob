#include "gtest/gtest.h"
#include "TestVariables.hpp"

class Distribution : public ::testing::Test {
 protected:
  virtual void SetUp()
  {
   for(A a = 0; a < A::extent();a++)
     for(B b = 0; b < B::extent();b++)
       for(C c = 0; c < C::extent();c++)
         for(D d = 0; d < D::extent();d++)
           pABgCD(a,b|c,d) = prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
             prob::read_index<C>::read(c) * prob::read_index<D>::read(d);

   for(A a = 0; a < A::extent();a++)
     for(B b = 0; b < B::extent();b++)
         pAB(a,b) = prob::read_index<A>::read(a) *B::extent() + prob::read_index<B>::read(b);
  }

  // virtual void TearDown() {}

 prob::distribution<double,A,B> pAB;
 prob::distribution<double,A,B, prob::given, C> pABgC;
 prob::distribution<double,A,B, prob::given, C,D> pABgCD;
};



TEST_F(Distribution, DistributionProperties)
{
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
TEST_F(Distribution, DistributionAssignment)
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
TEST_F(Distribution, DistributionReshape)
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
TEST_F(Distribution, ValueAccess)
{

  for(A a = 0; a < A::extent();a++)
      for(B b = 0; b < B::extent();b++)
        EXPECT_EQ(pAB(a,b), prob::read_index<A>::read(a) *B::extent() + prob::read_index<B>::read(b));

  for(A a = 0; a < A::extent();a++)
      for(B b = 0; b < B::extent();b++)
        EXPECT_EQ(pAB.prob_ref(a,b), prob::read_index<A>::read(a) *B::extent() + prob::read_index<B>::read(b));

  for(A a = 0; a < A::extent();a++)
      for(B b = 0; b < B::extent();b++)
        for(C c = 0; c < C::extent();c++)
          for(D d = 0; d < D::extent();d++)
            EXPECT_EQ(pABgCD(a,b|c,d), prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
              prob::read_index<C>::read(c) * prob::read_index<D>::read(d));

  for(A a = 0; a < A::extent();a++)
      for(B b = 0; b < B::extent();b++)
        for(C c = 0; c < C::extent();c++)
          for(D d = 0; d < D::extent();d++)
            EXPECT_EQ(pABgCD.prob_ref(a,b|c,d), prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
              prob::read_index<C>::read(c) * prob::read_index<D>::read(d));

  EXPECT_EQ(pAB.prob_read_or_zero(A(A::extent()),B(B::extent())), 0);
}


// Posterior access
TEST_F(Distribution, PosteriorAccess)
{
  for(A a = 0; a < A::extent();a++)
      for(B b = 0; b < B::extent();b++)
        for(C c = 0; c < C::extent();c++)
          for(D d = 0; d < D::extent();d++)
            EXPECT_EQ(pABgCD.posterior_distribution(c,d)(a,b), prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
              prob::read_index<C>::read(c) * prob::read_index<D>::read(d));

  for(A a = 0; a < A::extent();a++)
      for(B b = 0; b < B::extent();b++)
        for(C c = 0; c < C::extent();c++)
          for(D d = 0; d < D::extent();d++)
            EXPECT_EQ(pABgCD.posterior_distribution(c,d).prob_ref(a,b), prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
              prob::read_index<C>::read(c) * prob::read_index<D>::read(d));
}


TEST_F(Distribution, NormalizationAndSum)
{
  pABgCD.setConstant(1.0);

  for(A a = 0; a < A::extent();a++)
    for(B b = 0; b < B::extent();b++)
      for(C c = 0; c < C::extent();c++)
        for(D d = 0; d < D::extent();d++)
          EXPECT_EQ(pABgCD(a,b|c,d), 1.0);

  pABgCD.normalize();

  double n = 1.0 / (A::extent()*B::extent());

  for(A a = 0; a < A::extent();a++)
    for(B b = 0; b < B::extent();b++)
      for(C c = 0; c < C::extent();c++)
        for(D d = 0; d < D::extent();d++)
          EXPECT_EQ(pABgCD(a,b|c,d), n);

  pABgCD.setRandom();

  for(A a = 0; a < A::extent();a++)
    for(B b = 0; b < B::extent();b++)
      for(C c = 0; c < C::extent();c++)
        for(D d = 0; d < D::extent();d++)
          if(pABgCD(a,b|c,d) < 0)
            pABgCD(a,b|c,d) = -pABgCD(a,b|c,d);

  pABgCD.normalize();

  for(C c = 0; c < C::extent();c++)
    for(D d = 0; d < D::extent();d++)
      EXPECT_LT(abs(pABgCD.posterior_distribution(c,d).sum()-1), 1e-20);

  EXPECT_LT((pABgCD.sum_by_conditional() - Eigen::RowVectorXd::Ones(C::extent() * D::extent())).sum(), 1e-20);

  for(C c = 0; c < C::extent();c++)
    for(D d = 0; d < D::extent();d++)
      EXPECT_LT(abs(pABgCD.posterior_distribution(c,d).sum()-1), 1e-20);

}

// Map

TEST_F(Distribution, Map)
{
  pABgCD.map([] (double p) { return p*p; });

  pABgCD.each_index([&] (const A& a, const B& b, prob::given g, const C& c, const D& d)
      {
        double p = prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
                    prob::read_index<C>::read(c) * prob::read_index<D>::read(d);

        EXPECT_EQ(pABgCD(a,b|c,d), p*p );
      });

  prob::distribution<double,A,B, prob::given, C,D> qABgCD;

  qABgCD = pABgCD.map_copy([] (double p) { return sqrt(p); });

  qABgCD.each_index([&] (const A& a, const B& b, prob::given g, const C& c, const D& d)
      {
        double p = prob::read_index<A>::read(a) * prob::read_index<B>::read(b) *
                    prob::read_index<C>::read(c) * prob::read_index<D>::read(d);

        EXPECT_EQ(qABgCD(a,b|c,d), p );
      });
}


TEST_F(Distribution, MapConditional)
{
  pABgCD.setConstant(1.0);
  pABgCD.map_by_conditional([] (prob::distribution<double,A,B> &&d) { d(A(0),B(0)) = 0; d.normalize(); return d; });

  double p =  A::extent() * B::extent() - 1;

  pABgCD.each_index([&] (const A& a, const B& b, prob::given g, const C& c, const D& d)
      {
        if(a==0 && b==0)
          EXPECT_EQ(pABgCD(a,b|c,d), 0 );
        else
          EXPECT_LT(abs(pABgCD(a,b|c,d)-1/p), 1e-20);

      });

  prob::distribution<double,A,B, prob::given, C,D> qABgCD;

  pABgCD.setConstant(1.0);
  qABgCD = pABgCD.map_copy_by_conditional([] (prob::distribution<double,A,B> &&d) { d(A(0),B(0)) = 0; d.normalize(); return d; });

  qABgCD.each_index([&] (const A& a, const B& b, prob::given g, const C& c, const D& d)
      {
        EXPECT_EQ(pABgCD(a,b|c,d), 1 );

        if(a==0 && b==0)
          EXPECT_EQ(qABgCD(a,b|c,d), 0 );
        else
          EXPECT_LT(abs(qABgCD(a,b|c,d)-1/p), 1e-20);
      });
}

TEST_F(Distribution, GroupedMapSum)
{
  prob::distribution<double,X,Y,prob::given,Z,W> pXYgZW(X(2),Y(2)|Z(2),W(2));

  for(Z z = 0; z < 2; z++)
  {
    for(W w = 0; w < 2; w++)
    {
      pXYgZW(X(0),Y(0)|z,w) = 0.2;
      pXYgZW(X(1),Y(0)|z,w) = 0.3;
      pXYgZW(X(0),Y(1)|z,w) = 0.1;
      pXYgZW(X(1),Y(1)|z,w) = 0.4;
    }
  }

  auto pYgZW = pXYgZW.grouped_sum<1,2,3>();

  pYgZW.each_index([&pYgZW] (const Y& y, prob::given, const Z& z, const W& w)
  {
        EXPECT_EQ(pYgZW(y|z,w), 0.5);
  });

  auto pXgZW = pXYgZW.grouped_sum<0,2,3>();

  pXgZW.each_conditional_index([&pXgZW] (const Z& z, const W& w)
  {
        EXPECT_LT(abs(pXgZW(X(0)|z,w)-0.3), 1e-20);
        EXPECT_LT(abs(pXgZW(X(1)|z,w)-0.7), 1e-20);
  });

}

// Output / Input
TEST_F(Distribution, InputOutput)
{
  prob::distribution<double,X,Y,prob::given,Z,W> pXYgZW(X(2),Y(2)|Z(2),W(2));
  prob::distribution<double,X,Y,prob::given,Z,W> qXYgZW(X(2),Y(2)|Z(2),W(2));

  qXYgZW.setZero();

  pXYgZW(X(0), Y(0) | Z(0), W(0)) = 0.451703;
  pXYgZW(X(1), Y(0) | Z(0), W(0)) = 0.267738;
  pXYgZW(X(0), Y(1) | Z(0), W(0)) = 0.210834;
  pXYgZW(X(1), Y(1) | Z(0), W(0)) = 0.0697242;
  pXYgZW(X(0), Y(0) | Z(1), W(0)) = 0.392615;
  pXYgZW(X(1), Y(0) | Z(1), W(0)) = 0.170774;
  pXYgZW(X(0), Y(1) | Z(1), W(0)) = 0.338111;
  pXYgZW(X(1), Y(1) | Z(1), W(0)) = 0.0984997;
  pXYgZW(X(0), Y(0) | Z(0), W(1)) = 0.140682;
  pXYgZW(X(1), Y(0) | Z(0), W(1)) = 0.178517;
  pXYgZW(X(0), Y(1) | Z(0), W(1)) = 0.0761279;
  pXYgZW(X(1), Y(1) | Z(0), W(1)) = 0.604673;
  pXYgZW(X(0), Y(0) | Z(1), W(1)) = 0.297834;
  pXYgZW(X(1), Y(0) | Z(1), W(1)) = 0.379176;
  pXYgZW(X(0), Y(1) | Z(1), W(1)) = 0.0370919;
  pXYgZW(X(1), Y(1) | Z(1), W(1)) = 0.285898;

  std::string output = "Conditional Distribution\n\
X Y | Z W\n\
2 2 | 2 2\n\
0 0 | 0 0 : 0.451703\n\
1 0 | 0 0 : 0.267738\n\
0 1 | 0 0 : 0.210834\n\
1 1 | 0 0 : 0.0697242\n\
0 0 | 1 0 : 0.392615\n\
1 0 | 1 0 : 0.170774\n\
0 1 | 1 0 : 0.338111\n\
1 1 | 1 0 : 0.0984997\n\
0 0 | 0 1 : 0.140682\n\
1 0 | 0 1 : 0.178517\n\
0 1 | 0 1 : 0.0761279\n\
1 1 | 0 1 : 0.604673\n\
0 0 | 1 1 : 0.297834\n\
1 0 | 1 1 : 0.379176\n\
0 1 | 1 1 : 0.0370919\n\
1 1 | 1 1 : 0.285898\n";

  std::stringstream s;
  s << pXYgZW;
  EXPECT_EQ(output, s.str());

  qXYgZW = prob::distribution<double,X,Y,prob::given,Z,W>::load(s);

  EXPECT_EQ(qXYgZW, pXYgZW);
}


