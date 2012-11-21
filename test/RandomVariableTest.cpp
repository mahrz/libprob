/*
 * RandomVariableTest.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: harder
 */

#include "gtest/gtest.h"
#include "TestVariables.hpp"

class RandomVariable : public ::testing::Test {
 protected:
  virtual void SetUp()
  {
    x = 5; y = 6; z = 2;
    a = 5; b = 6; c = 2;
  }

  // virtual void TearDown() {}

  X x;
  Y y;
  Z z;

  A a;
  B b;
  C c;
};

TEST_F(RandomVariable, ConstructTest)
{
  // Dynamic
  EXPECT_EQ(x._val,5);

  // Static
  EXPECT_EQ(a._val,5);

  Y y(5);
  B b(5);

  // Dynamic
  EXPECT_EQ(y._val,5);

  // Static
  EXPECT_EQ(b._val,5);
}

TEST_F(RandomVariable, AssignTest)
{
  x = 5;
  EXPECT_EQ(x._val, 5);
  a = 5;
  EXPECT_EQ(a._val, 5);
}

TEST_F(RandomVariable, ArithmeticTest)
{

  EXPECT_EQ((x*y)._val, 30);
  EXPECT_EQ((y/z)._val, 3);
  EXPECT_EQ((x+y)._val, 11);
  EXPECT_EQ((y-z)._val, 4);
  EXPECT_EQ((x%z)._val, 1);

  EXPECT_EQ((a*b)._val, 30);
  EXPECT_EQ((b/c)._val, 3);
  EXPECT_EQ((a+b)._val, 11);
  EXPECT_EQ((b-c)._val, 4);
  EXPECT_EQ((a%c)._val, 1);
}

TEST_F(RandomVariable, ModificationArithmeticTest)
{
  EXPECT_EQ((x++)._val, 5);
  EXPECT_EQ((++x)._val, 7);

  EXPECT_EQ((a++)._val, 5);
  EXPECT_EQ((++a)._val, 7);

  EXPECT_EQ((x--)._val, 7);
  EXPECT_EQ((--x)._val, 5);

  EXPECT_EQ((a--)._val, 7);
  EXPECT_EQ((--a)._val, 5);

  x += y;
  x -= z;
  EXPECT_EQ(x._val, 9);

  a += b;
  a -= c;
  EXPECT_EQ(a._val, 9);

  x += 2;
  x -= 2;
  EXPECT_EQ(x._val, 9);

  a += 2;
  a -= 2;
  EXPECT_EQ(a._val, 9);

  x *= y;
  x /= z;
  EXPECT_EQ(x._val, 27);

  a *= b;
  a /= c;
  EXPECT_EQ(a._val, 27);

  x *= 2;
  x /= 2;
  EXPECT_EQ(x._val, 27);

  a *= 2;
  a /= 2;
  EXPECT_EQ(a._val, 27);
}

TEST_F(RandomVariable, CompareTest)
{
  EXPECT_TRUE(x < y);
  EXPECT_TRUE(x <= x);
  EXPECT_TRUE(x <= y);
  EXPECT_TRUE(x == x);
  EXPECT_TRUE(x != y);
  EXPECT_TRUE(x > z);
  EXPECT_TRUE(x >= x);
  EXPECT_TRUE(x >= z);

  EXPECT_TRUE(a < b);
  EXPECT_TRUE(a <= a);
  EXPECT_TRUE(a <= b);
  EXPECT_TRUE(a == a);
  EXPECT_TRUE(a != b);
  EXPECT_TRUE(a > c);
  EXPECT_TRUE(a >= a);
  EXPECT_TRUE(a >= c);
}

TEST_F(RandomVariable, CompareIntegerTest)
{
  EXPECT_TRUE(x < 6);
  EXPECT_TRUE(x <= 5);
  EXPECT_TRUE(x <= 6);
  EXPECT_TRUE(x == 5);
  EXPECT_TRUE(x != 6);
  EXPECT_TRUE(x > 2);
  EXPECT_TRUE(x >= 5);
  EXPECT_TRUE(x >= 2);

  EXPECT_TRUE(a < 6);
  EXPECT_TRUE(a <= 5);
  EXPECT_TRUE(a <= 6);
  EXPECT_TRUE(a == 5);
  EXPECT_TRUE(a != 6);
  EXPECT_TRUE(a > 2);
  EXPECT_TRUE(a >= 5);
  EXPECT_TRUE(a >= 2);
}

TEST_F(RandomVariable, ReadTest)
{
  EXPECT_EQ(prob::read_index<int>::read(5), 5);
  EXPECT_EQ(prob::read_index<X>::read(x), 5);
  EXPECT_EQ(prob::read_index<A>::read(a), 5);
}

TEST(StaticDynamic, IsStaticTest)
{
  bool is_static;

  is_static = prob::extents_assert<X,Y,Z>::is_static();
  EXPECT_TRUE(!is_static);

  is_static = prob::extents_assert<A,B,C>::is_static();
  EXPECT_TRUE(is_static);

  is_static = prob::extents_assert<X,Y,A>::is_static();
  EXPECT_TRUE(is_static);

  is_static = prob::extents_assert<X,prob::given,A>::is_static();
  EXPECT_TRUE(is_static);

  is_static = prob::extents_assert<A,prob::given,B>::is_static();
  EXPECT_TRUE(is_static);
}

TEST(StaticDynamic, StaticRowExtentsTest)
{
  int sizeA = prob::core::static_row_extents<A>::size();
  EXPECT_EQ(sizeA, 1);
  auto extentsA = prob::core::static_row_extents<A>::extents();
  EXPECT_EQ(extentsA, std::make_tuple<>());

  int sizeAgB = prob::core::static_row_extents<A, prob::given, B>::size();
  EXPECT_EQ(sizeAgB, B::extent());
  auto extentsAgB = prob::core::static_row_extents<A, prob::given, B>::extents();
  EXPECT_EQ(extentsAgB, std::make_tuple<int>(B::extent()));

  int sizeAgBCD = prob::core::static_row_extents<A, prob::given, B, C, D>::size();
  EXPECT_EQ(sizeAgBCD, B::extent()*C::extent()*D::extent());
  auto extentsAgBCD = prob::core::static_row_extents<A, prob::given, B, C, D>::extents();
  auto expected = std::make_tuple<int,int,int>(B::extent(), C::extent(), D::extent());
  EXPECT_EQ(extentsAgBCD, expected);
}

TEST(StaticDynamic, StaticColExtentsTest)
{
  int sizeA = prob::core::static_col_extents<A>::size();
  EXPECT_EQ(sizeA, A::extent());
  auto extentsA = prob::core::static_col_extents<A>::extents();
  EXPECT_EQ(extentsA, std::make_tuple<int>(A::extent()));

  int sizeAgB = prob::core::static_col_extents<A, prob::given, B>::size();
  EXPECT_EQ(sizeAgB, A::extent());
  auto extentsAgB = prob::core::static_col_extents<A, prob::given, B>::extents();
  EXPECT_EQ(extentsAgB, std::make_tuple<int>(A::extent()));

  int sizeABCDgE = prob::core::static_col_extents<A, B, C, D, prob::given, E>::size();
  EXPECT_EQ(sizeABCDgE, A::extent()*B::extent()*C::extent()*D::extent());
  std::tuple<int,int,int,int> extentsABCDgE = prob::core::static_col_extents<A, B, C, D, prob::given, E>::extents();
  auto expected = std::make_tuple<int,int,int,int>(A::extent(), B::extent(), C::extent(), D::extent());
  EXPECT_EQ(extentsABCDgE, expected);
}

TEST(Iterator, IndexIteratorTest)
{
  auto extents = std::make_tuple(A::extent(), B::extent(), 1, C::extent(), D::extent());

  int _a = 0;
  int _b = 0;
  int _c = 0;
  int _d = 0;

  auto f = [&_a, &_b, &_c, &_d] (const A& a, const B& b, const prob::given& g, const C& c, const D& d)
  {
    EXPECT_EQ(prob::read_index<A>::read(a), _a);
    EXPECT_EQ(prob::read_index<B>::read(b), _b);
    EXPECT_EQ(prob::read_index<C>::read(c), _c);
    EXPECT_EQ(prob::read_index<D>::read(d), _d);

    _d++;
    if(_d == D::extent())
    {
      _d = 0;
      _c++;
      if(_c == C::extent())
      {
        _c = 0;
        _b++;
        if(_b == B::extent())
        {
          _b = 0;
          _a++;
          if(_a == A::extent())
          {
            _a = 0;
          }
        }
      }
    }
  };

  prob::core::index_iterator< prob::core::vars<A, B, prob::given, C, D> >::apply_all(f,
              std::make_tuple(), extents);

}

TEST(Iterator, ReverseIndexIteratorTest)
{
  auto extents = std::make_tuple(A::extent(), B::extent(), 1, C::extent(), D::extent());

  int _a = 0;
  int _b = 0;
  int _c = 0;
  int _d = 0;

  auto f = [&_a, &_b, &_c, &_d] (const A& a, const B& b, const prob::given& g, const C& c, const D& d)
  {
    EXPECT_EQ(prob::read_index<A>::read(a), _a);
    EXPECT_EQ(prob::read_index<B>::read(b), _b);
    EXPECT_EQ(prob::read_index<C>::read(c), _c);
    EXPECT_EQ(prob::read_index<D>::read(d), _d);

    _a++;
    if(_a == A::extent())
    {
      _a = 0;
      _b++;
      if(_b == B::extent())
      {
        _b = 0;
        _c++;
        if(_c == C::extent())
        {
          _c = 0;
          _d++;
          if(_d == D::extent())
          {
            _d = 0;
          }
        }
      }
    }
  };

  prob::core::index_iterator_reverse< prob::core::vars<A, B, prob::given, C, D> >::apply_all(f,
              std::make_tuple(), extents);
}
