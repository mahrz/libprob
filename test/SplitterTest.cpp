/*
 * SplitterTest.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: harder
 */

#include "gtest/gtest.h"
#include "TestVariables.hpp"

TEST(Splitter, SplitterRowColumnCount)
{
  int posteriors, conditionals;
  posteriors = prob::core::splitter<X,Y>::posteriors();
  EXPECT_EQ(posteriors, 2);
  conditionals = prob::core::splitter<X,Y>::conditionals();
  EXPECT_EQ(conditionals, 0);

  posteriors = prob::core::splitter<X,Y, prob::given, Z,W>::posteriors();
  EXPECT_EQ(posteriors, 2);
  conditionals = prob::core::splitter<X,Y, prob::given, Z,W>::conditionals();
  EXPECT_EQ(conditionals, 2);

  posteriors = prob::core::splitter<A,X,Y, prob::given, C,Z,W>::posteriors();
  EXPECT_EQ(posteriors, 3);
  conditionals = prob::core::splitter<B,X,Y, prob::given, C,Z,W>::conditionals();
  EXPECT_EQ(conditionals, 3);

  posteriors = prob::core::splitter<X, prob::given, Z>::posteriors();
  EXPECT_EQ(posteriors, 1);
  conditionals = prob::core::splitter<X, prob::given, Z>::conditionals();
  EXPECT_EQ(conditionals, 1);
}

TEST(Splitter, IndexCheckTest)
{
  bool is_valid;
  is_valid = prob::core::check_indices<3,4,0,0,1,2,3,4,5,6>::valid();
  EXPECT_TRUE(is_valid);

  is_valid = prob::core::check_indices<3,4,0,1,7>::valid();
  EXPECT_TRUE(!is_valid);

  is_valid = prob::core::check_indices<3,4,0,1,3,2>::valid();
  EXPECT_TRUE(!is_valid);
}

/** @todo Tests for other internal mechanics */
