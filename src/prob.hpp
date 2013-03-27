/*
 * This file is part of:
 * LibInf - Information Theory Library for the Eigen Matrix Library
 *
 * Copyright 2010 Malte Harder <m.harder@me.com>
 * Based on the libit information theory library by Sander van Dijk
 *
 * LibInf is released under the BSD License 
 *
 * TODO: Put license here
 *
 * LibInf is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef _PROB_H_
#define _PROB_H_

#include <eigen3/Eigen/Eigen>

/**
 * @brief Epsilon value to assume probabilities are actually 0
 */
#define PROB_EPSILON 1e-15

#include "Util/MakeIndices.hpp"
#include "Util/TupleFunctions.hpp"
#include "Util/TypeTraits.hpp"
#include "Util/Formatters.hpp"
	 
#include "RandomVariable.hpp"
#include "Splitter.hpp"
#include "Distribution.hpp"
	  
#include "Algebra.hpp"
#include "Initializers.hpp"
#include "InformationTheory.hpp"
#include "InformationTheory/Decomposition.hpp"

#endif /* _PROB_H_ */
