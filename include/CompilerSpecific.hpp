// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//  
//  Copyright © 2015-2020: Samuel Leweke¹²
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

// CASEMA_STRONG_INLINE uses __forceinline on MSVC and Intel ICPC, but doens't use always_inline of GCC.
#if (defined _MSC_VER) || (defined __INTEL_COMPILER)
	#define CASEMA_STRONG_INLINE __forceinline
#else
	#define CASEMA_STRONG_INLINE inline
#endif

// CASEMA_ALWAYS_INLINE makes the function inline by force, use with care
#ifdef __GNUC__
	#define CASEMA_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
	#define CASEMA_ALWAYS_INLINE CASEMA_STRONG_INLINE
#endif

// CASEMA_DONT_INLINE explicitly prevents inlining
#if (defined __GNUC__)
	#define CASEMA_DONT_INLINE __attribute__((noinline))
#elif (defined _MSC_VER)
	#define CASEMA_DONT_INLINE __declspec(noinline)
#else
	#define CASEMA_DONT_INLINE
#endif

// Assume release build as default
#if defined(NDEBUG) || !defined(DEBUG)
	#ifndef CASEMA_NO_DEBUG
		#define CASEMA_NO_DEBUG
	#endif
	#ifdef CASEMA_DEBUG
		#undef CASEMA_DEBUG
	#endif
#endif
#if defined(DEBUG) || defined(CASEMA_DEBUG)
	#define CASEMA_DEBUG
	#undef CASEMA_NO_DEBUG
	#include <cassert>
#endif

#ifdef CASEMA_NO_DEBUG
	// Disable assert in release
	#define casema_assert_impl(x)
#else
	// Use standard assert
	#define casema_assert_impl(x) assert(x)
#endif

// Allow user supplied casema_assert
#ifndef casema_assert
	#define casema_assert(x) casema_assert_impl(x)
#endif

#ifdef CASEMA_NO_DEBUG
	#define CASEMA_ONLY_USED_FOR_DEBUG(x) (void)x
#else
	#define CASEMA_ONLY_USED_FOR_DEBUG(x)
#endif

// Improve branch prediction by marking likely and unlikely execution paths
// Example:
//     if (some_unlikely_condition) { ... }
// transforms to
//     if (casema_unlikely(some_unlikely_condition)) { ... }
#ifdef __GNUC__
	#define casema_likely(x) __builtin_expect(!!(x), 1)
	#define casema_unlikely(x) __builtin_expect(!!(x), 0)
#endif

#ifndef casema_likely
	#define casema_likely(x) (x)
	#define casema_unlikely(x) (x)
#endif
