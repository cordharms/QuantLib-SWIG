; Copyright (C) 2002, 2003 RiskMap srl
;
; This file is part of QuantLib, a free-software/open-source library
; for financial quantitative analysts and developers - http://quantlib.org/
;
; QuantLib is free software developed by the QuantLib Group; you can
; redistribute it and/or modify it under the terms of the QuantLib License;
; either version 1.0, or (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; QuantLib License for more details.
;
; You should have received a copy of the QuantLib License along with this
; program; if not, please email ferdinando@ametrano.net
;
; The QuantLib License is also available at http://quantlib.org/license.html
; The members of the QuantLib Group are listed in the QuantLib License

(require (lib "quantlib.ss" "quantlib"))
(load "unittest.scm")
(load "utilities.scm")

(load "calendars.scm")
(load "capfloor.scm")
(load "covariance.scm")
(load "date.scm")
(load "daycounters.scm")
(load "distributions.scm")
(load "europeanoption.scm")
(load "instruments.scm")
(load "marketelements.scm")
(load "operators.scm")
(load "piecewiseflatforward.scm")
(load "riskstatistics.scm")
(load "segmentintegral.scm")
(load "simpleswap.scm")
(load "solvers1d.scm")
(load "statistics.scm")
(load "swaption.scm")
(load "termstructures.scm")
; to be removed
(load "old_pricers.scm")

(define tests
  (list
   Calendar-suite
   CapFloor-suite
   Covariance-suite
   Date-suite
   DayCounter-suite
   Distribution-suite
   EuropeanOption-suite
   Instrument-suite
   MarketElement-suite
   Operator-suite
   PiecewiseFlatForward-suite
   RiskStatistics-suite
   SegmentIntegral-suite
   SimpleSwap-suite
   Solver1D-suite
   Statistics-suite
   Swaption-suite
   TermStructure-suite
   ; to be removed
   OldPricers-suite))

(define suite
  (apply make-test-suite "QuantLib-MzScheme test suite" (reverse tests)))

(test/text-ui suite)
