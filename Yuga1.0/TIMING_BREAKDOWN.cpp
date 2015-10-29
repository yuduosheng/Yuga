/******************************************************************************
 * This file is part of Cubica++.
 
 * Cubica++ is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.

 * Cubica++ is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 * details.

 * You should have received a copy of the GNU General Public License along
 * with Cubica++. If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/
#include "TIMING_BREAKDOWN.h"

map<string, double> TIMING_BREAKDOWN::_timingBreakdown;
vector<Real> TIMING_BREAKDOWN::_frameTime;
TIMER TIMING_BREAKDOWN::_totalTimer;
TIMER TIMING_BREAKDOWN::_currentTimer;

Real TIMING_BREAKDOWN::_totalTime = 0;
int TIMING_BREAKDOWN::_totalSteps = 0;
// to detect if tic is called multiple times before toc
int TIMING_BREAKDOWN::_ticCnt = 0;
