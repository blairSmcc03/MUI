/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file temporal_sampler_mean.h
 * @author B McCormick
 * @date 19 November 2025
 * @brief Temporal sampler that interpolates between the closest two dataframes
 */

#ifndef MUI_TEMPORAL_SAMPLER_LINEAR_H_
#define MUI_TEMPORAL_SAMPLER_LINEAR_H_

#include "../../general/util.h"
#include "../../mui_config.h"

namespace mui {

template<typename CONFIG=default_config> class temporal_sampler_linear {
public:
	using REAL       	= typename CONFIG::REAL;
	using INT        	= typename CONFIG::INT;
	using time_type  	= typename CONFIG::time_type;
	using iterator_type = typename CONFIG::iterator_type;
	
	temporal_sampler_linear(time_type dtNeighbour=time_type(0)) {
		dtNeighbour_ = dtNeighbour;
	}

	//- Filter based on time input
	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<std::pair<time_type,iterator_type>, TYPE> > &points ) const {
		// find closest dataframe before focus and after focus
		std::pair<std::pair<time_type,iterator_type>, TYPE> previousFrame, followingFrame;
		time_type  minimumFoundBefore = std::numeric_limits<double>::infinity();
		time_type maximumFoundAfter = -1* std::numeric_limits<double>::infinity();

		bool valueBeforeFound, valueAfterFound = false;
		for( auto i: points ) {

			time_type dt = focus - i.first.first;
			if(dt >= 0 && dt < minimumFoundBefore){
				minimumFoundBefore = dt;
				previousFrame = i;
				valueBeforeFound = true;
			}
			
			if(dt <= 0 && dt > maximumFoundAfter){
				maximumFoundAfter = dt;
				followingFrame = i; 
				valueAfterFound = true;
			}
		}

		

		if ( valueBeforeFound && valueAfterFound ){
			// Linear interpolation
		
			if(followingFrame.first.first - previousFrame.first.first > 0){
				TYPE gradient = (followingFrame.second - previousFrame.second)/(followingFrame.first.first - previousFrame.first.first);
				return previousFrame.second + (focus-previousFrame.first.first)*gradient;
			}else{
				return previousFrame.second;
			}
		}
		else {
			return TYPE(0);
		}
	}

	//- Filter based on time and iterator input - only time used
	template<typename TYPE>
	TYPE filter( std::pair<time_type,iterator_type> focus, const std::vector<std::pair<std::pair<time_type,iterator_type>, TYPE> > &points ) const {
		// find closest dataframe before focus and after focus
		std::pair<std::pair<time_type,iterator_type>, TYPE> previousFrame, followingFrame;
		time_type  minimumFoundBefore = std::numeric_limits<double>::infinity();
		time_type maximumFoundAfter = -1* std::numeric_limits<double>::infinity();

		bool valueBeforeFound, valueAfterFound = false;
		for( auto i: points ) {

			time_type dt = focus.first - i.first.first;
			if(dt >= 0 && dt < minimumFoundBefore){
				minimumFoundBefore = dt;
				previousFrame = i;
				valueBeforeFound = true;
			}
			
			if(dt <= 0 && dt > maximumFoundAfter){
				maximumFoundAfter = dt;
				followingFrame = i; 
				valueAfterFound = true;
			}
		}

		

		if ( valueBeforeFound && valueAfterFound ){
			// Linear interpolation
		
			if(followingFrame.first.first - previousFrame.first.first > 0){
				TYPE gradient = (followingFrame.second - previousFrame.second)/(followingFrame.first.first - previousFrame.first.first);
				return previousFrame.second + (focus.first-previousFrame.first.first)*gradient;
			}else{
				return previousFrame.second;
			}
		}
		else {
			return TYPE(0);
		}
	}

	time_type get_upper_bound( time_type focus ) const {
		if(dtNeighbour_ > 0){
			return std::ceil(focus/dtNeighbour_)*dtNeighbour_;
		}
		else{
			return std::numeric_limits<time_type>::max();
		}

	}

	time_type get_lower_bound( time_type focus ) const {
		if(dtNeighbour_ > 0){
			return std::floor(focus/dtNeighbour_)*dtNeighbour_;
		}
		else{
			return std::numeric_limits<time_type>::lowest();
		}
	}

	time_type get_barrier_time(time_type focus) const {
		return focus;
	}

protected:
	time_type dtNeighbour_;
};

}

#endif /* MUI_TEMPORAL_SAMPLER_SUM_H_ */
