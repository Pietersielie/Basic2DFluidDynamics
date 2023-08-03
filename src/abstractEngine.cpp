/**
 * Basic 2D Fluid Dynamics - A basic wind tunnel simulator
 * Copyright (C) 2023  Cilliers Pretorius
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "abstractEngine.h"

namespace BasicFluidDynamics {
	namespace Engine {
		//Default constructor
        AbstractEngine::AbstractEngine(float t, bool active) : dt(t), execTimePassed(0.0), simTimePassed(0.0), isActive(active){}
		
		//One step through for the engine
		void AbstractEngine::step(){
			if (isActive){
				startComputation();
				step_internal();
				endComputation();
			}
		}
		
		void AbstractEngine::startComputation(){
			execStartPoint = std::chrono::system_clock::now();
		}
		
		void AbstractEngine::endComputation(){
			execTimePassed += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - execStartPoint).count()/1000000000.0;
			simTimePassed += dt;
		}

		float AbstractEngine::getRunTimePassed() const{
			return execTimePassed;
		}

		float AbstractEngine::getSimTimePassed() const{
			return simTimePassed;
		}
		
		const float AbstractEngine::getDt() const{
			return dt;
		}
	}
}
