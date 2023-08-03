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

#ifndef BFD_ABSTRACT_ENGINE_H
#define BFD_ABSTRACT_ENGINE_H

#include <chrono>
namespace BasicFluidDynamics{
	namespace Engine{
		class AbstractEngine{
			private:
                //The timestep for this engine in seconds
				float dt;
			
				//Monitoring utilities as used in MWM
				std::chrono::system_clock::time_point execStartPoint;
				float execTimePassed, simTimePassed;
			protected:
				//This function must be overrode by any sub class to provide the actual computation.
                void startComputation();
                void endComputation();
                virtual void step_internal() = 0;

			public:
				//Bool variable to say if this engine is active (e.g., user interaction will start off false until triggered).
				bool isActive;
			
				//Default constructor
				AbstractEngine(float dt = 1.0, bool active = true);

				//Call to perform one step of the engine
				void step();
				//virtual void clear();
				float getRunTimePassed() const;
				float getSimTimePassed() const;
                const float getDt() const;
		};
	}
}
#endif //BFD_ABSTRACT_ENGINE_H
