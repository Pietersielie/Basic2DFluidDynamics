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
