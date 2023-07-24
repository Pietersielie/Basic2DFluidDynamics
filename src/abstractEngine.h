#ifndef PWM_ABSTRACT_ENGINE_H
#define PWM_ABSTRACT_ENGINE_H

#include <chrono>
namespace PWM{
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
			/*signals:
				void done();*/
		};
	}
}
#endif //PWM_ABSTRACT_ENGINE_H
