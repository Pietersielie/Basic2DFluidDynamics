#ifndef PWM_MODEL_AIR_LAYER_H
#define PWM_MODEL_AIR_LAYER_H

#include <algorithm>
#include "atmofuncs.h"
#include <cstddef>
#include <iostream>
#include <memory>
#include "planet.h"
//#include "SSGDataStructure.h"
#include <typeinfo>
#include <vector>

namespace PWM{
	namespace Model{
		template<typename T, typename V> class airLayer{
			private:
				//The planet this air layer is on
				std::shared_ptr<planet> Planet;
				
				// The altitude of this layer in meters
				V height;
			
				// The thickness of the layer in meters, centered on the altitude
				V thickness;
			
				// Storing all the various bits and bobs that make up an air layer
				std::shared_ptr<T> obstacles;//obstacles, i.e. land
				std::shared_ptr<T> velocityTheta;//velocity in the north-south axis (meters per second)
				std::shared_ptr<T> velocityPhi;//velocity in the east-west axis (meters per second)
				std::shared_ptr<T> temperature;//temperature in Kelvin
				std::shared_ptr<T> moisture;//the absolute humidity for the layer in kg per cubic metre
				std::shared_ptr<T> condensedWater;//the condensed water, i.e., clouds in this layer in kg per cubic metre
				std::shared_ptr<T> pressure;//the dynamic pressure in Pascal
                std::shared_ptr<T> particulates;//The amount of particulates (ash/dirt/dust) for the layer in kg per cubic metre
                std::shared_ptr<T> cloudTypes;//the type of clouds at each layer

				/* //Storing old versions of the info used in advection with swap.
				std::shared_ptr<T> oldVelTheta;
				std::shared_ptr<T> oldVelPhi;
				std::shared_ptr<T> oldTemp;
				std::shared_ptr<T> oldMoisture;
				std::shared_ptr<T> oldClouds;
				std::shared_ptr<T> oldPressure; */

                V tropopause = 10000;
			public:
				std::vector<std::shared_ptr<T>> scalarQuantities;
				std::vector<std::shared_ptr<T>> oldScalarQuantities;
				
				//Various constructors for different parameter inputs.
				airLayer(std::shared_ptr<planet> P, V layerHeight = 0, V layerThickness = 0);

				//Constructor meant to be used with square2DArray
				airLayer(std::shared_ptr<planet> P, V layerHeight, V layerThickness, size_t faceWidth, const V size = 50000.f);

                //Constructor meant to be used with flatStaggeredGrid
                airLayer(std::shared_ptr<planet> P, V layerHeight, V layerThickness, const size_t width, const size_t height, const V xSize, const V ySize);

				// Initialise the speed and pressure values per a random generator
				void randomInit();
				void init();

				/*// Release the GPU resources
				void releaseResources();*/

				//Function to swap all the scalar values with the new values after advection
				void swapScalars(std::vector<std::shared_ptr<T>>& buffer);

				//Function to swap velocities with new values after parallel modifications.
				void swapVels(std::shared_ptr<T>& oldVelTheta, std::shared_ptr<T>& oldVelPhi);

				void swapPres(std::shared_ptr<T>& oldPressure);

				void swapMoistures(std::shared_ptr<T>& oldMoisture);

				void swapClouds(std::shared_ptr<T>& oldClouds);

				void swapTemps(std::shared_ptr<T>& oldTemp);

                void swapParts(std::shared_ptr<T>& oldParticulates);

				//Function to solve the polar velocities as per Yang and Corse.
				//Defined here since both convection and advection change velocities.
				void solvePolarVelocities();
			
				// Uniformly set the temperature of a layer
				void setTemperature(V heat);
			
				// Uniformly set the moisture of a layer, with random variation if the <rnd> flag is set. For now, <rnd> flag removed.
				void setMoisture(V vap);
			
				// Compute the total moisture in the layer
				V getCumulativeMoisture() const;
			
				// Compute the average temperature in the layer
				V getMeanTemperature() const;
			
				/*// Compute the average velocity amplitude
				double getMeanVelocityAmplitude() const;*/
			
				// Const getters for each of the data fields
				const V getRadius() const;
				const std::shared_ptr<planet> getPlanet() const;
				const V& getHeight() const;
                const V getLayerTop() const;
                const V getLayerBot() const;
				const V& getThickness() const;
				
				T& getObstacles();
                std::shared_ptr<T>& getObsPtr();
				T& getVelocityTheta();
				T& getVelocityPhi();
				T& getTemperature();
				T& getMoisture();
				T& getCondensedWater();
				T& getPressure();
                T& getParticulates();
                T& getClouds();
				
				const V getObstacles(const size_t index) const;
				const V getObstacles(const size_t i, const size_t j) const;
				const V getObstacles(const size_t face, const size_t i, const size_t j) const;
				const V getVelocityTheta(const size_t index) const;
				const V getVelocityTheta(const size_t i, const size_t j) const;
				const V getVelocityTheta(const size_t face, const size_t i, const size_t j) const;
				const V getVelocityPhi(const size_t index) const;
				const V getVelocityPhi(const size_t i, const size_t j) const;
				const V getVelocityPhi(const size_t face, const size_t i, const size_t j) const;
				const V getTemperature(const size_t index) const;
				const V getTemperature(const size_t i, const size_t j) const;
				const V getTemperature(const size_t face, const size_t i, const size_t j) const;
				const V getMoisture(const size_t index) const;
				const V getMoisture(const size_t i, const size_t j) const;
                const V getMoisture(const size_t face, const size_t i, const size_t j) const;
                const V getCondensedWater(const size_t index) const;
                const V getCondensedWater(const size_t i, const size_t j) const;
                const V getCondensedWater(const size_t face, const size_t i, const size_t j) const;
                const V getClouds(const size_t index) const;
                const V getClouds(const size_t i, const size_t j) const;
                const V getClouds(const size_t face, const size_t i, const size_t j) const;
                const V getPressure(const size_t index) const;
                const V getPressure(const size_t i, const size_t j) const;
                const V getPressure(const size_t face, const size_t i, const size_t j) const;
                const V getParticulates(const size_t index) const;
                const V getParticulates(const size_t i, const size_t j) const;
                const V getParticulates(const size_t face, const size_t i, const size_t j) const;

				// Setters for each of the data fields
				void setHeight(V h);
				void setThickness(V t);
				void setObstacles(T& obs);
				void setVelocityTheta(T& vTheta);
				void setVelocityPhi(T& vPhi);
				void setTemperature(T& temps);
				void setMoisture(T& moist);
				void setCondensedWater(T& clouds);
				void setPressure(T& pres);
				void setParticulates(T& parts);
                void setClouds(T& cloudMask);

				void setObstacles(const size_t index, const V& val);
				void setObstacles(const size_t i, const size_t j, const V& val);
				void setObstacles(const size_t face, const size_t i, const size_t j, const V& val);
				void setVelocityTheta(const size_t index, const V& val);
				void setVelocityTheta(const size_t i, const size_t j, const V& val);
				void setVelocityTheta(const size_t face, const size_t i, const size_t j, const V& val);
				void setVelocityPhi(const size_t index, const V& val);
				void setVelocityPhi(const size_t i, const size_t j, const V& val);
				void setVelocityPhi(const size_t face, const size_t i, const size_t j, const V& val);
				void setTemperature(const size_t index, const V& val);
				void setTemperature(const size_t i, const size_t j, const V& val);
				void setTemperature(const size_t face, const size_t i, const size_t j, const V& val);
				void setMoisture(const size_t index, const V& val);
				void setMoisture(const size_t i, const size_t j, const V& val);
                void setMoisture(const size_t face, const size_t i, const size_t j, const V& val);
                void setCondensedWater(const size_t index, const V& val);
                void setCondensedWater(const size_t i, const size_t j, const V& val);
                void setCondensedWater(const size_t face, const size_t i, const size_t j, const V& val);
                void setClouds(const size_t index, const V& val);
                void setClouds(const size_t i, const size_t j, const V& val);
                void setClouds(const size_t face, const size_t i, const size_t j, const V& val);
				void setPressure(const size_t index, const V& val);
				void setPressure(const size_t i, const size_t j, const V& val);
				void setPressure(const size_t face, const size_t i, const size_t j, const V& val);
				void setParticulates(const size_t index, const V& val);
				void setParticulates(const size_t i, const size_t j, const V& val);
				void setParticulates(const size_t face, const size_t i, const size_t j, const V& val);

				// Various comparators for an air layer
				bool operator==(airLayer& other);
				bool operator!=(airLayer& other);

				template<typename X, typename Y>
                friend bool operator<(const airLayer<X, Y>& lhs, const airLayer<X, Y>& rhs);
                template<typename X, typename Y>
                friend bool operator<(std::shared_ptr<airLayer<X, Y>> lhs, std::shared_ptr<airLayer<X, Y>> rhs);
        };
		
		/**int i = 0; i < quantity.size(); ++i
		 * Implementing the templated functions above on general rules... fingers crossed.
		 */
		template<typename T, typename V>
        inline airLayer<T, V>::airLayer(std::shared_ptr<planet> P, V layerHeight, V layerThickness) : Planet(P), height(layerHeight), thickness(layerThickness){
			obstacles = std::make_shared<T>();
			velocityTheta = std::make_shared<T>();
			velocityPhi = std::make_shared<T>();
			temperature = std::make_shared<T>();
			moisture = std::make_shared<T>();
			condensedWater = std::make_shared<T>();
			pressure = std::make_shared<T>();
			particulates = std::make_shared<T>();
            cloudTypes = std::make_shared<T>();

			scalarQuantities.push_back(temperature);
			scalarQuantities.push_back(moisture);
			scalarQuantities.push_back(condensedWater);
			scalarQuantities.push_back(pressure);
			scalarQuantities.push_back(particulates);
			
			//default temp set based on height, with higher meaning cooler in the troposphere. Even if weird stuff happens higher, we don't care.
			V tempLess = height / 100.0;
			setTemperature(getPlanet()->getAveTemp() - tempLess);
		}
		
		/**
         * This constructor is specifically meant for using PWM::PWMDataStructure::square2DArray<V> for the various fields.
		 */
		template<typename T, typename V>
        inline airLayer<T, V>::airLayer(std::shared_ptr<planet> P, V layerHeight, V layerThickness, size_t faceWidth, const V size) : Planet(P), height(layerHeight), thickness(layerThickness){
			obstacles = std::make_shared<T>(faceWidth, size);
			velocityTheta = std::make_shared<T>(faceWidth, size);
			velocityPhi = std::make_shared<T>(faceWidth, size);
			temperature = std::make_shared<T>(faceWidth, size);
			moisture = std::make_shared<T>(faceWidth, size);
			condensedWater = std::make_shared<T>(faceWidth, size);
            pressure = std::make_shared<T>(faceWidth, size);
            particulates = std::make_shared<T>(faceWidth, size);
            cloudTypes = std::make_shared<T>(faceWidth, size);

			scalarQuantities.push_back(temperature);
			scalarQuantities.push_back(moisture);
			scalarQuantities.push_back(condensedWater);
			scalarQuantities.push_back(pressure);
			scalarQuantities.push_back(particulates);
			
			//default temp set based on height, with higher meaning cooler in the troposphere. Even if weird stuff happens higher, we don't care.
			//Lapse rate is equivalent to roughly 1 degree C per 100 m altitude in dry air, good enough for our initialisation.
            V tempLess = height * 0.01;
            if (height > tropopause){//See cloudUtils.h for tropopause value
                tempLess = 100 + (height - tropopause) * 0.002;
            }
            setTemperature(getPlanet()->getAveTemp() - tempLess);

			init();
		}

        /**
         * This constructor is specifically meant for using PWM::PWMDataStructure::flatStaggeredGrid<V> for the various fields.
         */
        template<typename T, typename V>
        inline airLayer<T, V>::airLayer(std::shared_ptr<planet> P, V layerHeight, V layerThickness, const size_t width, const size_t height, const V xSize, const V ySize) : Planet(P), height(layerHeight), thickness(layerThickness){
            obstacles = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            velocityTheta = std::make_shared<T>(width, height, 0.0, 0.5, xSize, ySize);
            velocityPhi = std::make_shared<T>(width, height, 0.5, 0.0, xSize, ySize);
            temperature = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            moisture = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            condensedWater = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            pressure = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            particulates = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            cloudTypes = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);

            scalarQuantities.push_back(temperature);
            scalarQuantities.push_back(moisture);
            scalarQuantities.push_back(condensedWater);
            scalarQuantities.push_back(pressure);
            scalarQuantities.push_back(particulates);

            //default temp set based on height, with higher meaning cooler in the troposphere. Even if weird stuff happens higher, we don't care.
            //Lapse rate is equivalent to roughly 1 degree C per 120 m altitude in dry air, good enough for our initialisation.
            V tempLess = layerHeight * 0.0085;
            if (layerHeight > tropopause){//See cloudUtils.h for tropopause value
                tempLess = tropopause * 0.0085 + (layerHeight - tropopause) * 0.002;
            }
            setTemperature(getPlanet()->getAveTemp() - tempLess + 10);
//            std::cout << "Layer at height " << layerHeight << " has average temperature of " << getPlanet()->getAveTemp() - tempLess + 10 << " degrees C." << std::endl;

            init();
        }

		template<typename T, typename V>
        inline void airLayer<T, V>::init(){
			*obstacles = 0;
			*velocityTheta = 0;
			*velocityPhi = 0;
			*moisture = 0;
			*condensedWater = 0;
			*pressure = 0;
			*particulates = 0;
		}
		
		template<typename T, typename V>
        inline void airLayer<T, V>::setTemperature(V heat){
			*temperature = heat;
		}
		
		template<typename T, typename V>
        inline void airLayer<T, V>::setMoisture(V vap){
			*moisture = vap;
		}
		
		template<typename T, typename V>
    inline V airLayer<T, V>::getCumulativeMoisture() const{
			return moisture->sum();
		}
		
		template<typename T, typename V>
    inline V airLayer<T, V>::getMeanTemperature() const{
			return temperature->mean();
		}

		template<typename T, typename V>
    inline const std::shared_ptr<planet> airLayer<T, V>::getPlanet() const{
			return Planet;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getRadius() const{
			return Planet->getRadius() + height;
		}

        template<typename T, typename V>
    inline const V& airLayer<T, V>::getHeight() const{
            return height;
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getLayerTop() const{
            return getHeight() + (getThickness() / 2);
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getLayerBot() const{
            return getHeight() - (getThickness() / 2);
        }
		
		template<typename T, typename V>
    inline const V& airLayer<T, V>::getThickness() const{
			return thickness;
		}
		
		template<typename T, typename V>
    inline T& airLayer<T, V>::getObstacles(){
            return *obstacles;
        }

        template<typename T, typename V>
    inline std::shared_ptr<T>& airLayer<T, V>::getObsPtr(){
            return obstacles;
        }

		template<typename T, typename V>
    inline const V airLayer<T, V>::getObstacles(const size_t index) const{
			return (*obstacles).getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getObstacles(const size_t i, const size_t j) const{
			return (*obstacles).getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getObstacles(const size_t face, const size_t i, const size_t j) const{
			return obstacles->getData(face, i, j);
		}

		template<typename T, typename V>
    inline T& airLayer<T, V>::getVelocityTheta(){
            return *velocityTheta;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getVelocityTheta(const size_t index) const{
            return velocityTheta->getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getVelocityTheta(const size_t i, const size_t j) const{
            return velocityTheta->getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getVelocityTheta(const size_t face, const size_t i, const size_t j) const{
            return velocityTheta->getData(face, i, j);
		}

		template<typename T, typename V>
    inline T& airLayer<T, V>::getVelocityPhi(){
			return *velocityPhi;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getVelocityPhi(const size_t index) const{
			return (*velocityPhi).getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getVelocityPhi(const size_t i, const size_t j) const{
			return (*velocityPhi).getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getVelocityPhi(const size_t face, const size_t i, const size_t j) const{
			return velocityPhi->getData(face, i, j);
		}

		template<typename T, typename V>
    inline T& airLayer<T, V>::getTemperature(){
			return *temperature;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getTemperature(const size_t index) const{
			return (*temperature).getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getTemperature(const size_t i, const size_t j) const{
			return (*temperature).getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getTemperature(const size_t face, const size_t i, const size_t j) const{
			return temperature->getData(face, i, j);
		}
		
		template<typename T, typename V>
    inline T& airLayer<T, V>::getMoisture(){
			return *moisture;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getMoisture(const size_t index) const{
			return (*moisture).getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getMoisture(const size_t i, const size_t j) const{
			return (*moisture).getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getMoisture(const size_t face, const size_t i, const size_t j) const{
			return moisture->getData(face, i, j);
		}
		
        template<typename T, typename V>
    inline T& airLayer<T, V>::getCondensedWater(){
            return *condensedWater;
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getCondensedWater(const size_t index) const{
            return (*condensedWater).getData(index);
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getCondensedWater(const size_t i, const size_t j) const{
            return (*condensedWater).getData(i, j);
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getCondensedWater(const size_t face, const size_t i, const size_t j) const{
            return condensedWater->getData(face, i, j);
        }

        template<typename T, typename V>
    inline T& airLayer<T, V>::getClouds(){
            return *cloudTypes;
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getClouds(const size_t index) const{
            return cloudTypes->getData(index);
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getClouds(const size_t i, const size_t j) const{
            return cloudTypes->getData(i, j);
        }

        template<typename T, typename V>
    inline const V airLayer<T, V>::getClouds(const size_t face, const size_t i, const size_t j) const{
            return cloudTypes->getData(face, i, j);
        }
		
		template<typename T, typename V>
    inline T& airLayer<T, V>::getPressure(){
			return *pressure;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getPressure(const size_t index) const{
            return pressure->getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getPressure(const size_t i, const size_t j) const{
            return pressure->getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getPressure(const size_t face, const size_t i, const size_t j) const{
			return pressure->getData(face, i, j);
		}

		template<typename T, typename V>
    inline T& airLayer<T, V>::getParticulates(){
			return *particulates;
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getParticulates(const size_t index) const{
			return (*particulates).getData(index);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getParticulates(const size_t i, const size_t j) const{
			return (*particulates).getData(i, j);
		}

		template<typename T, typename V>
    inline const V airLayer<T, V>::getParticulates(const size_t face, const size_t i, const size_t j) const{
			return particulates->getData(face, i, j);
		}

		/*template<typename T, typename V>
		T& airLayer<T, V>::getOldPressure() const{
			return *oldPressure;
		}

		template<typename T, typename V>
		T& airLayer<T, V>::getOldVelTheta() const{
			return *oldVelTheta;
		}

		template<typename T, typename V>
		T& airLayer<T, V>::getOldVelPhi() const{
			return *oldVelPhi;
		}

		template<typename T, typename V>
		T& airLayer<T, V>::getOldMoisture() const{
			return *oldMoisture;
		}

		template<typename T, typename V>
		T& airLayer<T, V>::getOldTemperature() const{
			return *oldTemp;
		}

		template<typename T, typename V>
		T& airLayer<T, V>::getOldClouds() const{
			return *oldClouds;
		}*/

		template<typename T, typename V>
    inline void airLayer<T, V>::setHeight(V h){
			height = h;
            V tempLess = height / 100.0;
            setTemperature(getPlanet()->getAveTemp() - tempLess);
		}
		
		template<typename T, typename V>
    inline void airLayer<T, V>::setThickness(V t){
			thickness = t;
		}
		
		template<typename T, typename V>
    inline void airLayer<T, V>::setObstacles(T& obs){
			*obstacles = obs;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setObstacles(const size_t index, const V& val){
			(*obstacles).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setObstacles(const size_t i, const size_t j, const V& val){
			(*obstacles).setData(i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setObstacles(const size_t face, const size_t i, const size_t j, const V& val){
			obstacles->setData(face, i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityTheta(T& velTheta){
			*velocityTheta = velTheta;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityTheta(const size_t index, const V& val){
			(*velocityTheta).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityTheta(const size_t i, const size_t j, const V& val){
			(*velocityTheta).setData(i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityTheta(const size_t face, const size_t i, const size_t j, const V& val){
			velocityTheta->setData(face, i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityPhi(T& velPhi){
			*velocityPhi = velPhi;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityPhi(const size_t index, const V& val){
			(*velocityPhi).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityPhi(const size_t i, const size_t j, const V& val){
			(*velocityPhi).setData(i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setVelocityPhi(const size_t face, const size_t i, const size_t j, const V& val){
			velocityPhi->setData(face, i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setTemperature(T& temps){
			*temperature = temps;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setTemperature(const size_t index, const V& val){
			(*temperature).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setTemperature(const size_t i, const size_t j, const V& val){
			(*temperature).setData(i, j, val);
		}
		
		template<typename T, typename V>
    inline void airLayer<T, V>::setTemperature(const size_t face, const size_t i, const size_t j, const V& val){
			temperature->setData(face, i, j, val);
		}

        template<typename T, typename V>
    inline void airLayer<T, V>::setMoisture(T& moist){
            *moisture = moist;
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setMoisture(const size_t index, const V& val){
            (*moisture).setData(index, val);
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setMoisture(const size_t i, const size_t j, const V& val){
            (*moisture).setData(i, j, val);
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setMoisture(const size_t face, const size_t i, const size_t j, const V& val){
            moisture->setData(face, i, j, val);
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setClouds(T& cloudMask){
            *cloudTypes = cloudMask;
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setClouds(const size_t index, const V& val){
            cloudTypes->setData(index, val);
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setClouds(const size_t i, const size_t j, const V& val){
            cloudTypes->setData(i, j, val);
        }

        template<typename T, typename V>
    inline void airLayer<T, V>::setClouds(const size_t face, const size_t i, const size_t j, const V& val){
            cloudTypes->setData(face, i, j, val);
        }

		template<typename T, typename V>
    inline void airLayer<T, V>::setCondensedWater(T& clouds){
			*condensedWater = clouds;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setCondensedWater(const size_t index, const V& val){
			(*condensedWater).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setCondensedWater(const size_t i, const size_t j, const V& val){
			(*condensedWater).setData(i, j, val);
		}
		
		template<typename T, typename V>
    inline void airLayer<T, V>::setCondensedWater(const size_t face, const size_t i, const size_t j, const V& val){
			condensedWater->setData(face, i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setPressure(T& pres){
			*pressure = pres;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setPressure(const size_t index, const V& val){
			(*pressure).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setPressure(const size_t i, const size_t j, const V& val){
			(*pressure).setData(i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setPressure(const size_t face, const size_t i, const size_t j, const V& val){
			pressure->setData(face, i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setParticulates(T& parts){
			*particulates = parts;
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setParticulates(const size_t index, const V& val){
			(*particulates).setData(index, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setParticulates(const size_t i, const size_t j, const V& val){
			(*particulates).setData(i, j, val);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::setParticulates(const size_t face, const size_t i, const size_t j, const V& val){
			particulates->setData(face, i, j, val);
		}

		/*template<typename T, typename V>
		void airLayer<T, V>::setOldPressure(T& oldPres){
			*oldPressure = oldPres;
		}*/

		template<typename T, typename V>
    inline bool airLayer<T, V>::operator==(airLayer<T, V>& other){
			if (getHeight() != other.getHeight())
				return false;
			if (getThickness() != other.getThickness())
				return false;
			if (getObstacles() != other.getObstacles())
				return false;
			if (getVelocityTheta() != other.getVelocityTheta())
				return false;
			if (getVelocityPhi() != other.getVelocityPhi())
				return false;
			if (getTemperature() != other.getTemperature())
				return false;
			if (getMoisture() != other.getMoisture())
				return false;
			if (getCondensedWater() != other.getCondensedWater())
				return false;
			if (getPressure() != other.getPressure())
				return false;
			if (getParticulates() != other.getParticulates())
				return false;
			if (getPlanet() != other.getPlanet()){
				return false;
			}
			return true;
		}

		template<typename T, typename V>
    inline bool airLayer<T, V>::operator!=(airLayer<T, V>& other){
			return !(*this == other);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapScalars(std::vector<std::shared_ptr<T>>& buffer){
			if (buffer.size() != scalarQuantities.size()){
				std::cerr << "Error! Can't swap scalars with vectors of unequal size!" << std::endl;
				return;
			}
			for (int i = 0; i < scalarQuantities.size(); ++i){
				std::swap((*(scalarQuantities[i])), (*(buffer[i])));
			}
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapVels(std::shared_ptr<T>& oldVelTheta, std::shared_ptr<T>& oldVelPhi){
			std::swap(velocityTheta, oldVelTheta);
			std::swap(velocityPhi, oldVelPhi);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapPres(std::shared_ptr<T>& oldPressure){
			std::swap(pressure, oldPressure);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapMoistures(std::shared_ptr<T>& oldMoisture){
			std::swap(moisture, oldMoisture);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapClouds(std::shared_ptr<T>& oldClouds){
			std::swap(condensedWater, oldClouds);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapTemps(std::shared_ptr<T>& oldTemp){
			std::swap(temperature, oldTemp);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::swapParts(std::shared_ptr<T>& oldParticulates){
			std::swap(particulates, oldParticulates);
		}

		template<typename T, typename V>
    inline void airLayer<T, V>::randomInit(){
            velocityTheta->randomInit(-Planet->getMaxWindSpeed() / 2, Planet->getMaxWindSpeed() / 2);
			//oldVelTheta->copy(*velocityTheta);
			
            velocityPhi->randomInit(-Planet->getMaxWindSpeed() / 2, Planet->getMaxWindSpeed() / 2);
			//oldVelPhi->copy(*velocityPhi);
			
			moisture->randomInit(0, Planet->getMaxVapourDensity());
			//oldMoisture->copy(*moisture);
			
			condensedWater->randomInit(0, Planet->getMaxVapourDensity());
			//oldClouds->copy(*condensedWater);

			V dph = PWM::Utils::altitudeAdjustedPressure(height, Planet);
			V aveDiff = Planet->getAverageSealevelPressure() - dph;
			pressure->randomInit(dph - aveDiff, 108000 - aveDiff);//Variations based on known highest and lowest pressures recorded on earth.
			//oldPressure->copy(*pressure);
		}

        template<typename T, typename V>
    inline bool operator<(const airLayer<T, V>& lhs, const airLayer<T, V>& rhs){
            return lhs.getHeight() < rhs.getHeight();
        }

        template<typename T, typename V>
    inline bool operator<(std::shared_ptr<airLayer<T, V>> lhs, std::shared_ptr<airLayer<T, V>> rhs){
            return lhs->getHeight() < rhs->getHeight();
        }

		/* template<typename T, typename V>
		void airLayer<T, V>::solvePolarVelocities(){
			std::cerr << "Error! Not implemented for general template, use a specialized template or define airLayer::solvePolarVelocities() for this template!" << std::endl;
		} */

		/*template<>
        void airLayer<PWM::PWMDataStructure::SSGDataStructure<double>, double>::solvePolarVelocities(){
            std::pair<double, double> NorthPoleBuffer = std::make_pair(0., 0.);
            std::pair<double, double> SouthPoleBuffer = std::make_pair(0., 0.);
            
            size_t minTheta = 0;
            size_t maxTheta = getVelocityTheta().getTheta() - 1;
            #pragma omp parallel for
                for (int i = 0; i < getVelocityTheta().getPhi(); ++i){
                    #if __cplusplus > 201703L
                        double phi = i * (2 * std::numbers::pi / getVelocityTheta().getPhi());
                    #else
                        double phi = i * (2 * 3.1415926535 / getVelocityTheta().getPhi());
                    #endif
                    double ootBeltvPhi = (getVelocityPhi(minTheta, i) + getVelocityPhi(minTheta, i + 1) / 2.);
                    double totBeltvPhi = (getVelocityPhi(minTheta + 1, i) + getVelocityPhi(minTheta + 1, i + 1) / 2.);
                    double vPhiLateralLine = (ootBeltvPhi + totBeltvPhi) / 2.;
                    double vThetaLateralLine = getVelocityTheta(1, i);
                    
                    NorthPoleBuffer.first += vThetaLateralLine * std::cos(phi) - vPhiLateralLine * std::sin(phi);
                    NorthPoleBuffer.second += -vThetaLateralLine * std::sin(phi) + vPhiLateralLine * std::cos(phi);

                    ootBeltvPhi = (getVelocityPhi(maxTheta, i) + getVelocityPhi(maxTheta, i + 1) / 2.);
                    totBeltvPhi = (getVelocityPhi(maxTheta - 1, i) + getVelocityPhi(maxTheta - 1, i + 1) / 2.);
                    vPhiLateralLine = (ootBeltvPhi + totBeltvPhi) / 2.;
                    vThetaLateralLine = getVelocityTheta(maxTheta - 1, i);
                    
                    SouthPoleBuffer.first += -vThetaLateralLine * std::cos(phi) - vPhiLateralLine * std::sin(phi);
                    SouthPoleBuffer.second += -vThetaLateralLine * std::sin(phi) + vPhiLateralLine * std::cos(phi);
                }
            #pragma omp barrier

            //Average out the projected components
            NorthPoleBuffer.first = NorthPoleBuffer.first / getVelocityTheta().getPhi();
            NorthPoleBuffer.second = NorthPoleBuffer.second / getVelocityTheta().getPhi();
            SouthPoleBuffer.first = SouthPoleBuffer.first / getVelocityTheta().getPhi();
            SouthPoleBuffer.second = SouthPoleBuffer.second / getVelocityTheta().getPhi();

            #pragma omp parallel for
                for (int i = 0; i < getVelocityTheta().getPhi(); ++i){
                    #if __cplusplus > 201703L
                        double phi = i * (2 * std::numbers::pi / getVelocityTheta().getPhi());
                    #else
                        double phi = i * (2 * 3.1415926535 / getVelocityTheta().getPhi());
                    #endif

                    double northPole = NorthPoleBuffer.first * std::cos(phi) + NorthPoleBuffer.second * std::sin(phi);
                    double southPole = -SouthPoleBuffer.first * std::cos(phi) - SouthPoleBuffer.second * std::sin(phi);

                    setVelocityTheta(minTheta, i, northPole);
                    setVelocityTheta(maxTheta, i, southPole);
                }
            #pragma omp barrier
        }*/
	}
}
#endif //PWM_MODEL_AIR_LAYER_H
