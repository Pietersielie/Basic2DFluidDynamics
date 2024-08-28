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

#ifndef BFD_MODEL_FLUID_LAYER_H
#define BFD_MODEL_FLUID_LAYER_H

#include <algorithm>
#include <cstddef>
#include "flatStaggeredGrid.h"
#include <iostream>
#include <memory>
#include <vector>

namespace BasicFluidDynamics{
	namespace Model{
        template<typename T, typename V> class fluidLayer{
            private:
			
				// Storing all the various bits and bobs that make up an air layer
				std::shared_ptr<T> obstacles;//obstacles, i.e. land
                std::shared_ptr<T> velocityY;//velocity in the north-south axis (meters per second)
                std::shared_ptr<T> velocityX;//velocity in the east-west axis (meters per second)
                std::shared_ptr<T> colour;//temperature in Kelvin
                std::shared_ptr<T> pressure;//modified pressure, not true pascal values

                V density = 1.204; //density in kg per cubic metre -- 1.204 for air, 1000 for water
			public:
				std::vector<std::shared_ptr<T>> scalarQuantities;
				std::vector<std::shared_ptr<T>> oldScalarQuantities;
				
				//Various constructors for different parameter inputs.
                fluidLayer();

                //Constructor meant to be used with flatStaggeredGrid
                fluidLayer(const size_t width, const size_t height, const V xSize, const V ySize);

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

                void swapCols(std::shared_ptr<T>& oldCols);

                // Uniformly set the colour of a layer
                void setColour(V heat);
                // Compute the average colour in the layer
                V getMeanColour() const;

                void setDensity(const V rho);
                V getDensity() const;

				T& getObstacles();
                std::shared_ptr<T>& getObsPtr();
                T& getVelocityY();
                T& getVelocityX();
                T& getColour();
                T& getPressure();

				const V getObstacles(const size_t index) const;
				const V getObstacles(const size_t i, const size_t j) const;
                const V getVelocityY(const size_t index) const;
                const V getVelocityY(const size_t i, const size_t j) const;
                const V getVelocityX(const size_t index) const;
                const V getVelocityX(const size_t i, const size_t j) const;
                const V getColour(const size_t index) const;
                const V getColour(const size_t i, const size_t j) const;
                const V getPressure(const size_t index) const;
                const V getPressure(const size_t i, const size_t j) const;

				// Setters for each of the data fields
                void setObstacles(T& obs);
                void setVelocityY(T& vTheta);
                void setVelocityX(T& vPhi);
                void setColour(T& temps);
                void setPressure(T& pres);

				void setObstacles(const size_t index, const V& val);
				void setObstacles(const size_t i, const size_t j, const V& val);
                void setVelocityY(const size_t index, const V& val);
                void setVelocityY(const size_t i, const size_t j, const V& val);
                void setVelocityX(const size_t index, const V& val);
                void setVelocityX(const size_t i, const size_t j, const V& val);
                void setColour(const size_t index, const V& val);
                void setColour(const size_t i, const size_t j, const V& val);
                void setPressure(const size_t index, const V& val);
				void setPressure(const size_t i, const size_t j, const V& val);

				// Various comparators for an air layer
                bool operator==(fluidLayer& other);
                bool operator!=(fluidLayer& other);

				template<typename X, typename Y>
                friend bool operator<(const fluidLayer<X, Y>& lhs, const fluidLayer<X, Y>& rhs);
                template<typename X, typename Y>
                friend bool operator<(std::shared_ptr<fluidLayer<X, Y>> lhs, std::shared_ptr<fluidLayer<X, Y>> rhs);
        };
		
		/**int i = 0; i < quantity.size(); ++i
		 * Implementing the templated functions above on general rules... fingers crossed.
		 */
		template<typename T, typename V>
        inline fluidLayer<T, V>::fluidLayer(){
			obstacles = std::make_shared<T>();
            velocityY = std::make_shared<T>();
            velocityX = std::make_shared<T>();
            colour = std::make_shared<T>();
            pressure = std::make_shared<T>();

            scalarQuantities.push_back(colour);
        }

        /**
         * This constructor is specifically meant for using flatStaggeredGrid<V> for the various fields.
         */
        template<typename T, typename V>
        inline fluidLayer<T, V>::fluidLayer(const size_t width, const size_t height, const V xSize, const V ySize){
            obstacles = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            velocityY = std::make_shared<T>(width, height, 0.0, 0.5, xSize, ySize);
            velocityX = std::make_shared<T>(width, height, 0.5, 0.0, xSize, ySize);
            colour = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);
            pressure = std::make_shared<T>(width, height, 0.5, 0.5, xSize, ySize);

            scalarQuantities.push_back(colour);
            init();
        }

		template<typename T, typename V>
        inline void fluidLayer<T, V>::init(){
			*obstacles = 0;
            *velocityY = 0;
            *velocityX = 0;
            *colour = 0;
		}
		
		template<typename T, typename V>
        inline void fluidLayer<T, V>::setColour(V heat){
            *colour = heat;
		}
		
		template<typename T, typename V>
        inline V fluidLayer<T, V>::getMeanColour() const{
            return colour->mean();
        }

        template<typename T, typename V>
        inline void fluidLayer<T, V>::setDensity(const V rho){
            density = rho;
        }

        template<typename T, typename V>
        inline V fluidLayer<T, V>::getDensity() const{
            return density;
        }

        template<typename T, typename V>
        inline T& fluidLayer<T, V>::getObstacles(){
            return *obstacles;
        }

        template<typename T, typename V>
        inline std::shared_ptr<T>& fluidLayer<T, V>::getObsPtr(){
            return obstacles;
        }

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getObstacles(const size_t index) const{
			return (*obstacles).getData(index);
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getObstacles(const size_t i, const size_t j) const{
			return (*obstacles).getData(i, j);
		}

        template<typename T, typename V>
        inline T& fluidLayer<T, V>::getVelocityY(){
            return *velocityY;
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getVelocityY(const size_t index) const{
            return velocityY->getData(index);
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getVelocityY(const size_t i, const size_t j) const{
            return velocityY->getData(i, j);
		}

        template<typename T, typename V>
        inline T& fluidLayer<T, V>::getVelocityX(){
            return *velocityX;
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getVelocityX(const size_t index) const{
            return (*velocityX).getData(index);
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getVelocityX(const size_t i, const size_t j) const{
            return (*velocityX).getData(i, j);
		}

        template<typename T, typename V>
        inline T& fluidLayer<T, V>::getColour(){
            return *colour;
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getColour(const size_t index) const{
            return (*colour).getData(index);
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getColour(const size_t i, const size_t j) const{
            return (*colour).getData(i, j);
		}

		template<typename T, typename V>
        inline T& fluidLayer<T, V>::getPressure(){
			return *pressure;
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getPressure(const size_t index) const{
            return pressure->getData(index);
		}

		template<typename T, typename V>
        inline const V fluidLayer<T, V>::getPressure(const size_t i, const size_t j) const{
            return pressure->getData(i, j);
		}
		
		template<typename T, typename V>
        inline void fluidLayer<T, V>::setObstacles(T& obs){
			*obstacles = obs;
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setObstacles(const size_t index, const V& val){
			(*obstacles).setData(index, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setObstacles(const size_t i, const size_t j, const V& val){
			(*obstacles).setData(i, j, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setVelocityY(T& velTheta){
            *velocityY = velTheta;
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setVelocityY(const size_t index, const V& val){
            (*velocityY).setData(index, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setVelocityY(const size_t i, const size_t j, const V& val){
            (*velocityY).setData(i, j, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setVelocityX(T& velPhi){
            *velocityX = velPhi;
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setVelocityX(const size_t index, const V& val){
            (*velocityX).setData(index, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setVelocityX(const size_t i, const size_t j, const V& val){
            (*velocityX).setData(i, j, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setColour(T& temps){
            *colour = temps;
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setColour(const size_t index, const V& val){
            (*colour).setData(index, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setColour(const size_t i, const size_t j, const V& val){
            (*colour).setData(i, j, val);
		}
		
		template<typename T, typename V>
        inline void fluidLayer<T, V>::setPressure(T& pres){
            *pressure = pres;
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setPressure(const size_t index, const V& val){
			(*pressure).setData(index, val);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::setPressure(const size_t i, const size_t j, const V& val){
			(*pressure).setData(i, j, val);
		}

		template<typename T, typename V>
        inline bool fluidLayer<T, V>::operator==(fluidLayer<T, V>& other){
			if (getObstacles() != other.getObstacles())
				return false;
            if (getVelocityY() != other.getVelocityY())
				return false;
            if (getVelocityX() != other.getVelocityX())
				return false;
            if (getColour() != other.getColour())
                return false;
			return true;
		}

		template<typename T, typename V>
        inline bool fluidLayer<T, V>::operator!=(fluidLayer<T, V>& other){
			return !(*this == other);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::swapScalars(std::vector<std::shared_ptr<T>>& buffer){
			if (buffer.size() != scalarQuantities.size()){
				std::cerr << "Error! Can't swap scalars with vectors of unequal size!" << std::endl;
				return;
			}
			for (int i = 0; i < scalarQuantities.size(); ++i){
				std::swap((*(scalarQuantities[i])), (*(buffer[i])));
			}
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::swapVels(std::shared_ptr<T>& oldVelTheta, std::shared_ptr<T>& oldVelPhi){
            std::swap(velocityY, oldVelTheta);
            std::swap(velocityX, oldVelPhi);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::swapPres(std::shared_ptr<T>& oldPressure){
			std::swap(pressure, oldPressure);
		}

		template<typename T, typename V>
        inline void fluidLayer<T, V>::swapCols(std::shared_ptr<T>& oldTemp){
            std::swap(colour, oldTemp);
        }

        template<typename T, typename V>
        inline bool operator<(const fluidLayer<T, V>& lhs, const fluidLayer<T, V>& rhs){
            return lhs.getHeight() < rhs.getHeight();
        }

        template<typename T, typename V>
        inline bool operator<(std::shared_ptr<fluidLayer<T, V>> lhs, std::shared_ptr<fluidLayer<T, V>> rhs){
            return lhs->getHeight() < rhs->getHeight();
        }
	}
}
#endif //PWM_MODEL_AIR_LAYER_H
