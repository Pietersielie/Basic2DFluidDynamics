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

#ifndef BFD_ADVECTION_ENGINE_H
#define BFD_ADVECTION_ENGINE_H

//#define DEBUG

#include "abstractEngine.h"
#include "fluidLayer.h"
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include "flatStaggeredGrid.h"
#include "mathUtils.h"
#include <omp.h>
#include <vector>

#define MAX_LAYER 10

namespace BasicFluidDynamics{
    namespace Engine{
        template<typename T, typename V> class advectionEngine : public AbstractEngine{
            private:

                //Advect
                void advect(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);

                //Apply the coriolis force to all cells based on latitude
                void applyCoriolis(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);

                //Data structure to hold reference to each of the layers in the simulation
                std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>> fluidLayers;

                //Buffers for parallel advection:
                std::shared_ptr<T> bufferVelY;
                std::shared_ptr<T> bufferVelX;
                std::shared_ptr<T> bufferCol;

                //A list of the scalar buffers, used for simultaneous advection of all scalar values
                std::vector<std::shared_ptr<T>> scalarBuffers;

                inline void init(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);

                bool coriolis;
                int stepCount = 0;
            protected:
                //The overridden function from AbstractEngine that does the actual advection
                void step_internal() override;
            public:
                //Default constructor
                advectionEngine(float dt = 1.0, bool active = true, bool cor = true);

                advectionEngine(const size_t width, const float dt = 1.0, bool active = true, bool cor = true);

                advectionEngine(const size_t width, const V actualSize, const float dt = 1.0, bool active = true, bool cor = true);

                advectionEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, const float dt, bool active, bool cor = true);

//                //Calculate divergence and pressure (without solving) to simulate uplift/downdrafts to some extent
//                void calcPressures();

//                //Calculate divergence and apply a coarse smoothing to the pressure field
//                void smoothPressure();

                //Function to add an air layer to the engine
                void addLayer(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);
                const std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>& getfluidLayers() const;
                bool getCoriolis() const;
                void setCoriolis(bool cor);
        };

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelY = std::make_shared<T>();
            bufferVelX = std::make_shared<T>();
            bufferCol = std::make_shared<T>();

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelY = std::make_shared<T>(width);
            bufferVelX = std::make_shared<T>(width);
            bufferCol = std::make_shared<T>(width);

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const V actualSize, float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelY = std::make_shared<T>(width, actualSize);
            bufferVelX = std::make_shared<T>(width, actualSize);
            bufferCol = std::make_shared<T>(width, actualSize);

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelY = std::make_shared<T>(width, height, 0.0, 0.5, xWrldSize, yWrldSize);
            bufferVelX = std::make_shared<T>(width, height, 0.5, 0.0, xWrldSize, yWrldSize);
            bufferCol = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);

//            div = Eigen::VectorXd(width * height);
//            pressure = Eigen::VectorXd(width * height);
//            PrevPressure = Eigen::VectorXd(width * height);

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::init(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::init() for this template!" << std::endl;
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::addLayer(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            fluidLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>& advectionEngine<T, V>::getfluidLayers() const{
            return fluidLayers;
        }

        template<typename T, typename V>
        inline bool advectionEngine<T, V>::getCoriolis() const{
            return coriolis;
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::setCoriolis(bool cor){
            coriolis = cor;
        }

        /**
         * The body forces that act upon a rotating frame of reference are the Euler force, coriolis force, and centrifugal force.
         * Equation is derived from F = ma. Full equation is explained at https://en.wikipedia.org/wiki/Coriolis_force.
         * Since the planet rotates uniformly based on distance from rotational axis with no acceleration, we can ignore Euler force.
         * Similarly, since we're only dealing with the lower atmosphere at relatively low speeds, centrifugal forces can be ignored.
         * Thus, final equation is F = ma + 2m(ω * v). We ignore all but horizontal terms for our system (given the neglible
         * vertical velocities in earth's fluid systems).
         *
         * Therefore, our final equation is deltaV = (-vPhi, vTheta) * 2ωsinφ, where ω is the spin rate, and φ is the latitude.
         * deltaV is the change in velocity in vector form, with the first element being the north-south direction, and second the east west direction.
         * This force is only applied away from the poles, since there is no rotation at the poles.
         */
        template<typename T, typename V>
        inline void advectionEngine<T, V>::applyCoriolis(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::applyCoriolis(BasicFluidDynamics::Model::fluidLayer<T, V>& l) for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::applyCoriolis(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double omega = 0;
//            double omega = l->getPlanet()->getAngularVelocity();//Omega is the angular velocity of our rotating frame of reference.
            auto ObsRef = l->getObsPtr();
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        //Latitude is converted into radians for the sin() function.
                        double lat = BasicFluidDynamics::Utils::degToRad(l->getObstacles().getCoordinates(i, j).getLatitude());
                        double F = 2 * omega * std::sin(lat);//F is the Coriolis parameter

                        auto VelYLoc = l->getVelocityY().getWorldLoc(i, j);
                        auto VelXLoc = l->getVelocityX().getWorldLoc(i, j);

                        //Calculate change in velocity due to acceleration from coriolis force
                        double changeVelY = l->getVelocityY(i, j) + (F * (-l->getVelocityX().sampleAt(VelYLoc, ObsRef)) * getDt());
                        double changeVelX = l->getVelocityX(i, j) + (F * l->getVelocityY().sampleAt(VelXLoc, ObsRef) * getDt());

                        //Write data to buffers
                        bufferVelY->setData(i, j, changeVelY);
                        bufferVelX->setData(i, j, changeVelX);
                    }
                }
            #pragma omp barrier
            l->swapVels(bufferVelY, bufferVelX);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::advect(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::advect() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::advect(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            auto VelYOffset = l->getVelocityY().getOffset();
            auto VelXOffset = l->getVelocityX().getOffset();
            auto centreOffset = l->getColour().getOffset();
            double zero = 0;
            auto dt = getDt();
            auto ObsRef = l->getObsPtr();
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){

                        //Advect the velocity field in the x direction (left right)
                        auto VelXLoc = l->getVelocityX().getWorldLoc(i, j);
                        double deltaVelXI = l->getVelocityY().sampleAt(VelXLoc, ObsRef) * -dt / cellSizeX;
                        double deltaVelXJ = l->getVelocityX().sampleAt(VelXLoc, ObsRef) * -dt / cellSizeY;
                        double newVelXI = (i + deltaVelXI) + VelXOffset.first;
                        double newVelXJ = (j + deltaVelXJ) + VelXOffset.second;
                        auto newVelXLoc = std::make_pair(newVelXI, newVelXJ);
                        bufferVelX->setData(i, j, l->getVelocityX().sampleAt(newVelXLoc, ObsRef));

                        //Advect the velocity field in the y direction (top bottom)
                        auto VelYLoc = l->getVelocityY().getWorldLoc(i, j);
                        double deltaVelYI = l->getVelocityY().sampleAt(VelYLoc, ObsRef) * -dt / cellSizeX;
                        double deltaVelYJ = l->getVelocityX().sampleAt(VelYLoc, ObsRef) * -dt / cellSizeY;
                        double newVelYI = (i + deltaVelYI) + VelYOffset.first;
                        double newVelYJ = (j + deltaVelYJ) + VelYOffset.second;
                        auto newVelYLoc = std::make_pair(newVelYI, newVelYJ);
                        bufferVelY->setData(i, j, l->getVelocityY().sampleAt(newVelYLoc, ObsRef));

                        //Check if the cell contents can be advected
                        if (l->getObstacles(i, j)){
                            bufferCol->setData(i, j, l->getColour(i, j));
                            continue;
                        }
                        //Advect the cell contents
                        auto centreLoc = l->getObstacles().getWorldLoc(i, j);
                        double deltaCentreI = l->getVelocityY().sampleAt(centreLoc, ObsRef) * -dt / cellSizeX;
                        double deltaCentreJ = l->getVelocityX().sampleAt(centreLoc, ObsRef) * -dt / cellSizeY;
                        double newCentreI = (i + deltaCentreI) + centreOffset.first;
                        double newCentreJ = (j + deltaCentreJ) + centreOffset.second;
                        auto newCentreLoc = std::make_pair(newCentreI, newCentreJ);

                        double newCol = l->getColour().sampleAt(newCentreLoc, ObsRef);
                        bufferCol->setData(i, j, newCol);
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelY, bufferVelX);
            l->swapCols(bufferCol);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::step_internal(){
            for(auto l : fluidLayers){
                // std::cout << "X velocities prior to advection:\n" << l->getVelocityY().print() << std::endl;
                advect(l);

                // std::cout << "X velocities after advection:\n" << l->getVelocityY().print() << std::endl;
                applyCoriolis(l);
            }
            ++stepCount;
        }
    }
}
#endif //BFD_ADVECTION_ENGINE_H
