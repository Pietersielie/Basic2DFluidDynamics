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
                std::shared_ptr<T> bufferVelTheta;
                std::shared_ptr<T> bufferVelPhi;
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

            bufferVelTheta = std::make_shared<T>();
            bufferVelPhi = std::make_shared<T>();
            bufferCol = std::make_shared<T>();

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>(width);
            bufferVelPhi = std::make_shared<T>(width);
            bufferCol = std::make_shared<T>(width);

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const V actualSize, float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>(width, actualSize);
            bufferVelPhi = std::make_shared<T>(width, actualSize);
            bufferCol = std::make_shared<T>(width, actualSize);

            scalarBuffers.push_back(bufferCol);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>(width, height, 0.0, 0.5, xWrldSize, yWrldSize);
            bufferVelPhi = std::make_shared<T>(width, height, 0.5, 0.0, xWrldSize, yWrldSize);
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

                        auto velThetaLoc = l->getVelocityTheta().getWorldLoc(i, j);
                        auto velPhiLoc = l->getVelocityPhi().getWorldLoc(i, j);

                        //Calculate change in velocity due to acceleration from coriolis force
                        double changeVelTheta = l->getVelocityTheta(i, j) + (F * (-l->getVelocityPhi().sampleAt(velThetaLoc, ObsRef)) * getDt());
                        double changeVelPhi = l->getVelocityPhi(i, j) + (F * l->getVelocityTheta().sampleAt(velPhiLoc, ObsRef) * getDt());

                        //Write data to buffers
                        bufferVelTheta->setData(i, j, changeVelTheta);
                        bufferVelPhi->setData(i, j, changeVelPhi);
                    }
                }
            #pragma omp barrier
            l->swapVels(bufferVelTheta, bufferVelPhi);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::advect(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::advect() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::advect(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            auto velThetaOffset = l->getVelocityTheta().getOffset();
            auto velPhiOffset = l->getVelocityPhi().getOffset();
            auto centreOffset = l->getColour().getOffset();
            double zero = 0;
            auto dt = getDt();
            auto ObsRef = l->getObsPtr();
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){

                        //Advect the velocity field in the x direction (left right)
                        auto velPhiLoc = l->getVelocityPhi().getWorldLoc(i, j);
                        double deltaVelPhiI = l->getVelocityTheta().sampleAt(velPhiLoc, ObsRef) * -dt / cellSizeX;
                        double deltaVelPhiJ = l->getVelocityPhi().sampleAt(velPhiLoc, ObsRef) * -dt / cellSizeY;
                        double newVelPhiI = (i + deltaVelPhiI) + velPhiOffset.first;
                        double newVelPhiJ = (j + deltaVelPhiJ) + velPhiOffset.second;
                        auto newVelPhiLoc = std::make_pair(newVelPhiI, newVelPhiJ);
                        bufferVelPhi->setData(i, j, l->getVelocityPhi().sampleAt(newVelPhiLoc, ObsRef));

                        //Advect the velocity field in the y direction (top bottom)
                        auto velThetaLoc = l->getVelocityTheta().getWorldLoc(i, j);
                        double deltaVelThetaI = l->getVelocityTheta().sampleAt(velThetaLoc, ObsRef) * -dt / cellSizeX;
                        double deltaVelThetaJ = l->getVelocityPhi().sampleAt(velThetaLoc, ObsRef) * -dt / cellSizeY;
                        double newVelThetaI = (i + deltaVelThetaI) + velThetaOffset.first;
                        double newVelThetaJ = (j + deltaVelThetaJ) + velThetaOffset.second;
                        auto newVelThetaLoc = std::make_pair(newVelThetaI, newVelThetaJ);
                        bufferVelTheta->setData(i, j, l->getVelocityTheta().sampleAt(newVelThetaLoc, ObsRef));

                        //Check if the cell contents can be advected
                        if (l->getObstacles(i, j)){
                            bufferCol->setData(i, j, l->getColour(i, j));
                            continue;
                        }
                        //Advect the cell contents
                        auto centreLoc = l->getObstacles().getWorldLoc(i, j);
                        double deltaCentreI = l->getVelocityTheta().sampleAt(centreLoc, ObsRef) * -dt / cellSizeX;
                        double deltaCentreJ = l->getVelocityPhi().sampleAt(centreLoc, ObsRef) * -dt / cellSizeY;
                        double newCentreI = (i + deltaCentreI) + centreOffset.first;
                        double newCentreJ = (j + deltaCentreJ) + centreOffset.second;
                        auto newCentreLoc = std::make_pair(newCentreI, newCentreJ);

                        double newCol = l->getColour().sampleAt(newCentreLoc, ObsRef);
                        bufferCol->setData(i, j, newCol);
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
            l->swapCols(bufferCol);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::step_internal(){
            for(auto l : fluidLayers){
                // std::cout << "X velocities prior to advection:\n" << l->getVelocityTheta().print() << std::endl;
                advect(l);

                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
                applyCoriolis(l);
            }
            ++stepCount;
        }
    }
}
#endif //PWM_ADVECTION_ENGINE_H
