#ifndef PWM_ADVECTION_ENGINE_H
#define PWM_ADVECTION_ENGINE_H

//#define DEBUG

#include "abstractEngine.h"
#include "airLayer.h"
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include "flatStaggeredGrid.h"
#include "mathUtils.h"
#include <omp.h>
#include "square2DArray.h"
#include <vector>
#include "vizUtils.h"

#define MAX_LAYER 10

namespace PWM{
    namespace Engine{
        template<typename T, typename V> class advectionEngine : public AbstractEngine{
            private:

                //Advect
                void advect(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

//                //Calculate and apply advection for all velocities (offset from cell centers)
//                void advectVelocities(PWM::Model::airLayer<T, V>& l);
//
//                //Calculate and apply advection for all scalar quantities (stored in cell)
//                void advectScalarQuantities(PWM::Model::airLayer<T, V>& l);

                //Apply the coriolis force to all cells based on latitude
                void applyCoriolis(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

//                //Apply the geometric terms resulting from the spherical application of advection
//                void applyGeometricTerms(PWM::Model::airLayer<T, V>& l);

//                //Function that calculates the divergence of the layer
//                void calcDiv(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

//                void multiGridSolve(T& pressure, const T& divergence, int iteration = 0);

//                //Function that calculates the pressures needed to correct for nabla u = 0.
//                void solvePressureMatrix(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

//                //Solve the spherical Poisson equation using the Intel Math Kernel Library
//                void sphericalPoissonSolver(PWM::Model::airLayer<T, V>& l);

//                //Calculate and apply the divergence correction of velocities
//                void solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<T, V>>& l, V scaleFactor = 1);

                //Data structure to hold reference to each of the layers in the simulation
                std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>> airLayers;

                //Buffers for parallel advection:
                std::shared_ptr<T> bufferVelTheta;
                std::shared_ptr<T> bufferVelPhi;
                std::shared_ptr<T> bufferTemp;
                std::shared_ptr<T> bufferMoisture;
                std::shared_ptr<T> bufferClouds;
                std::shared_ptr<T> bufferPressure;
                std::shared_ptr<T> bufferParticulates;

                std::shared_ptr<T> divergence;

//                Eigen::VectorXd div;
//                Eigen::VectorXd pressure;
//                Eigen::VectorXd PrevPressure;
//                Eigen::MatrixXd ATest;

                //A list of the scalar buffers, used for simultaneous advection of all scalar values
                std::vector<std::shared_ptr<T>> scalarBuffers;

                inline void init(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);

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
                void addLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l);
                const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& getAirLayers() const;
                bool getCoriolis() const;
                void setCoriolis(bool cor);
        };

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>();
            bufferVelPhi = std::make_shared<T>();
            bufferTemp = std::make_shared<T>();
            bufferMoisture = std::make_shared<T>();
            bufferClouds = std::make_shared<T>();
            bufferPressure = std::make_shared<T>();
            bufferParticulates = std::make_shared<T>();
            divergence = std::make_shared<T>();

            scalarBuffers.push_back(bufferTemp);
            scalarBuffers.push_back(bufferMoisture);
            scalarBuffers.push_back(bufferClouds);
            scalarBuffers.push_back(bufferPressure);
            scalarBuffers.push_back(bufferParticulates);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>(width);
            bufferVelPhi = std::make_shared<T>(width);
            bufferTemp = std::make_shared<T>(width);
            bufferMoisture = std::make_shared<T>(width);
            bufferClouds = std::make_shared<T>(width);
            bufferPressure = std::make_shared<T>(width);
            bufferParticulates = std::make_shared<T>(width);
            divergence = std::make_shared<T>(width);

            scalarBuffers.push_back(bufferTemp);
            scalarBuffers.push_back(bufferMoisture);
            scalarBuffers.push_back(bufferClouds);
            scalarBuffers.push_back(bufferPressure);
            scalarBuffers.push_back(bufferParticulates);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const V actualSize, float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>(width, actualSize);
            bufferVelPhi = std::make_shared<T>(width, actualSize);
            bufferTemp = std::make_shared<T>(width, actualSize);
            bufferMoisture = std::make_shared<T>(width, actualSize);
            bufferClouds = std::make_shared<T>(width, actualSize);
            bufferPressure = std::make_shared<T>(width, actualSize);
            bufferParticulates = std::make_shared<T>(width, actualSize);
            divergence = std::make_shared<T>(width, actualSize);

            scalarBuffers.push_back(bufferTemp);
            scalarBuffers.push_back(bufferMoisture);
            scalarBuffers.push_back(bufferClouds);
            scalarBuffers.push_back(bufferPressure);
            scalarBuffers.push_back(bufferParticulates);
        }

        template<typename T, typename V>
        inline advectionEngine<T, V>::advectionEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, float dt, bool active, bool cor): AbstractEngine(dt, active), coriolis(cor) {
            airLayers = std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>(width, height, 0.0, 0.5, xWrldSize, yWrldSize);
            bufferVelPhi = std::make_shared<T>(width, height, 0.5, 0.0, xWrldSize, yWrldSize);
            bufferTemp = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);
            bufferMoisture = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);
            bufferClouds = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);
            bufferPressure = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);
            bufferParticulates = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);
            divergence = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);

//            div = Eigen::VectorXd(width * height);
//            pressure = Eigen::VectorXd(width * height);
//            PrevPressure = Eigen::VectorXd(width * height);

            scalarBuffers.push_back(bufferTemp);
            scalarBuffers.push_back(bufferMoisture);
            scalarBuffers.push_back(bufferClouds);
            scalarBuffers.push_back(bufferPressure);
            scalarBuffers.push_back(bufferParticulates);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::init(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::init() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::init(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l){
            size_t faceWidth = l->getObstacles().getWidth();

            bufferVelTheta = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            bufferVelPhi = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            bufferTemp = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            bufferMoisture = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            bufferClouds = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            bufferPressure = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            bufferParticulates = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);
            divergence = std::make_shared<PWM::PWMDataStructure::square2DArray<double>>(faceWidth);

            scalarBuffers.clear();
            scalarBuffers.push_back(bufferTemp);
            scalarBuffers.push_back(bufferMoisture);
            scalarBuffers.push_back(bufferClouds);
            scalarBuffers.push_back(bufferPressure);
            scalarBuffers.push_back(bufferParticulates);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::addLayer(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            airLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<PWM::Model::airLayer<T, V>>>& advectionEngine<T, V>::getAirLayers() const{
            return airLayers;
        }

        template<typename T, typename V>
        inline bool advectionEngine<T, V>::getCoriolis() const{
            return coriolis;
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::setCoriolis(bool cor){
            coriolis = cor;
        }

        /*template<typename T, typename V>
        void advectionEngine<T, V>::calcDiv(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::calcDiv(PWM::Model::airLayer<T, V>& l) for this template!" << std::endl;
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::calcDiv(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l){
            int cellNum = l->getObstacles().getWidth();
            double cellSize = 1.0f / cellNum;
            double velSolid = 0;
            #pragma omp parallel for
                for (int i = 0; i < cellNum; ++i){
                    for (int j = 0; j < cellNum; ++j){
                        if (l->getObstacles(i, j)){
                            divergence->setData(i, j, 0);
                            continue;
                        }
                        double velBot, velTop, velRight, velLeft;
                        if (l->getObstacles(i + 1, j))
                            velBot = velSolid;
                        else
                            velBot = l->getVelocityTheta(i + 1, j);
                        if (l->getObstacles(i - 1, j))
                            velTop = velSolid;
                        else
                            velTop = l->getVelocityTheta(i - 1, j);
                        if (l->getObstacles(i, j - 1))
                            velLeft = velSolid;
                        else
                            velLeft = l->getVelocityPhi(i, j - 1);
                        if (l->getObstacles(i, j + 1))
                            velRight = velSolid;
                        else
                            velRight = l->getVelocityPhi(i, j + 1);
                        double div = ((velBot - velTop) + (velRight - velLeft)) / 2.f / cellSize;
                        divergence->setData(i, j, div);
                    }
                }
            #pragma omp barrier
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::calcDiv(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double velSolid = 0;
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        if (l->getObstacles(i, j)){
                            divergence->setData(i, j, 0);
                            continue;
                        }
                        double velBot, velTop, velRight, velLeft;
                        velTop = (l->getObstacles(i - 1, j)) ? velSolid : l->getVelocityTheta(i, j);
                        velBot = (l->getObstacles(i + 1, j)) ? velSolid : l->getVelocityTheta(i + 1, j);
                        velLeft = (l->getObstacles(i, j - 1)) ? velSolid : l->getVelocityPhi(i, j);
                        velRight = (l->getObstacles(i, j + 1)) ? velSolid : l->getVelocityPhi(i, j + 1);
                        double div = (velBot - velTop) / cellSizeY + (velRight - velLeft) / cellSizeX;
                        divergence->setData(i, j, div);
                    }
                }
        }*/

/*        template<typename T, typename V>
        void advectionEngine<T, V>::advectScalarQuantities(PWM::Model::airLayer<T, V>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::advectScalarQuantities(PWM::Model::airLayer<T, V>& l) for this template!" << std::endl;
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::advectVelocities(PWM::Model::airLayer<T, V>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::advectVelocities(PWM::Model::airLayer<T, V>& l) for this template!" << std::endl;
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::applyGeometricTerms(PWM::Model::airLayer<T, V>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::applyGeometricTerms(PWM::Model::airLayer<T, V>& l) for this template!" << std::endl;
        }*/

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
        inline void advectionEngine<T, V>::applyCoriolis(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::applyCoriolis(PWM::Model::airLayer<T, V>& l) for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::applyCoriolis(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l){
            double omega = l->getPlanet()->getAngularVelocity();//Omega is the angular velocity of our rotating frame of reference.
            #pragma omp parallel for
                for (int i = 0; i < l->getObstacles().getWidth(); ++i){
                    for (int j = 0; j < l->getObstacles().getWidth(); ++j){
                        //Latitude is converted into radians for the sin() function.
                        #if __cplusplus > 201703L
                            double lat = (l->getObstacles().getCoordinates(i, j).getLatitude() * std::numbers::pi) / 180;
                        #else
                            double lat = (l->getObstacles().getCoordinates(i, j).getLatitude() * 3.1415926535) / 180;
                        #endif

                        double f = 2 * omega * std::sin(lat);//F is the Coriolis parameter

                        //Calculate change in velocity due to acceleration from coriolis force
                        double changeVelTheta = l->getVelocityTheta(i, j) + (f * (-l->getVelocityPhi(i, j)) * getDt());
                        double changeVelPhi = l->getVelocityPhi(i, j) + (f * l->getVelocityTheta(i, j) * getDt());

                        //Write data to buffers
                        l->getVelocityTheta().setData(i, j, changeVelTheta);
                        l->getVelocityPhi().setData(i, j, changeVelPhi);
                    }
                }
            #pragma omp barrier
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::applyCoriolis(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double omega = l->getPlanet()->getAngularVelocity();//Omega is the angular velocity of our rotating frame of reference.
            auto ObsRef = l->getObsPtr();
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        //Latitude is converted into radians for the sin() function.
                        #if __cplusplus > 201703L
                            double lat = (l->getObstacles().getCoordinates(i, j).getLatitude() * std::numbers::pi) / 180;
                        #else
                            double lat = (l->getObstacles().getCoordinates(i, j).getLatitude() * 3.1415926535) / 180;
                        #endif

                        double f = 2 * omega * std::sin(lat);//F is the Coriolis parameter

                        auto velThetaLoc = l->getVelocityTheta().getWorldLoc(i, j);
                        auto velPhiLoc = l->getVelocityPhi().getWorldLoc(i, j);

                        //Calculate change in velocity due to acceleration from coriolis force
                        double changeVelTheta = l->getVelocityTheta(i, j) + (f * (-l->getVelocityPhi().sampleAt(velThetaLoc, ObsRef)) * getDt());
                        double changeVelPhi = l->getVelocityPhi(i, j) + (f * l->getVelocityTheta().sampleAt(velPhiLoc, ObsRef) * getDt());

                        //Write data to buffers
                        bufferVelTheta->setData(i, j, changeVelTheta);
                        bufferVelPhi->setData(i, j, changeVelPhi);
                    }
                }
            #pragma omp barrier
            l->swapVels(bufferVelTheta, bufferVelPhi);
        }

        template<typename T, typename V>
        inline void advectionEngine<T, V>::advect(std::shared_ptr<PWM::Model::airLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::advect() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::advect(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l){
            int cellNum = l->getObstacles().getWidth();
            double cellSize = l->getObstacles().getSize() / cellNum;
            double zero = 0;
            #pragma omp parallel for
                for (int i = 0; i < cellNum; ++i){
                    for (int j = 0; j < cellNum; ++j){
                        if (l->getObstacles(i, j)){
                            bufferVelTheta->setData(i, j, l->getVelocityTheta(i, j));
                            bufferVelPhi->setData(i, j, l->getVelocityPhi(i, j));
                            bufferTemp->setData(i, j, l->getTemperature(i, j));
                            bufferMoisture->setData(i, j, l->getMoisture(i, j));
                            bufferClouds->setData(i, j, l->getCondensedWater(i, j));
                            bufferParticulates->setData(i, j, l->getParticulates(i, j));
                            continue;
                        }
                        double deltaX = l->getVelocityTheta(i, j) * getDt() / cellSize;
                        double deltaY = l->getVelocityPhi(i, j) * getDt() / cellSize;
                        double newX = (i - deltaX);// / cellNum;
                        double newY = (j - deltaY);// / cellNum;

                        bufferVelTheta->setData(i, j, l->getVelocityTheta().getInterpolated(newX, newY));
                        bufferVelPhi->setData(i, j, l->getVelocityPhi().getInterpolated(newX, newY));
                        bufferTemp->setData(i, j, l->getTemperature().getInterpolated(newX, newY));
                        bufferMoisture->setData(i, j, l->getMoisture().getInterpolated(newX, newY));
                        bufferClouds->setData(i, j, l->getCondensedWater().getInterpolated(newX, newY));
                        bufferParticulates->setData(i, j, l->getParticulates().getInterpolated(newX, newY));
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
            l->swapTemps(bufferTemp);
            l->swapMoistures(bufferMoisture);
            l->swapClouds(bufferClouds);
            l->swapParts(bufferParticulates);
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::advect(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>& l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            auto velThetaOffset = l->getVelocityTheta().getOffset();
            auto velPhiOffset = l->getVelocityPhi().getOffset();
            auto centreOffset = l->getTemperature().getOffset();
            double zero = 0;
            auto dt = getDt();
            auto ObsRef = l->getObsPtr();
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
//                        if (i == 167 && j == 467)
//                            std::cerr << "Debug point reached" << std::endl;

                        //Advect the velocity field in the y direction (top bottom)
                        auto velThetaLoc = l->getVelocityTheta().getWorldLoc(i, j);
                        double deltaVelThetaI = l->getVelocityTheta().sampleAt(velThetaLoc, ObsRef) * -dt / cellSizeX;
                        double deltaVelThetaJ = l->getVelocityPhi().sampleAt(velThetaLoc, ObsRef) * -dt / cellSizeY;
                        double newVelThetaI = (i + deltaVelThetaI) + velThetaOffset.first;
                        double newVelThetaJ = (j + deltaVelThetaJ) + velThetaOffset.second;
                        auto newVelThetaLoc = std::make_pair(newVelThetaI, newVelThetaJ);
                        bufferVelTheta->setData(i, j, l->getVelocityTheta().sampleAt(newVelThetaLoc, ObsRef));

                        //Advect the velocity field in the x direction (left right)
                        auto velPhiLoc = l->getVelocityPhi().getWorldLoc(i, j);
                        double deltaVelPhiI = l->getVelocityTheta().sampleAt(velPhiLoc, ObsRef) * -dt / cellSizeX;
                        double deltaVelPhiJ = l->getVelocityPhi().sampleAt(velPhiLoc, ObsRef) * -dt / cellSizeY;
                        double newVelPhiI = (i + deltaVelPhiI) + velPhiOffset.first;
                        double newVelPhiJ = (j + deltaVelPhiJ) + velPhiOffset.second;
                        auto newVelPhiLoc = std::make_pair(newVelPhiI, newVelPhiJ);
                        bufferVelPhi->setData(i, j, l->getVelocityPhi().sampleAt(newVelPhiLoc, ObsRef));

                        //Check if the cell contents can be advected
                        if (l->getObstacles(i, j)){
                            bufferTemp->setData(i, j, l->getTemperature(i, j));
                            bufferMoisture->setData(i, j, l->getCondensedWater(i, j));
                            bufferClouds->setData(i, j, l->getMoisture(i, j));
                            bufferParticulates->setData(i, j, l->getParticulates(i, j));
                            continue;
                        }
                        //Advect the cell contents
                        auto centreLoc = l->getObstacles().getWorldLoc(i, j);
                        double deltaCentreI = l->getVelocityTheta().sampleAt(centreLoc, ObsRef) * -dt / cellSizeX;
                        double deltaCentreJ = l->getVelocityPhi().sampleAt(centreLoc, ObsRef) * -dt / cellSizeY;
                        double newCentreI = (i + deltaCentreI) + centreOffset.first;
                        double newCentreJ = (j + deltaCentreJ) + centreOffset.second;
                        auto newCentreLoc = std::make_pair(newCentreI, newCentreJ);

                        double newTemp = l->getTemperature().sampleAt(newCentreLoc, ObsRef);
                        bufferTemp->setData(i, j, newTemp);
                        double newMois = l->getMoisture().sampleAt(newCentreLoc, ObsRef);
                        bufferMoisture->setData(i, j, newMois);
                        double newCloud = l->getCondensedWater().sampleAt(newCentreLoc, ObsRef);
                        bufferClouds->setData(i, j, newCloud);
                        double newPart = l->getParticulates().sampleAt(newCentreLoc, ObsRef);
                        bufferParticulates->setData(i, j, newPart);
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
            l->swapTemps(bufferTemp);
            l->swapMoistures(bufferMoisture);
            l->swapClouds(bufferClouds);
            l->swapParts(bufferParticulates);
        }

        /*template<typename T, typename V>
        void advectionEngine<T, V>::multiGridSolve(T& pressure, const T& divergence, int iteration){
            static T* fine_r[MAX_LAYER];
            static T* fine_e[MAX_LAYER];
            static T* coarse_r[MAX_LAYER];
            static T* coarse_e[MAX_LAYER];

            static float initialized = false;
            if (!initialized){
                for (int i = 0; i < MAX_LAYER; ++i){
                    fine_r[i] = nullptr;
                    fine_e[i] = nullptr;
                    coarse_r[i] = nullptr;
                    coarse_e[i] = nullptr;
                }
                initialized = true;
            }
            auto w = pressure.getWidth();
            auto s = pressure.getSize();
            if (! fine_r[iteration])
                fine_r [iteration] = new T(w, s);
            if (! fine_e[iteration])
                fine_e [iteration] = new T(w, s);
            if (! coarse_r[iteration])
                coarse_r [iteration] = new T(w / 2, s);
            if (! coarse_e[iteration])
                coarse_e [iteration] = new T(w / 2, s);

            if (! pressure.compareSizes(*fine_r[iteration]))
                fine_r[iteration]->resize(pressure);

            *fine_r[iteration] = 0;
            *fine_e[iteration] = 0;
            *coarse_r[iteration] = 0;
            *coarse_e[iteration] = 0;

            pressure.gaussSeidelSmoothBS(divergence, 4);

            //Compute residual
            pressure.laplacian(*fine_r[iteration]);
            fine_r[iteration]->mul(-1);
            fine_r[iteration]->add(divergence);

            if (pressure.getWidth() <= 2){
                coarse_e[iteration]->gaussSeidelSmoothBS(*coarse_r[iteration], 10, divergence.getWidth());
            }
            else{
                multiGridSolve(*coarse_e[iteration], *coarse_r[iteration], iteration + 1);
            }

            // Interpolate?
            coarse_e[iteration]->expand(*fine_e[iteration]);

            //Apply correction (p = p + e)
            pressure += *fine_e[iteration];

            pressure.gaussSeidelSmoothBS(divergence, 4);
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::solvePressureMatrix(std::shared_ptr<PWM::Model::airLayer<T, V> > &l){
           std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::solvePressureMatrix() for this template!" << std::endl;
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::solvePressureMatrix(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double> > &l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            int nCell = nX * nY;
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double zero = 0, negOne = -1, velSolid = 0;
            auto dt = getDt();
            double rho = PWM::Utils::altitudeAdjustedDensity(l->getHeight(), l->getPlanet());
            double dtRho = dt / rho;
            Eigen::SparseMatrix<double> A(nCell, nCell);
            A.reserve(Eigen::VectorXi::Constant(nCell, 5));

//            ATest = Eigen::MatrixXd::Zero(nCell, nCell);

            div.setZero(nCell);
            pressure.setZero(nCell);

            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        size_t cidx = PWM::Utils::convert2Dto1DUtil(nY, nX, i, j);
                        if (l->getObstacles(i, j)){
                            div[cidx] = 0;
                            continue;
                        }
                        int pijCoeff = 4;
                        double velBot = 0, velTop = 0, velRight = 0, velLeft = 0;
                        double coeffAbove = 0, coeffBelow = 0, coeffLeft = 0, coeffRight = 0;

                        //Check above the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i - 1, j)){
                            velTop = velSolid;
                            --pijCoeff;
                            coeffAbove = zero;
                        }
                        else{
                            velTop = l->getVelocityTheta(i, j);
                            coeffAbove = negOne;
                        }

                        //Check below the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i + 1, j)){
                            velBot = velSolid;
                            --pijCoeff;
                            coeffBelow = zero;
                        }
                        else{
                            velBot = l->getVelocityTheta(i + 1, j);
                            coeffBelow = negOne;
                        }

                        //Check left of the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i, j - 1)){
                            velLeft = velSolid;
                            --pijCoeff;
                            coeffLeft = zero;
                        }
                        else{
                            velLeft = l->getVelocityPhi(i, j);
                            coeffLeft = negOne;
                        }

                        //Check right of the cell for obstacles and adjust formula as needed
                        if (l->getObstacles(i, j + 1)){
                            velRight = velSolid;
                            --pijCoeff;
                            coeffRight = zero;
                        }
                        else{
                            velRight = l->getVelocityPhi(i, j + 1);
                            coeffRight = negOne;
                        }

                        double dive = (velBot - velTop) / cellSizeY + (velRight - velLeft) / cellSizeX;//Calculate divergence
                        double adjustedDive = (dive / dtRho) * (cellSizeX * cellSizeY);//Isolate coefficients to be either 4/3/2/1 for ij or -1 (or 0) for the neighbours.
                        size_t aIdx, bIdx, lIdx, rIdx;
                        aIdx = Utils::convert2Dto1DUtil(nY, nX, i - 1, j);
                        bIdx = Utils::convert2Dto1DUtil(nY, nX, i + 1, j);
                        lIdx = Utils::convert2Dto1DUtil(nY, nX, i, j - 1);
                        rIdx = Utils::convert2Dto1DUtil(nY, nX, i, j + 1);

                        //Set b value
                        div[cidx] = -adjustedDive;
                        //Set values in A (fill in relevant values in row ij)
                        A.insert(cidx, aIdx) = coeffAbove;
                        A.insert(cidx, bIdx) = coeffBelow;
                        A.insert(cidx, lIdx) = coeffLeft;
                        A.insert(cidx, rIdx) = coeffRight;
                        A.insert(cidx, cidx) = pijCoeff;

//                        ATest(cidx, aIdx) = coeffAbove;
//                        ATest(cidx, bIdx) = coeffBelow;
//                        ATest(cidx, lIdx) = coeffLeft;
//                        ATest(cidx, rIdx) = coeffRight;
//                        ATest(cidx, cidx) = pijCoeff;
                    }
                }
            #pragma omp barrier

//            Eigen::LDLT<Eigen::MatrixXd> ALLT(ATest);

//            Eigen::MatrixXd denseA(A);
//            Eigen::SparseMatrix<double> sparseA = ATest.sparseView();

//            if (!ATest.isApprox(denseA))
//                std::cerr << "A and ATest are not the same!" << std::endl;

//            if (!ALLT.isPositive()){
//                std::cerr << "Eigen thinks A is not positive!" << std::endl;
//                std::cerr << "A has determinant " << ATest.determinant() << std::endl;
//            }

//            Eigen::EigenSolver<Eigen::MatrixXd> es(ATest, false);
//            std::cerr << "A has eigenvalues:\n" << es.eigenvalues() << std::endl;

//            if (!ATest.isApprox(ATest.transpose())){
//                std::cerr << "A is not symmetric!" << std::endl;
//                std::cerr << "A:\n" << ATest << std::endl;
//                std::cerr << "\nA^T:\n" << ATest.transpose() << std::endl;
//            }
//            if (ALLT.info() == Eigen::NumericalIssue)
//                std::cerr << "A is not positive definite!" << std::endl;

//            Eigen::VectorXd x(nCell);
//            auto cphhSolve = ATest.colPivHouseholderQr();
//            x = cphhSolve.solve(div);
//            if (!(div.isApprox(A * x)))
//                std::cerr << "colPivHouseholderQr does not solve Ap = div correctly." << std::endl;

            A.makeCompressed();
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> solver;
            solver.setTolerance(solverTolerance);
            solver.setMaxIterations(solverIterations);
            solver.compute(A);
            if (solver.info() != Eigen::Success)
                std::cerr << "Eigen solver decomposition failed!" << std::endl;

            #pragma omp parallel for
                for (int i = 0; i < nCell; ++i){
                    PrevPressure[i] = l->getPressure(i);
                }
            #pragma omp barrier

            pressure = solver.solveWithGuess(div, PrevPressure);
            if (solver.info() != Eigen::Success)
                std::cerr << "Eigen solver solve failed after " << solver.iterations() << " steps!" << std::endl;
            if (!(solver.error() < solverTolerance))
                std::cerr << "Solver error: " << solver.error() << std::endl;
//            if (x.isApprox(pressure)){
//                std::cerr << "ColPivHouseholder and CG give same results! Hooray!" << std::endl;
//            }
//            if (div.isApprox(A * pressure)){
//                std::cerr << "ConjugateGradient correctly solves Ap = div in " << solver.iterations() << " iterations." << std::endl;
//            }

            #pragma omp parallel for
                for (int i = 0; i < nCell; ++i){
                    l->setPressure(i, pressure[i]);
                }
            #pragma omp barrier
            #ifdef DEBUG
                std::cerr << "Solver used " << solver.iterations() << " iterations to reach an error of " << solver.error() << "." << std::endl;
                std::cerr << "Pressure has a range of " << l->getPressure().range() << ", a mean of " << l->getPressure().mean() << ", and a standard deviation of " << l->getPressure().stdDev() << std::endl;
            #endif
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::smoothPressure(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::smoothPressure() for this template!" << std::endl;
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::smoothPressure(){
            startComputation();
            for (auto l : airLayers){
                calcDiv(l);
                l->getPressure().gaussSeidelSmoothBS(*divergence, 4);
            }
            endComputation();
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<T, V>>& l, V scaleFactor){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<T, V>>& l, V scaleFactor) for this template!" << std::endl;
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::square2DArray<double>, double>>& l, double scaleFactor){
            int cellNum = l->getObstacles().getWidth();
            auto dt = getDt();
            auto rho = PWM::Utils::altitudeAdjustedDensity(l->getHeight(), l->getPlanet());
            double cellSize = 1.0f / cellNum;
            #pragma omp parallel for
                for (int i = 0; i < cellNum; ++i){
                    for (int j = 0; j < cellNum; ++j){
                        double gradX = scaleFactor * (dt / rho) * (l->getPressure(i + 1, j) - l->getPressure(i - 1, j)) / cellSize;
                        double gradY = scaleFactor * (dt / rho) * (l->getPressure(i, j + 1) - l->getPressure(i, j - 1)) / cellSize;
                        bufferVelTheta->setData(i, j, (l->getVelocityTheta(i, j) - gradX));
                        bufferVelPhi->setData(i, j, (l->getVelocityPhi(i, j) - gradY));
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::solvePressureProjection(std::shared_ptr<PWM::Model::airLayer<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>>& l, double scaleFactor){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double velSolid = 0;
            auto dt = getDt();
            auto rho = PWM::Utils::altitudeAdjustedDensity(l->getHeight(), l->getPlanet());
            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
//                        if (i == 256 && j == 256){
//                            std::cerr << "Debug point reached" << std::endl;
//                        }
                        //Check if the central cell is solid
                        if (l->getObstacles(i, j)){
                            bufferVelTheta->setData(i, j, velSolid);
                            bufferVelPhi->setData(i, j, velSolid);
                            continue;
                        }
                        //If not, check if cell above is solid
                        if (l->getObstacles(i - 1, j)){
                            bufferVelTheta->setData(i, j, velSolid);
                        }
                        else{//Pressure update as per normal formula
                            double gradX = scaleFactor * (dt / rho) * (l->getPressure(i, j) - l->getPressure(i - 1, j)) / cellSizeX;
                            double newX = l->getVelocityTheta(i, j) - gradX;
                            bufferVelTheta->setData(i, j, newX);
                        }
                        //Also check if cell left is solid
                        if (l->getObstacles(i, j - 1)){
                            bufferVelPhi->setData(i, j, velSolid);
                        }
                        else{//Pressure update as per normal formula
                            double gradY = scaleFactor * (dt / rho) * (l->getPressure(i, j) - l->getPressure(i, j - 1)) / cellSizeY;
                            double newY = l->getVelocityPhi(i, j) - gradY;
                            bufferVelPhi->setData(i, j, newY);
                        }
                    }
                }
            #pragma omp barrier

            //swap
            l->swapVels(bufferVelTheta, bufferVelPhi);
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::correctDivergence(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::correctDivergence() for this template!" << std::endl;
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::correctDivergence(){
            startComputation();
            for (auto l : airLayers){
                calcDiv(l);
                // std::cout << "Divergence:\n" << divergence->print() << std::endl;
                l->getPressure().mul(0);
                multiGridSolve(l->getPressure(), *divergence);
                //std::cout << "Pressure:\n" << l->getPressure().print() << std::endl;
                solvePressureProjection(l);
            }
            endComputation();
        }

        template<>
        void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::correctDivergence(){
            startComputation();
            for (auto l : airLayers){
                solvePressureMatrix(l);
                solvePressureProjection(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
            }
            endComputation();
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::calcPressures(){
            startComputation();
            for (auto l : airLayers){
                solvePressureMatrix(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
            }
            endComputation();
        }

        template<typename T, typename V>
        void advectionEngine<T, V>::solvePressureProjection(V scaleFactor){
            startComputation();
            for (auto l : airLayers){
                solvePressureProjection(l, scaleFactor);
            }
            endComputation();
        }*/

        template<typename T, typename V>
        inline void advectionEngine<T, V>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define advectionEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::square2DArray<double>, double>::step_internal(){
            for(auto l : airLayers){
                // std::cout << "X velocities prior to advection:\n" << l->getVelocityTheta().print() << std::endl;
                advect(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
                applyCoriolis(l);
            }
        }

        template<>
        inline void advectionEngine<PWM::PWMDataStructure::flatStaggeredGrid<double>, double>::step_internal(){
            for(auto l : airLayers){
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
