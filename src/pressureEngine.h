#ifndef BFD_PRESSURE_ENGINE_H
#define BFD_PRESSURE_ENGINE_H

#define peEpsilon
#define solverTolerance 0.1
#define solverIterations 1000

#include "abstractEngine.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "flatStaggeredGrid.h"
#include "fluidLayer.h"
#include <omp.h>
#include <memory>

namespace BasicFluidDynamics{
    namespace Engine{
        template<typename T, typename V> class pressureEngine : public AbstractEngine{
            private:
                //Data structure to hold reference to each of the layers in the simulation
                std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>> fluidLayers;

                //Function that calculates the divergence of the layer
                void calcDiv(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);

                void multiGridSolve(T& pressure, const T& divergence, int iteration = 0);

                //Function that calculates the pressures needed to correct for nabla u = 0.
                void solvePressureMatrix(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);

                //Solve the spherical Poisson equation using the Intel Math Kernel Library
                void sphericalPoissonSolver(BasicFluidDynamics::Model::fluidLayer<T, V>& l);

                //Calculate and apply the divergence correction of velocities
                void solvePressureProjection(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l, V scaleFactor = 1);

                std::shared_ptr<T> bufferVelTheta;
                std::shared_ptr<T> bufferVelPhi;
                std::shared_ptr<T> divergence;

                Eigen::VectorXd div;
                Eigen::VectorXd pressure;
                Eigen::VectorXd PrevPressure;
            protected:
                void step_internal() override;
            public:
                pressureEngine(float dt, bool active);
                pressureEngine(const size_t width, const V WrldSizeSize, const float dt, bool active);
                pressureEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, const float dt, bool active);

                //Calculate divergence and pressure (without solving) to simulate uplift/downdrafts to some extent
                void calcPressures();

                //Calculate the velocity adjustments resulting from the pressure field
                void solvePressureProjection();

                //Function to add a fluid layer to the engine
                void addLayer(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l);
                const std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>& getFluidLayers() const;
        };

        template<typename T, typename V>
        inline pressureEngine<T, V>::pressureEngine(float dt, bool active) : AbstractEngine(dt, active){
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();

            bufferVelTheta = std::make_shared<T>();
            bufferVelPhi = std::make_shared<T>();
            divergence = std::make_shared<T>();
            div = Eigen::VectorXd();
            pressure = Eigen::VectorXd();
            PrevPressure = Eigen::VectorXd();
        }

        template<typename T, typename V>
        inline pressureEngine<T, V>::pressureEngine(const size_t width, const V WrldSize, const float dt, bool active) : AbstractEngine(dt, active){
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();
            bufferVelTheta = std::make_shared<T>(width, WrldSize);
            bufferVelPhi = std::make_shared<T>(width, WrldSize);
            divergence = std::make_shared<T>(width, WrldSize);
        }

        template<typename T, typename V>
        inline pressureEngine<T, V>::pressureEngine(const size_t width, const size_t height, const V xWrldSize, const V yWrldSize, const float dt, bool active) : AbstractEngine(dt, active){
            fluidLayers = std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>();
            bufferVelTheta = std::make_shared<T>(width, height, 0.0, 0.5, xWrldSize, yWrldSize);
            bufferVelPhi = std::make_shared<T>(width, height, 0.5, 0.0, xWrldSize, yWrldSize);
            divergence = std::make_shared<T>(width, height, 0.5, 0.5, xWrldSize, yWrldSize);

            div = Eigen::VectorXd(width * height);
            pressure = Eigen::VectorXd(width * height);
            PrevPressure = Eigen::VectorXd(width * height);
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::addLayer(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            fluidLayers.push_back(l);
        }

        template<typename T, typename V>
        inline const std::vector<std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>>& pressureEngine<T, V>::getFluidLayers() const{
            return fluidLayers;
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::calcDiv(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::calcDiv(BasicFluidDynamics::Model::fluidLayer<T, V>& l) for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::calcDiv(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>>& l){
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
            #pragma omp barrier
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::solvePressureMatrix(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V> > &l){
           std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::solvePressureMatrix() for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::solvePressureMatrix(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double> > &l){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            int nCell = nX * nY;
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double zero = 0, negOne = -1, velSolid = 0;
            auto dt = getDt();
            double rho = l->getDensity();
            double dtRho = dt / rho;
            Eigen::SparseMatrix<double> A(nCell, nCell);
            A.reserve(Eigen::VectorXi::Constant(nCell, 5));

//            ATest = Eigen::MatrixXd::Zero(nCell, nCell);

            div.setZero(nCell);
            pressure.setZero(nCell);

            std::vector<Eigen::Triplet<double>> trips;

//            #pragma omp parallel for
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        size_t cidx = BasicFluidDynamics::Utils::convert2Dto1DUtil(nY, nX, i, j);
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

                        //A.insert(cidx, aIdx) = coeffAbove;
                        //A.insert(cidx, bIdx) = coeffBelow;
                        //A.insert(cidx, lIdx) = coeffLeft;
                        //A.insert(cidx, rIdx) = coeffRight;
                        //A.insert(cidx, cidx) = pijCoeff;

                        trips.emplace_back(Eigen::Triplet<double>(cidx, aIdx, coeffAbove));
                        trips.emplace_back(Eigen::Triplet<double>(cidx, bIdx, coeffBelow));
                        trips.emplace_back(Eigen::Triplet<double>(cidx, lIdx, coeffLeft));
                        trips.emplace_back(Eigen::Triplet<double>(cidx, rIdx, coeffRight));
                        trips.emplace_back(Eigen::Triplet<double>(cidx, cidx, pijCoeff));

//                        ATest(cidx, aIdx) = coeffAbove;
//                        ATest(cidx, bIdx) = coeffBelow;
//                        ATest(cidx, lIdx) = coeffLeft;
//                        ATest(cidx, rIdx) = coeffRight;
//                        ATest(cidx, cidx) = pijCoeff;
                    }
                }
//             #pragma omp barrier

            A.setFromTriplets(trips.begin(), trips.end());

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
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper/*, Eigen::IncompleteCholesky<double, Eigen::Lower|Eigen::Upper>*/> solver;
            solver.setTolerance(solverTolerance);
            solver.setMaxIterations(solverIterations);
            solver.compute(A);
            if (solver.info() != Eigen::Success)
                std::cerr << "Eigen solver decomposition failed! Solver instead has " << solver.info() << std::endl;

            #pragma omp parallel for
                for (int i = 0; i < nCell; ++i){
                    PrevPressure[i] = l->getPressure(i);
                }
            #pragma omp barrier

            pressure = solver.solveWithGuess(div, PrevPressure);
            if (solver.info() != Eigen::Success){
                std::cerr << "Eigen solver solve failed!" << std::endl;
            }
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
        inline void pressureEngine<T, V>::solvePressureProjection(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l, V scaleFactor){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::solvePressureProjection(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<T, V>>& l, V scaleFactor) for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::solvePressureProjection(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>>& l, double scaleFactor){
            int nX = l->getObstacles().getX(), nY = l->getObstacles().getY();
            double cellSizeX = l->getObstacles().cellSizeX(), cellSizeY = l->getObstacles().cellSizeY();
            double velSolid = 0;
            auto dt = getDt();
            double rho = l->getDensity();
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
        inline void pressureEngine<T, V>::solvePressureProjection(){
            startComputation();
            for (auto l : fluidLayers){
                solvePressureProjection(l);
            }
            endComputation();
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::calcPressures(){
            startComputation();
            for (auto l : fluidLayers){
                solvePressureMatrix(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
            }
            endComputation();
        }

        template<typename T, typename V>
        inline void pressureEngine<T, V>::step_internal(){
            std::cout << "Error! Not implemented for general template, use a specialized template or define pressureEngine::step_internal() for this template!" << std::endl;
        }

        template<>
        inline void pressureEngine<BasicFluidDynamics::Data::flatStaggeredGrid<double>, double>::step_internal(){
            for (auto l : fluidLayers){
                solvePressureMatrix(l);
                solvePressureProjection(l);
                // std::cout << "X velocities after advection:\n" << l->getVelocityTheta().print() << std::endl;
            }
        }
    }
}
#endif // PRESSUREENGINE_H
