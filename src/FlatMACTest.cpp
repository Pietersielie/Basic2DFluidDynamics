/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
 */

#include "advectionEngine.h"
#include "fluidLayer.h"
#include <filesystem>
#include "flatStaggeredGrid.h"
#include <iostream>
#include "pressureEngine.h"
#include "vizUtils.h"

typedef double valType;
typedef BasicFluidDynamics::Data::flatStaggeredGrid<valType> dsType;

void testConstructors(){

}

void testComparators(){

}

void testGetters(){

}

void testSetters(){

}

void testMiscellany(const size_t nX, const size_t nY, const int steps, const int frac){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && ((dir_Entry.path().extension().string() == ".png") || (dir_Entry.path().extension().string() == ".ppm")))
            std::filesystem::remove(dir_Entry.path());
    }

    auto l1 = std::make_shared<BasicFluidDynamics::Model::fluidLayer<dsType, valType>>(BasicFluidDynamics::Model::fluidLayer<dsType, valType>(nX, nY, 10000, 2500));
    l1->getObstacles().copy(0);
    for (int i = 0; i < nX; ++i){//For width
        l1->setObstacles(0, i, 1);
        l1->setObstacles(nY - 1, i, 1);
    }
    auto radius = nY / 6;
    std::pair<size_t, size_t> hillCentre = std::make_pair(nY / 2, nX / 5);
    int iMin = hillCentre.first - radius, iMax = hillCentre.first + radius;
    int jMin = hillCentre.second - radius, jMax = hillCentre.second + radius;
    for (int i = iMin; i < iMax; ++i){
        for (int j = jMin; j < jMax; ++j){
            valType dist = BasicFluidDynamics::Utils::calcCartesianDistance<valType, int, int>(hillCentre, std::make_pair(i, j));
            if (dist <= radius){
//                l1->setObstacles(i, j, 1.);
//                l1->setObstacles(i, j + (nX / 2), 1.);
                l1->setObstacles(i /*+ (nY / 3.2)*/, j + (nX / 4), 1.);
                l1->setTemperature(i, j + (nX / 3), 0.);
//                l1->setObstacles(i - (nY / 3.2), j + (nX / 4), 1.);
            }
        }
    }
    std::stringstream fO;
    fO << "../ImageOutput/Layer_Obstacles.ppm";
    BasicFluidDynamics::Utils::writeTerrElevImage(fO.str(), l1->getObstacles());
    l1->getVelocityPhi().copy(80);
    l1->getVelocityTheta().copy(0);
    int hotWidth = nY / 10;
    int startj = nY / 2, startk = 0;
    double highT = 1000.0;
    for (int i = -hotWidth; i < hotWidth; ++i){
        for (int j = -hotWidth; j < hotWidth; ++j){
            l1->setTemperature(i + startj, j + startk, highT);
        }
    }

    auto x = BasicFluidDynamics::Engine::advectionEngine<dsType, valType>(nX, nY, 10000, 2500, 0.05f, true, false);
    auto y = BasicFluidDynamics::Engine::pressureEngine<dsType, valType>(nX, nY, 10000, 2500, 0.05f, true);
    x.addLayer(l1);
    y.addLayer(l1);
    double previousStep = 0;
    for (int i = 0; i <= steps; ++i){
//        if (i == steps / 2){
//            for (int j = 0; j < nY; ++j){
//                l1->setObstacles(j, 0, 1);
//                l1->setObstacles(j, nX - 1, 1);
//            }
//            std::stringstream fO1;
//            fO1 << "../ImageOutput/FlatMACTest/Layer_Obstacles_Step_" << i << ".ppm";
//            BasicFluidDynamics::Utils::writeTerrElevImage(fO1.str(), l1->getObstacles());
//        }
        if (i % frac == 0){
            std::stringstream fT, fVY, fVX, fP;
            fT << "../ImageOutput/FlatMACTest/Layer_Temperature_Step_" << i + 1 << ".ppm";
            fVX << "../ImageOutput/FlatMACTest/Layer_VelX_Step_" << i + 1 << ".ppm";
            fVY << "../ImageOutput/FlatMACTest/Layer_VelY_Step_" << i + 1 << ".ppm";
            fP << "../ImageOutput/FlatMACTest/Layer_Pressure_Step_" << i + 1 << ".ppm";
            BasicFluidDynamics::Utils::writeTempImage(fT.str(), l1->getColour());
            BasicFluidDynamics::Utils::writeVelImage(fVX.str(), l1->getVelocityTheta());
            BasicFluidDynamics::Utils::writeVelImage(fVY.str(), l1->getVelocityPhi());
//            BasicFluidDynamics::Utils::writeTerrElevImage(fP.str(), l1->getPressure());
        }
        previousStep = x.getRunTimePassed();
        x.step();
        y.step();
        double currentStep = x.getRunTimePassed();
        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;
    }
}

int main(int argc, char** argv){
    std::cout << "\nTesting advection with \033[1;34m'MAC grid'\033[0m:" << std::endl;

    std::cout << "\n\033[1;37mTesting miscellany:\033[0m" << std::endl;
    testMiscellany(1024, 256, 6400, 10);

    std::cout << "\nTesting advection with \033[1;34m'MAC grid'\033[0m complete." << std::endl;
}
