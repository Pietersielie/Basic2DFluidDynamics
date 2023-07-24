/**
 * Template for testing some condition with formatting included:
 * if (condition)
 * 	   std::cout << "\033[1;32m<function> works as expected.\033[0m" << std::endl;
 * else
 *	   std::cout << "\033[1;31mError! <function> does not perform as expected!\033[0m" << std::endl;
 */

#include "advectionEngine.h"
#include "airLayer.h"
#include <filesystem>
#include "flatStaggeredGrid.h"
#include <iostream>
#include "planet.h"
#include "pressureEngine.h"
#include "vizUtils.h"

typedef double valType;
typedef PWM::PWMDataStructure::flatStaggeredGrid<valType> dsType;

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
    currDir.remove_filename().remove_filename().concat("ImageOutput/FlatMACTest");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && ((dir_Entry.path().extension().string() == ".png") || (dir_Entry.path().extension().string() == ".ppm")))
            std::filesystem::remove(dir_Entry.path());
    }

    auto planetos = std::make_shared<PWM::Model::planet>(PWM::Model::planet("../resources/PlanetEarth.json"));
    auto l1 = std::make_shared<PWM::Model::airLayer<dsType, valType>>(PWM::Model::airLayer<dsType, valType>(planetos, 100, 100, nX, nY, 10000, 2500));
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
            valType dist = PWM::Utils::calcCartesianDistance<valType, int, int>(hillCentre, std::make_pair(i, j));
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
    fO << "../ImageOutput/FlatMACTest/Layer_Obstacles.ppm";
    PWM::Utils::writeTerrElevImage(fO.str(), l1->getObstacles());
//    l1->getVelocityPhi().randomInit(-20, 20);
//    l1->getVelocityTheta().randomInit(-20, 20);
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

    auto x = PWM::Engine::advectionEngine<dsType, valType>(nX, nY, 10000, 2500, 0.05f, true, false);
    auto y = PWM::Engine::pressureEngine<dsType, valType>(nX, nY, 10000, 2500, 0.05f, true);
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
//            PWM::Utils::writeTerrElevImage(fO1.str(), l1->getObstacles());
//        }
        if (i % frac == 0){
            std::stringstream fT, fVY, fVX, fP;
            fT << "../ImageOutput/FlatMACTest/Layer_Temperature_Step_" << i + 1 << ".ppm";
            fVX << "../ImageOutput/FlatMACTest/Layer_VelX_Step_" << i + 1 << ".ppm";
            fVY << "../ImageOutput/FlatMACTest/Layer_VelY_Step_" << i + 1 << ".ppm";
            fP << "../ImageOutput/FlatMACTest/Layer_Pressure_Step_" << i + 1 << ".ppm";
            PWM::Utils::writeTempImage(fT.str(), l1->getTemperature());
            PWM::Utils::writeVelImage(fVX.str(), l1->getVelocityTheta());
            PWM::Utils::writeVelImage(fVY.str(), l1->getVelocityPhi());
//            PWM::Utils::writeTerrElevImage(fP.str(), l1->getPressure());
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
