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

#include <opencv2/opencv.hpp>

typedef double valType;
typedef BasicFluidDynamics::Data::flatStaggeredGrid<valType> dsType;

constexpr int xBuffer = 30;
constexpr valType smokeColour = 800;
constexpr valType baseColour = 200;
constexpr valType cellSize = 20;

std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<dsType, valType>> constructAirLayer(const std::string file){
    cv::Mat Image;
    cv::threshold(cv::imread(file, cv::IMREAD_GRAYSCALE), Image, 127, 255, cv::ThresholdTypes::THRESH_BINARY_INV);
    int imWidth = Image.cols;
    int imHeight = Image.rows;

    int nX = imWidth + 4 * xBuffer;
    int nY = imHeight + 2;
    auto res = std::make_shared<BasicFluidDynamics::Model::fluidLayer<dsType, valType>>(BasicFluidDynamics::Model::fluidLayer<dsType, valType>(nX, nY, cellSize * nX, cellSize * nY));
    res->getObstacles().copy(0);
    res->getColour().copy(baseColour);
    //Set the top and bottom of wind tunnel
    for (int i = 0; i < nX; ++i){//For width
        res->setObstacles(0, i, 1);
        res->setObstacles(nY - 1, i, 1);
    }
    //Set the obstacles/shapes in the wind tunnel
    for (int i = 0; i < imHeight; ++i){
        for (int j = 0; j < imWidth; ++j){
            uchar val = Image.at<uchar>(i, j);
            res->setObstacles(i + 1, j + (2 * xBuffer), val > 0 ? 1 : 0);
        }
    }
    return res;
}

//Set conditions at start and end of windtunnel (i.e., colour and speed)
void setWindTunnel(std::shared_ptr<BasicFluidDynamics::Model::fluidLayer<dsType, valType>>& layer, const valType speed){
    int nX = layer->getObstacles().getX(), nY = layer->getObstacles().getY();
    int bufferWidth = xBuffer / 2;
    for (int i = 1; i < nY - 1; ++i){
        for (int j = 0; j < bufferWidth; ++j){
            layer->setVelocityX(i, j, speed);
            layer->setVelocityX(i, nX - j - 1, speed);
            layer->setVelocityY(i, j, 0.0);
            layer->setVelocityY(i, nX - j - 1, 0.0);
            layer->setColour(i, j, baseColour);
        }
    }

    int hotWidth = xBuffer;
    int startj = nY / 2, startk = 0;
    for (int i = -hotWidth; i < hotWidth; ++i)
        for (int j = 0; j < xBuffer; ++j)
            layer->setColour(i + startj, j + startk, smokeColour);
}

void testMiscellany(const std::string shapeFile, const valType windSpeed, const float simTime, const size_t frac){
    auto currDir = std::filesystem::current_path();
    currDir.remove_filename().remove_filename().concat("ImageOutput/");
    for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
        if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && ((dir_Entry.path().extension().string() == ".png") || (dir_Entry.path().extension().string() == ".ppm")))
            std::filesystem::remove(dir_Entry.path());
    }

    float dT = (windSpeed * 1.4) / cellSize;

    auto l1 = constructAirLayer(shapeFile);
    setWindTunnel(l1, windSpeed);
//    l1->getVelocityX().copy(windSpeed);
    std::stringstream fT, fVY, fVX, fP;
    fT << "../ImageOutput/Init_Colour.png";
    fVX << "../ImageOutput/Init_VelX.png";
    fVY << "../ImageOutput/Init_VelY.png";
    fP << "../ImageOutput/Init_Pressure.png";
    BasicFluidDynamics::Utils::writeTempImage(fT.str(), l1->getColour(), l1->getObstacles(), 0);
    BasicFluidDynamics::Utils::writeVelImage(fVX.str(), l1->getVelocityY(), l1->getObstacles(), 0);
    BasicFluidDynamics::Utils::writeVelImage(fVY.str(), l1->getVelocityX(), l1->getObstacles(), 0);
    BasicFluidDynamics::Utils::writePresImage(fP.str(), l1->getPressure(), l1->getObstacles(), 0);

    int nX = l1->getObstacles().getX(), nY = l1->getObstacles().getY();
    valType xS = l1->getObstacles().getXWorldLength(), yS = l1->getObstacles().getYWorldLength();
    auto x = BasicFluidDynamics::Engine::advectionEngine<dsType, valType>(nX, nY, xS, yS, dT, true, false);
    auto y = BasicFluidDynamics::Engine::pressureEngine<dsType, valType>(nX, nY, xS, yS, dT, true);
    x.addLayer(l1);
    y.addLayer(l1);

    int steps = std::ceil(simTime / dT);

    std::cout << "Testing shape " << shapeFile << " in a wind tunnel at speed " << windSpeed << " metres per second, for " << simTime << " seconds and a timestep of " << dT << " seconds." << std::endl;
    std::cout << "The grid is " << nY << " cells high and " << nX << " cells wide." << std::endl;
    double previousStep = 0;
    for (int i = 0; i <= steps; ++i){
        if (i % frac == 0){
            std::stringstream fT, fVY, fVX, fP;
            fT << "../ImageOutput/Colour_Step_" << i + 1 << ".png";
            fVX << "../ImageOutput/VelX_Step_" << i + 1 << ".png";
            fVY << "../ImageOutput/VelY_Step_" << i + 1 << ".png";
            fP << "../ImageOutput/Pressure_Step_" << i + 1 << ".png";
            BasicFluidDynamics::Utils::writeTempImage(fT.str(), l1->getColour(), l1->getObstacles(), xBuffer);
            BasicFluidDynamics::Utils::writeVelImage(fVX.str(), l1->getVelocityY(), l1->getObstacles(), xBuffer);
            BasicFluidDynamics::Utils::writeVelImage(fVY.str(), l1->getVelocityX(), l1->getObstacles(), xBuffer);
            BasicFluidDynamics::Utils::writePresImage(fP.str(), l1->getPressure(), l1->getObstacles(), xBuffer);
        }
        previousStep = x.getRunTimePassed();
        x.step();
        y.step();
        double currentStep = x.getRunTimePassed();
        std::cout << "Step " << i + 1 << " took " << currentStep - previousStep << " seconds." << std::endl;
        setWindTunnel(l1, windSpeed);
    }
}

int main(int argc, char** argv){
    std::string fil = "../resources/Mask1.png";
    valType speed = 10;
    float time = 300;
    int interval = 5;
    if (argc > 1){
        fil = argv[1];
    }
    if (argc > 2){
        speed = std::stod(argv[2]);
    }
    if (argc > 3){
        time = std::stof(argv[3]);
    }
    if (argc > 4){
        interval = std::stoi(argv[4]);
    }
    if (argc == 1 || argc > 5){
        std::cout << "Input is available in the form of CLI parameters as follows:\n\n    [shapeFileName] [WindSpeed (m.s-1)] [simulationTime (s)] [interval (integer)]\n\n"
                  << "Default values are:"
                  << "\n  Shape file: " << fil
                  << "\n  Wind speed: " << speed
                  << "\n  Simulation time: " << time
                  << "\n  Interval: " << interval << std::endl;
    }

    testMiscellany(fil, speed, time, interval);
}
