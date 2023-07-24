#ifndef VIZUTILS_H
#define VIZUTILS_H

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include "flatStaggeredGrid.h"
#include "planet.h"
#include "square2DArray.h"
#include "terrain_structure.hpp"
#include <tuple>
#include <vector>
#include "world.h"

namespace PWM {
    namespace Utils {
        inline std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> cloudColourMap = {
            std::make_tuple(0, 0, 0), //EMPTY
            std::make_tuple(0, 255, 0),//STRATUS,
            std::make_tuple(0, 153, 51),//ALTOSTRATUS,
            std::make_tuple(0, 204, 255),//CIRROSTRATUS,
            std::make_tuple(153, 255, 51),//NIMBOSTRATUS,
            std::make_tuple(0, 0, 255),//CIRRUS,
            std::make_tuple(255, 0, 0),//CUMULUS,
            std::make_tuple(255, 255, 0),//STRATOCUMULUS,
            std::make_tuple(204, 0, 102),//ALTOCUMULUS,
            std::make_tuple(255, 0, 255),//CIRROCUMULUS,
            std::make_tuple(255, 153, 0),//CUMULONIMBUS,
            std::make_tuple(255, 255, 255)//CTYPEEND
        };

        inline int createGifImage(std::string fileName, std::string filePattern){
            try {
                std::string s = "convert -delay 30 -loop 0 " + filePattern + " " + fileName + ".gif";
                char c[s.size() + 1];
                std::strcpy(c, s.c_str());
                return std::system(c);
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write GIF image to file " << fileName << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeAirGifImages(std::string directory, const std::shared_ptr<PWM::Model::world<PWM::PWMDataStructure::square2DArray<T>, PWM::PWMDataStructure::square2DArray<std::string>, T, std::string>>& wld){
            for (int j = 0; j < wld->getAirLayers().size(); ++j){
                std::stringstream fT, fP, fV1, fV2, fCV;
                fT << "../" << directory << "/Layer_" << j << "_\\(" << wld->getAirLayer(j)->getHeight() << "m\\)_Temperature";
                fP << "../" << directory << "/Layer_" << j << "_\\(" << wld->getAirLayer(j)->getHeight() << "m\\)_Moisture";
                fV1 << "../" << directory << "/Layer_" << j << "_\\(" << wld->getAirLayer(j)->getHeight() << "m\\)_VelocityX";
                fV2 << "../" << directory << "/Layer_" << j << "_\\(" << wld->getAirLayer(j)->getHeight() << "m\\)_VelocityY";
                fCV << "../" << directory << "/ConvectionLayer_" << j << "_Vertical_Velocity";
                createGifImage(fT.str(), fT.str() + "_Step_*.ppm");
                createGifImage(fP.str(), fP.str() + "_Step_*.ppm");
                createGifImage(fV1.str(), fV1.str() + "_Step_*.ppm");
                createGifImage(fV2.str(), fV2.str() + "_Step_*.ppm");
                if (j < wld->getAirLayers().size() - 1) createGifImage(fCV.str(), fCV.str() + "_Step_*.ppm");;
            }

            auto currDir = std::filesystem::current_path().remove_filename().concat(directory);
            for (auto dir_Entry : std::filesystem::directory_iterator(currDir)){
                if (!dir_Entry.path().empty() && dir_Entry.path().has_filename() && dir_Entry.path().extension().string() == ".ppm")
                    std::filesystem::remove(dir_Entry.path());
            }
            return 0;
        }

        template<typename T>
        inline int writeMoisImage(std::string file, std::shared_ptr<PWM::Model::planet> planetos, const PWM::PWMDataStructure::square2DArray<T>& data, std::pair<T, T>& minMax){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = minMax;
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    float moisNormed = (std::abs(spread) < 1 ? 0 : ((data.getData(i) - mM.first) / spread) * 255);
                    int r = 255 - (255 - moisNormed);
                    int g = 0;
                    int b = 255 - moisNormed;
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeMoisImage(std::string file, std::shared_ptr<PWM::Model::planet> planetos, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(0.2, 1.0);
                T max = 2.0;
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    float moisNormed, mois = data.getData(i);
                    if (mois <= 0){
                        r = g = b = 0;
                    }
                    else if (mois < mM.first){
                        moisNormed = ((mois) / mM.first) * 255;
                        r = 0;
                        g = 0;
                        b = 255 - (255 - moisNormed);
                    }
                    else if (mois < mM.second){
                        moisNormed = ((mois - mM.first) / spread) * 255;
                        r = 0;
                        g = 255 - (255 - moisNormed);
                        b = 255 - moisNormed;
                    }
                    else if (mois < max){
                        moisNormed = (((mois - mM.second) / (max - mM.second)) * 255);
                        r = 255 - (255 - moisNormed);
                        g = 255 - moisNormed;
                        b = 0;
                    }
                    else{
                        r = g = b = 255;
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeMoisImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(0.1, 1.0);
                T max = 2.0;
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    float moisNormed, mois = data.getData(i);
                    if (mois <= 0){
                        r = g = b = 0;
                    }
                    else if (mois < mM.first){
                        moisNormed = ((mois) / mM.first) * 255;
                        r = 0;
                        g = 0;
                        b = 255 - (255 - moisNormed);
                    }
                    else if (mois < mM.second){
                        moisNormed = ((mois - mM.first) / spread) * 255;
                        r = 0;
                        g = 255 - (255 - moisNormed);
                        b = 255 - moisNormed;
                    }
                    else if (mois < max){
                        moisNormed = (((mois - mM.second) / (max - mM.second)) * 255);
                        r = 255 - (255 - moisNormed);
                        g = 255 - moisNormed;
                        b = 0;
                    }
                    else{
                        r = g = b = 255;
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTempImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(240, 350);
                T spread = mM.second - mM.first;
                T max = 1200;

                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;

                    T t = data.getData(i);

                    // JG - change of viz to help build intuition
                    // colour band approach
                    if(t < 213.0f) // < -60C
                    {
                        r = 0; g = 0; b = 0; // black
                    }
                    else if(t < 263.15f) // [-60C, -10C]
                    {
                        r = 0; g = 0; b = 255; // blue
                    }
                    else if(t < 283.15f) // [-10C,10C]
                    {
                        r = 255; g = 0; b = 255; // purple
                    }
                    else if(t < 373.15)
                    {
                        r = 255; g = (t-283.15) / 90.0f; b = (t-283.15) / 90.0f; // red shading to white
                    }
                    else
                    {
                        r = 255; g = 255; b = 255; // white
                    }

                    /*
                            float tNormed, t = data.getData(i);
                            if (t < mM.first){
                                tNormed = (std::abs(spread) < 1 ? 0 : ((t) / mM.first) * 255);
                                r = 0;
                                g = 0;
                                b = 255 - (255 - tNormed);
                            }
                            else if (t < mM.second){
                                tNormed = (std::abs(spread) < 1 ? 0 : ((t - mM.first) / spread) * 255);
                                r = 255 - (255 - tNormed);
                                g = 0;
                                b = 255 - tNormed;
                            }
                            else if (t < max){
                                tNormed = (((t - mM.second) / (max - mM.second)) * 255);
                                r = 255;
                                g = 255 - (255 - tNormed);
                                b = 255 - (255 - tNormed);
                            }
                            else{
                                r = g = b = 255;
                            } */
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTempImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(240, 350);
                T spread = mM.second - mM.first;
                T max = 1200;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    T t = data.getData(i);

                    // JG - change of viz to help build intuition
                    // colour band approach
                    if(t < 213.0f) // < -60C
                    {
                        r = 0; g = 0; b = 0; // black
                    }
                    else if(t < 263.15f) // [-60C, -10C]
                    {
                        r = 0; g = 0; b = 255; // blue
                    }
                    else if(t < 283.15f) // [-10C,10C]
                    {
                        r = 255; g = 0; b = 255; // purple
                    }
                    else if(t < 373.15)
                    {
                        r = 255; g = (t-283.15) / 90.0f; b = (t-283.15) / 90.0f; // red shading to white
                    }
                    else
                    {
                        r = 255; g = 255; b = 255; // white
                    }
                    /*
                            int tNormed;
                            T t = data.getData(i);
                            if (t < mM.first){
                                tNormed = (std::abs(spread) < 1 ? 0 : ((t) / mM.first) * 255);
                                r = 0;
                                g = 0;
                                b = 255 - (255 - tNormed);
                            }
                            else if (t < mM.second){
                                tNormed = (std::abs(spread) < 1 ? 0 : ((t - mM.first) / spread) * 255);
                                r = 255 - (255 - tNormed);
                                g = 0;
                                b = 255 - tNormed;
                            }
                            else if (t < max){
                                tNormed = (((t - mM.second) / (max - mM.second)) * 255);
                                r = 255;
                                g = 255 - (255 - tNormed);
                                b = 255 - (255 - tNormed);
                            }
                            else{
                                r = g = b = 255;
                            }
         */
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeAshImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(240, 350);
                T spread = mM.second - mM.first;
                T max = 1;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    T rho = data.getData(i);
                    if (rho > max)
                    {
                        r,g=0;
                        b = 255;
                    }
                    else
                    {
                        r=255 - 255*rho;
                        g=255 - 255*rho;
                        b=255;
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeAshImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                auto mM = std::make_pair<T, T>(240, 350);
                T spread = mM.second - mM.first;
                T max = 1;
                for (int i = 0; i < data.size(); ++i){
                    int r, g, b;
                    T rho = data.getData(i);
                    if (rho > max)
                    {
                        r,g=0;
                        b = 255;
                    }
                    else
                    {
                        r=255 - 255*rho;
                        g=255 - 255*rho;
                        b=255;
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeVelImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                T max = 40;
                int r, g, b;
                for (int i = 0; i < data.size(); ++i){
                    T t = data.getData(i);
                    int tNormed = ((std::abs(t) / max) * 255);
                    r = 0;
                    if (PWM::Utils::getSign(t) > 0){
                        if (t < max){
                            g = 255 - std::min(255 - tNormed, 255);
                            b = 0;
                        }
                        else{
                            g = 255;
                            b = 0;
                        }
                    }
                    else{
                        if (t > -max){
                            g = 0;
                            b = 255 - std::min(255 - tNormed, 255);

                        }
                        else{
                            g = 0;
                            b = 255;
                        }
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeVelImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                T max = 140;
                int r, g, b;
                for (int i = 0; i < data.size(); ++i){
                    int tNormed = ((std::abs(data.getData(i)) / max) * 255);
                    r = 0;
                    if (PWM::Utils::getSign(data.getData(i)) > 0){
                        g = 255 - std::min(255 - tNormed, 255);
                        b = 0;
                    }
                    else{
                        g = 0;
                        b = 255 - std::min(255 - tNormed, 255);
                    }
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeCloudImages(std::string outputDirectory, int i, int j, float height, const PWM::PWMDataStructure::square2DArray<T>& mask){
            std::stringstream fC;
            fC << outputDirectory << "/Layer_" << j << "_(" << height << "m)_Clouds_Step_" << i << ".ppm";
            std::ofstream fil;
            try {
                //                std::cerr << "Writing layer " << j << " to " << fC.str() << "." << std::endl;
                fil.open(fC.str());
                //Header
                fil << "P3\n" << mask.getWidth() << " " << mask.getWidth() << "\n255" << std::endl;

                //Image data
                for (int x = 0; x < mask.getWidth(); ++x){
                    for (int y = 0; y < mask.getWidth(); ++y){
                        auto cloudType = static_cast<int>(mask.getData(x, y));
                        auto colour = cloudColourMap[cloudType];
                        fil << std::get<0>(colour) << " " << std::get<1>(colour) << " " << std::get<2>(colour) << " ";
                    }
                }
                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << fC.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeCloudImages(std::string outputDirectory, int i, int j, float height, const PWM::PWMDataStructure::flatStaggeredGrid<T>& mask){
            std::stringstream fC;
            fC << outputDirectory << "/Layer_" << j << "_(" << height << "m)_Clouds_Step_" << i << ".ppm";
            std::ofstream fil;
            size_t nX = mask.getX(), nY = mask.getY();
            try {
                //                std::cerr << "Writing layer " << j << " to " << fC.str() << "." << std::endl;
                fil.open(fC.str());
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                for (int x = 0; x < nY; ++x){
                    for (int y = 0; y < nX; ++y){
                        auto cloudType = static_cast<int>(mask.getData(x, y));
                        auto colour = cloudColourMap[cloudType];
                        fil << std::get<0>(colour) << " " << std::get<1>(colour) << " " << std::get<2>(colour) << " ";
                    }
                }
                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << fC.str() << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTerrElevImage(std::string file, const PWM::PWMDataStructure::square2DArray<T>& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.getWidth() << " " << data.getWidth() << "\n255" << std::endl;

                //Image data
                auto mM = data.minmax();
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.size(); ++i){
                    float tNormed = ((data.getData(i) - mM.first) / spread) * 255;
                    int r = 255 - (255 - tNormed);
                    int g = 255 - (255 - tNormed);
                    int b = 255 - (255 - tNormed);
                    fil << r << " " << g << " " << b << " ";
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTerrElevImage(std::string file, const PWM::PWMDataStructure::flatStaggeredGrid<T>& data){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << nX << " " << nY << "\n255" << std::endl;

                //Image data
                auto mM = data.minmax();
                T spread = mM.second - mM.first;
                for (int i = 0; i < nY; ++i){
                    for (int j = 0; j < nX; ++j){
                        float tNormed = ((data.getData(i, j) - mM.first) / spread) * 255;
                        int r = 255 - (255 - tNormed);
                        int g = 255 - (255 - tNormed);
                        int b = 255 - (255 - tNormed);
                        fil << r << " " << g << " " << b << " ";
                    }
                    fil << '\n';
                }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int writeTerrElevImage(std::string file, const terrain_structure& data){
            std::ofstream fil;
            try {
                fil.open(file);
                //Header
                fil << "P3\n" << data.field_size << " " << data.field_size << "\n255" << std::endl;

                //Image data
                T min = std::numeric_limits<T>::max(), max = std::numeric_limits<T>::min();
                for (int i = 0; i < data.field_size; ++i)
                    for (int j = 0; j < data.field_size; ++j){
                        double val = data.height_field[i][j];
                        if (val < min)
                            min = val;
                        if (val > max)
                            max = val;
                    }
                auto mM = std::make_pair(min, max);
                T spread = mM.second - mM.first;
                for (int i = 0; i < data.field_size; ++i)
                    for (int j = 0; j < data.field_size; ++j){
                        float tNormed = (std::abs(spread) < 1 ? 0 : ((data.height_field[i][j] - mM.first) / spread) * 255);
                        int r = 255 - (255 - tNormed);
                        int g = 255 - (255 - tNormed);
                        int b = 255 - (255 - tNormed);
                        fil << r << " " << g << " " << b << " ";
                    }

                fil.close();
                return 0;
            } catch (std::exception& e) {
                std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
                std::cerr << "Error: " << e.what() << std::endl;
                std::cerr << "\033[1;37mImage not writtenn.\033[0m" << std::endl;
                return -1;
            }
        }

        template<typename T>
        inline int temperatureStatistics(const PWM::PWMDataStructure::flatStaggeredGrid<T>& data)
        {
            float meant = 0.0f; float maxt = std::numeric_limits<float>::min();

            for (int i = 0; i < data.size(); ++i)
            {
                float t = (float) data.getData(i);
                meant += t;
                if(t > maxt)
                    maxt = t;
            }
            meant /= (float) data.size();
            std::cerr << "Temperature: mean = " << meant << " max = " << maxt << std::endl;
            return 0;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline int writeAirImages(std::string outputDirectory, int i, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            for (int j = 0; j < wld->getAirLayers().size(); ++j){
                int hght = (int) wld->getAirLayer(j)->getHeight();
                std::stringstream fT, fM, fCW, fV1, fV2, fCV, fCV1, fCV2, fP, fA;
                fT << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Temperature_Step_" << i << ".ppm";
                fM << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Moisture_Step_" << i << ".ppm";
                fCW << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_CondensedWater_Step_" << i << ".ppm";
                fV1 << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_VelocityY_Step_" << i << ".ppm";
                fV2 << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_VelocityX_Step_" << i << ".ppm";
                //                fP << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Pressure_Step_" << i << ".ppm";
                fA << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Particulates_Step_" << i << ".ppm";
                fCV << outputDirectory << "/ConvectionLayer_" << j << "_Vertical_Velocity_Step_" << i << ".ppm";
                fCV1 << outputDirectory << "/ConvectionLayer_" << j << "_Rainfall_Step_" << i << ".ppm";
                fCV2 << outputDirectory << "/ConvectionLayer_" << j << "_Ashfall_Step_" << i << ".ppm";
                writeTempImage<V>(fT.str(), wld->getAirLayer(j)->getTemperature());
                if (j < wld->getAirLayers().size() - 1){
//                    writeVelImage(fCV.str(), wld->getConvectionLayer(j)->getVerticalVelocities());
//                    writeMoisImage(fCV1.str(), wld->getConvectionLayer(j)->getRainfall());
//                    writeAshImage(fCV2.str(), wld->getConvectionLayer(j)->getAshfall());
                }

                writeMoisImage<V>(fCW.str(), wld->getAirLayer(j)->getCondensedWater());
                writeMoisImage<V>(fM.str(), wld->getAirLayer(j)->getMoisture());
//                writeVelImage<V>(fV1.str(), wld->getAirLayer(j)->getVelocityTheta());
//                writeVelImage<V>(fV2.str(), wld->getAirLayer(j)->getVelocityPhi());
//                writeAshImage(fA.str(), wld->getAirLayer(j)->getParticulates());
                writeCloudImages<V>(outputDirectory, i, j, wld->getAirLayer(j)->getHeight(), wld->getAirLayer(j)->getClouds());
            }
            return 0;
        }

        template<typename T, typename TT, typename V, typename VV>
        inline int writeAirImagesSelect(std::string outputDirectory, int i, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            for (int j = 0; j < wld->getAirLayers().size(); j+=2){ // ++j
                // JG change
                int hght = (int) wld->getAirLayer(j)->getHeight();
                std::stringstream fT, fM, fCW, fV1, fV2, fCV, fCV1, fCV2, fP, fA;
                fT << outputDirectory << "/Layer_" << j << "_(" << hght  << "m)_Temperature_Step_" << i << ".ppm";
                fM << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_Moisture_Step_" << i << ".ppm";
                fCW << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_CondensedWater_Step_" << i << ".ppm";
                fV1 << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_VelocityX_Step_" << i << ".ppm";
                fV2 << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_VelocityY_Step_" << i << ".ppm";
                //                fP << outputDirectory << "/Layer_" << j << "_(" << wld->getAirLayer(j)->getHeight() << "m)_Pressure_Step_" << i << ".ppm";
                fA << outputDirectory << "/Layer_" << j << "_(" << hght << "m)_Particulates_Step_" << i << ".ppm";
                fCV << outputDirectory << "/ConvectionLayer_" << j << "_Vertical_Velocity_Step_" << i << ".ppm";
                fCV1 << outputDirectory << "/ConvectionLayer_" << j << "_Rainfall_Step_" << i << ".ppm";
                fCV2 << outputDirectory << "/ConvectionLayer_" << j << "_Ashfall_Step_" << i << ".ppm";
                writeTempImage<V>(fT.str(), wld->getAirLayer(j)->getTemperature());
                if (j < wld->getAirLayers().size() - 1){
                    writeVelImage(fCV.str(), wld->getConvectionLayer(j)->getVerticalVelocities());
                    // writeMoisImage(fCV1.str(), wld->getConvectionLayer(j)->getRainfall());
                    // writeAshImage(fCV2.str(), wld->getConvectionLayer(j)->getAshfall());
                }
                writeMoisImage<V>(fCW.str(), wld->getAirLayer(j)->getCondensedWater());
                // writeMoisImage<V>(fM.str(), wld->getAirLayer(j)->getMoisture());
                // writeVelImage<V>(fV1.str(), wld->getAirLayer(j)->getVelocityTheta());
                // writeVelImage<V>(fV2.str(), wld->getAirLayer(j)->getVelocityPhi());
                writeAshImage(fA.str(), wld->getAirLayer(j)->getParticulates());
                // writeCloudImages<V>(outputDirectory, i, j, wld->getAirLayer(j)->getHeight(), wld->getAirLayer(j)->getClouds());
                // JG change to handle memory
            }
            return 0;
        }
    
        template<typename T, typename TT, typename V, typename VV>
        inline void airImagesAnalysis(int i, int cx, int cy, float r, const std::shared_ptr<PWM::Model::world<T, TT, V, VV>>& wld){
            
            // report proportion of clouds within a specified radius of the vent across all layers, and ashfall over the entire domain
            
            int cpixelcount = 0, cloudcount = 0, ashfallcount = 0, apixelcount = 0;
            for (int j = 0; j < wld->getAirLayers().size(); j++){
                auto & cmask = wld->getAirLayer(j)->getClouds();
                size_t nx = cmask.getX(), ny = cmask.getY();
               
                // look for clouds within a pixel radius of the vent
                for (int x = 0; x < nx; ++x){
                    for (int y = 0; y < ny; ++y){
                        float d = (cx - x) * (cx - x) + (cy - y) * (cy - y);
                        if(d < r*r)
                        {
                            auto cloudType = static_cast<int>(cmask.getData(x, y));
                            if(cloudType > 0)
                                cloudcount++;
                            cpixelcount++;
                        }
                    }
                }
                
                if(j < wld->getAirLayers().size()-1)
                {
                    // PWM::PWMDataStructure::square2DArray<T>&
                    auto & amask = wld->getConvectionLayer(j)->getAshfall();
                    nx = amask.getX(), ny = amask.getY();
                    
                    // search for ashfall over the entire domain
                    for (int x = 0; x < nx; ++x){
                        for (int y = 0; y < ny; ++y){
                            auto ashval = static_cast<int>(amask.getData(x, y));
                            if(ashval > 0.0000001f)
                                ashfallcount++;
                            apixelcount++;
                        }
                    }
                }
             
            }
            std::cerr << "TIMESTEP " << i << std::endl;
            std::cerr << "Cloud proportion = " << (float) cloudcount / (float) cpixelcount << std::endl;
            std::cerr << "Ashfall proportion = " << (float) ashfallcount / (float) apixelcount << std::endl;
            if(ashfallcount > 0 || cloudcount > 0)
                std::cerr << "*********************" << std::endl;
        }
    }
}
#endif // VIZUTILS_H
