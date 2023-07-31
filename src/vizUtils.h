#ifndef VIZUTILS_H
#define VIZUTILS_H

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include "flatStaggeredGrid.h"
#include <opencv2/imgcodecs.hpp>
#include <tuple>
#include <vector>

namespace BasicFluidDynamics {
    namespace Utils {
        template<typename T>
        inline int writeMoisImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T>& data){
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

        template<typename T1, typename T2>
        inline int writeMoisImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T1>& data, const BasicFluidDynamics::Data::flatStaggeredGrid<T2>& obs, const size_t xBuffer){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX - xBuffer - xBuffer, CV_8UC4);
            //Image data
            auto mM = std::make_pair<T1, T1>(0.1, 1.0);
            T1 max = 2.0;
            T1 spread = mM.second - mM.first;
            for (int i = 0; i < nY; ++i){
                for (int j = xBuffer; j < nX - xBuffer; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = g = b = 0;
                    }
                    else{
                        float moisNormed, mois = data.getData(i, j);
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
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j - xBuffer);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeTempImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T>& data){
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
                    int tNormed;
                    T t = data.getData(i);
                    if (t < mM.first){
                        tNormed = ((t) / mM.first) * 255;
                        r = 0;
                        g = 0;
                        b = 255 - (255 - tNormed);
                    }
                    else if (t < mM.second){
                        tNormed = ((t - mM.first) / spread) * 255;
                        r = 255 - (255 - tNormed);
                        g = 0;
                        b = 255 - tNormed;
                    }
                    else if (t < max){
                        tNormed = ((t - mM.second) / (max - mM.second)) * 255;
                        r = 255;
                        g = 255 - (255 - tNormed);
                        b = 255 - (255 - tNormed);
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

        template<typename T1, typename T2>
        inline int writeTempImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T1>& data, const BasicFluidDynamics::Data::flatStaggeredGrid<T2>& obs, const size_t xBuffer){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX - xBuffer - xBuffer, CV_8UC4);
            //Image data
            auto mM = std::make_pair<T1, T1>(240, 350);
            T1 spread = mM.second - mM.first;
            T1 max = 1200;
            for (int i = 0; i < nY; ++i){
                for (int j = xBuffer; j < nX - xBuffer; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = 100;
                        g = 182;
                        b = 106;
                    }
                    else{
                        int tNormed;
                        T1 t = data.getData(i, j);
                        if (t < mM.first){
                            tNormed = ((t) / mM.first) * 255;
                            r = 0;
                            g = 0;
                            b = 255 - (255 - tNormed);
                        }
                        else if (t < mM.second){
                            tNormed = ((t - mM.first) / spread) * 255;
                            r = 255 - (255 - tNormed);
                            g = 0;
                            b = 255 - tNormed;
                        }
                        else if (t < max){
                            tNormed = ((t - mM.second) / (max - mM.second)) * 255;
                            r = 255;
                            g = 255 - (255 - tNormed);
                            b = 255 - (255 - tNormed);
                        }
                        else{
                            r = g = b = 255;
                        }
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j - xBuffer);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writeVelImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T>& data){
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
                    if (BasicFluidDynamics::Utils::getSign(data.getData(i)) > 0){
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

        template<typename T, typename T2>
        inline int writeVelImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T>& data, const BasicFluidDynamics::Data::flatStaggeredGrid<T2>& obs, const size_t xBuffer){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX - xBuffer - xBuffer, CV_8UC4);
            //Image data
            T max = 140;
            int r, g, b;
            for (int i = 0; i < nY; ++i){
                for (int j = xBuffer; j < nX - xBuffer; ++j){
                    if (obs.getData(i, j)){
                        r = 166;
                        g = b = 39;
                    }
                    else{
                        T t = data.getData(i, j);
                        int tNormed = ((std::abs(t) / max) * 255);
                        r = 0;
                        if (BasicFluidDynamics::Utils::getSign(t) > 0){
                            g = 255 - std::min(255 - tNormed, 255);
                            b = 0;
                        }
                        else{
                            g = 0;
                            b = 255 - std::min(255 - tNormed, 255);
                        }
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j - xBuffer);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }

        template<typename T>
        inline int writePresImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T>& data){
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

        template<typename T, typename T2>
        inline int writePresImage(std::string file, const BasicFluidDynamics::Data::flatStaggeredGrid<T>& data, const BasicFluidDynamics::Data::flatStaggeredGrid<T2>& obs, const size_t xBuffer){
            std::ofstream fil;
            size_t nX = data.getX(), nY = data.getY();
            cv::Mat img(nY, nX - xBuffer - xBuffer, CV_8UC4);
            //Image data
            auto mM = data.minmax();
            T spread = mM.second - mM.first;
            for (int i = 0; i < nY; ++i){
                for (int j = xBuffer; j < nX - xBuffer; ++j){
                    int r, g, b;
                    if (obs.getData(i, j)){
                        r = g = 206;
                        b = 25;
                    }
                    else{
                        float tNormed = ((data.getData(i, j) - mM.first) / spread) * 255;
                        r = 255 - (255 - tNormed);
                        g = 255 - (255 - tNormed);
                        b = 255 - (255 - tNormed);
                    }
                    cv::Vec4b& bgra = img.at<cv::Vec4b>(i, j - xBuffer);
                    bgra[0] = b;
                    bgra[1] = g;
                    bgra[2] = r;
                    bgra[3] = UCHAR_MAX;
                }
            }
            std::vector<int> compression_params;
            compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
            compression_params.push_back(9);
            bool res = false;
            try {
                res = cv::imwrite(file, img, compression_params);
            } catch (std::exception& e) {
                std::cerr << "Exception converting image to PNG format: " << e.what() << std::endl;
            }
            if (res) return 0;

            std::cerr << "\033[1;31mError! Can't write image to file " << file << "!\033[0m" << std::endl;
            std::cerr << "\033[1;37mImage not written.\033[0m" << std::endl;
            return -1;
        }
    }
}
#endif // VIZUTILS_H
