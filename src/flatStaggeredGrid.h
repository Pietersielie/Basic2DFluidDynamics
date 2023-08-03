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

#ifndef BFD_FLATSTAGGEREDGRID_H
#define BFD_FLATSTAGGEREDGRID_H

#include "AbstractDataStructure.h"
#include <memory>

namespace BasicFluidDynamics {
    namespace Data {
        template<typename T>
        class flatStaggeredGrid;

        template<typename T>
        void swap(flatStaggeredGrid<T>& i, flatStaggeredGrid<T>& j);

        template<typename T>
        class flatStaggeredGrid : public AbstractDataStructure<T> {
            private:
                size_t nX; //Dimensions in the X direction (width)
                size_t nY; //Dimensions in the Y direction (height)

                double x_offset; //Offset in the X direction, in cell width units
                double y_offset; //Offset in the Y direction, in cell width units

                double worldXLength; //Length of the world in the X direction, in metres.
                double worldYLength; //Length of the world in the Y direction, in metres.

                Coordinate coord = Coordinate(0, 0);
            public:
                inline flatStaggeredGrid();
                inline flatStaggeredGrid(const flatStaggeredGrid& other);
                inline flatStaggeredGrid(const size_t width, const size_t height, double xOff = 0.5f, double yOff = 0.5f, double wrldX = -1, double wrldY = -1);

                inline const size_t getX() const;
                inline const size_t getY() const;
                inline const double getXOffset() const;
                inline const double getYOffset() const;
                inline const std::pair<double, double> getOffset() const;
                inline const double getXWorldLength() const;
                inline const double getYWorldLength() const;

                inline const size_t convert2Dto1D(const size_t i, const size_t j) const;
                inline const std::pair<size_t, size_t> convert1Dto2D(const size_t index) const;

                inline const double gridLength() const override;
                inline const double cellSizeX() const;
                inline const double cellSizeY() const;
                inline const Coordinate getCoordinates() const;
                inline const Coordinate getCoordinates(const size_t index) const override;
                inline const Coordinate getCoordinates(const size_t i, const size_t j) const;

                inline const std::pair<double, double> getWorldLoc(const size_t i, const size_t j) const;
                inline const std::pair<double, double> getWorldLoc(const size_t index) const;

                inline const T sampleAt(const std::pair<double, double> xy) const;
                inline const T sampleAt(const std::pair<double, double> xy, std::shared_ptr<flatStaggeredGrid<T>>& obs) const;
                inline const T sampleAt(const size_t i, const size_t j) const;
                inline const T sampleAt(const size_t i, const size_t j, std::shared_ptr<flatStaggeredGrid<T>>& obs) const;
                inline const T getInterpolated(const double x, const double y) const override;

                inline const std::string print() const override;

                inline const T getData(const size_t index) const;
                inline const T getData(const size_t i, const size_t j) const;
                inline void setData(const size_t index, const T& val);
                inline void setData(const size_t i, const size_t j, const T& val);
                inline flatStaggeredGrid& operator=(const T& val);
                flatStaggeredGrid operator-() const;
                inline void movingAverageSmoothing(std::shared_ptr<flatStaggeredGrid<T>>& dest, std::shared_ptr<flatStaggeredGrid<T>>& obs, const float window);
                inline void movingAverageSmoothing(std::shared_ptr<flatStaggeredGrid<T>>& dest, const float window) const;

                inline bool operator==(flatStaggeredGrid<T>& other) const;
                inline bool operator!=(flatStaggeredGrid<T>& other) const;
        };

        template<typename T>
        flatStaggeredGrid<T>::flatStaggeredGrid() : AbstractDataStructure<T>(){
            worldXLength = 0;
            worldYLength = 0;
            nX = 0;
            nY = 0;
            x_offset = 0;
            y_offset = 0;
        }

        template<typename T>
        flatStaggeredGrid<T>::flatStaggeredGrid(const flatStaggeredGrid &other) : AbstractDataStructure<T>(other.getX() * other.getY()), nX(other.getX()), nY(other.getY()), x_offset(other.getXOffset()), y_offset(other.getYOffset()), worldXLength(other.getXWorldLength()), worldYLength(other.getYWorldLength()){
            this->copy(other);
        }

        template<typename T>
        flatStaggeredGrid<T>::flatStaggeredGrid(const size_t width, const size_t height, double xOff, double yOff, double wrldX, double wrldY) : AbstractDataStructure<T>(width * height), nX(width), nY(height), x_offset(xOff), y_offset(yOff){
            if (wrldX <= 0)
                wrldX = width;
            if (wrldY <= 0)
                wrldY = height;
            worldXLength = wrldX;
            worldYLength = wrldY;
        }

        template<typename T>
        const size_t flatStaggeredGrid<T>::getX() const{
            return nX;
        }

        template<typename T>
        const size_t flatStaggeredGrid<T>::getY() const{
            return nY;
        }

        template<typename T>
        const double flatStaggeredGrid<T>::getXOffset() const{
           return x_offset;
        }

        template<typename T>
        const double flatStaggeredGrid<T>::getYOffset() const{
            return y_offset;
        }

        template<typename T>
        const std::pair<double, double> flatStaggeredGrid<T>::getOffset() const{
            return std::make_pair(x_offset, y_offset);
        }

        template<typename T>
        const double flatStaggeredGrid<T>::getXWorldLength() const{
           return worldXLength;
        }

        template<typename T>
        const double flatStaggeredGrid<T>::getYWorldLength() const{
            return worldYLength;
        }

        template<typename T>
        const size_t flatStaggeredGrid<T>::convert2Dto1D(const size_t i, const size_t j) const{
            int i1, j1;
            i1 = (i < 0) ? i + getY() : i;
            j1 = (j < 0) ? j + getX() : j;
            return (i1 % getY()) * getX() + (j1 % getX());
        }

        template<typename T>
        const std::pair<size_t, size_t> flatStaggeredGrid<T>::convert1Dto2D(const size_t index) const{
            size_t i, j;
            i = index / getX();
            j = index % getX();
            return std::make_pair(i, j);
        }

        template<typename T>
        const double flatStaggeredGrid<T>::gridLength() const{
            return worldXLength / nX;
        }

        template<typename T>
        const double flatStaggeredGrid<T>::cellSizeX() const{
            return worldXLength / nX;
        }

        template<typename T>
        const double flatStaggeredGrid<T>::cellSizeY() const{
            return worldYLength / nY;
        }

        template<typename T>
        const Coordinate flatStaggeredGrid<T>::getCoordinates() const{
            return coord;
        }

        template<typename T>
        const Coordinate flatStaggeredGrid<T>::getCoordinates(const size_t index) const{
            return getCoordinates();
        }

        template<typename T>
        const Coordinate flatStaggeredGrid<T>::getCoordinates(const size_t i, const size_t j) const{
            return getCoordinates(convert2Dto1D(i, j));
        }

        template<typename T>
        const std::pair<double, double> flatStaggeredGrid<T>::getWorldLoc(const size_t i, const size_t j) const{
            return std::make_pair(i + getXOffset(), j + getYOffset());
        }

        template<typename T>
        const std::pair<double, double> flatStaggeredGrid<T>::getWorldLoc(const size_t index) const{
            auto ij = convert1Dto2D(index);
            return getWorldLoc(ij.first, ij.second);
        }

        template<typename T>
        const T flatStaggeredGrid<T>::sampleAt(const std::pair<double, double> xy) const{
            double newX = xy.first - getXOffset();
            double newY = xy.second - getYOffset();
            int xIndex = (int) std::floor(newX);
            int yIndex = (int) std::floor(newY);
            double alphaX = newX - xIndex, alphaY = newY - yIndex;
            int x1, x2, y1, y2;
            x1 = xIndex % getY();
            x2 = (x1 + 1) % getY();
            y1 = yIndex % getX();
            y2 = (y1 + 1) % getX();

            T val11, val12, val21, val22;
            val11 = this->getData(x1, y1);
            val12 = this->getData(x1, y2);
            val21 = this->getData(x2, y1);
            val22 = this->getData(x2, y2);
            return BasicFluidDynamics::Utils::interpolate(val11, val12, val21, val22, alphaY, alphaX);
        }

        template<typename T>
        const T flatStaggeredGrid<T>::sampleAt(const std::pair<double, double> xy, std::shared_ptr<flatStaggeredGrid<T>>& obs) const{
            double newX = xy.first - getXOffset();
            double newY = xy.second - getYOffset();
            int xIndex = (int) std::floor(newX);
            int yIndex = (int) std::floor(newY);
            double alphaX = newX - xIndex, alphaY = newY - yIndex;
            int x1, x2, y1, y2;
            x1 = xIndex % getY();
            x2 = (x1 + 1) % getY();
            y1 = yIndex % getX();
            y2 = (y1 + 1) % getX();

            int obs11, obs12, obs21, obs22;
            obs11 = std::round(obs->sampleAt(x1, y1));
            obs12 = std::round(obs->sampleAt(x1, y2));
            obs21 = std::round(obs->sampleAt(x2, y1));
            obs22 = std::round(obs->sampleAt(x2, y2));

            T val11, val12, val21, val22;
            val11 = this->getData(x1, y1);
            val12 = this->getData(x1, y2);
            val21 = this->getData(x2, y1);
            val22 = this->getData(x2, y2);

            T val1, val2;
            if (obs11 && obs12){
                if (obs21 && obs22){
                    return 0;
                }
                else if (obs21){
                    return val22;
                }
                else if (obs22){
                    return val21;
                }
                else{
                    return BasicFluidDynamics::Utils::interpolate(val21, val22, alphaY);
                }
            }
            else if (obs11){
                if (obs21 && obs22){
                    return val12;
                }
                else if (obs21){
                    return BasicFluidDynamics::Utils::interpolate(val12, val22, alphaX);
                }
                else if (obs22){
                    return BasicFluidDynamics::Utils::interpolate(val12, val21, alphaX);
                }
                else {
                    val1 = val12;
                    val2 = BasicFluidDynamics::Utils::interpolate(val21, val22, alphaY);
                    return BasicFluidDynamics::Utils::interpolate(val1, val2, alphaX);
                }
            }
            else if (obs12){
                if (obs21 && obs22){
                    return val11;
                }
                else if (obs21){
                    return BasicFluidDynamics::Utils::interpolate(val11, val22, alphaX);
                }
                else if (obs22){
                    return BasicFluidDynamics::Utils::interpolate(val11, val21, alphaX);
                }
                else {
                    val1 = val11;
                    val2 = BasicFluidDynamics::Utils::interpolate(val21, val22, alphaY);
                    return BasicFluidDynamics::Utils::interpolate(val1, val2, alphaX);
                }
            }
            else{
                return BasicFluidDynamics::Utils::interpolate(val11, val12, val21, val22, alphaY, alphaX);
            }
        }

        template<typename T>
        const T flatStaggeredGrid<T>::sampleAt(const size_t i, const size_t j) const{
            return sampleAt(getWorldLoc(i, j));
        }

        template<typename T>
        const T flatStaggeredGrid<T>::sampleAt(const size_t i, const size_t j, std::shared_ptr<flatStaggeredGrid<T>>& obs) const{
            return sampleAt(getWorldLoc(i, j), obs);
        }

        template<typename T>
        const T flatStaggeredGrid<T>::getInterpolated(const double x, const double y) const{
            return sampleAt(std::make_pair(x, y));
        }

        template<>
        inline const std::string flatStaggeredGrid<std::string>::getInterpolated(const double x, const double y) const{
            int i = std::round(x), j = std::round(y);
            return getData(i, j);
        }

        template<typename T>
        const std::string flatStaggeredGrid<T>::print() const{
            std::ostringstream r;
            r.precision(5);
            for (int i = 0; i < getY(); ++i){
                r << getData(i, 0);
                for (int j = 1; j < getX(); ++j){
                    r << '\t' << getData(i, j);
                }
                r << '\n';
            }
            return r.str();
        }

        template<typename T>
        const T flatStaggeredGrid<T>::getData(const size_t index) const{
            return AbstractDataStructure<T>::getData(index);
        }

        template<typename T>
        const T flatStaggeredGrid<T>::getData(const size_t i, const size_t j) const{
            return AbstractDataStructure<T>::getData(convert2Dto1D(i, j));
        }

        template<typename T>
        void flatStaggeredGrid<T>::setData(const size_t index, const T& val){
            AbstractDataStructure<T>::setData(index, val);
        }

        template<typename T>
        void flatStaggeredGrid<T>::setData(const size_t i, const size_t j, const T& val){
            AbstractDataStructure<T>::setData(convert2Dto1D(i, j), val);
        }

        template<typename T>
        flatStaggeredGrid<T>& flatStaggeredGrid<T>::operator=(const T& val){
            this->copy(val);
            return *this;
        }

        template<typename T>
        flatStaggeredGrid<T> flatStaggeredGrid<T>::operator-() const{
            flatStaggeredGrid<T> res(*this);
            res.mul(-1);
            return res;
        }

        template<typename T>
        void flatStaggeredGrid<T>::movingAverageSmoothing(std::shared_ptr<flatStaggeredGrid<T>> &dest, std::shared_ptr<flatStaggeredGrid<T>> &obs, const float window){
            int winInt = std::ceil(window);
            if (winInt < 1){
                dest->copy(*this);
                return;
            }
            #pragma omp parallel for
                for (int i = 0; i < getY(); ++i) {
                    for (int j = 0; j < getX(); ++j){
                        if (obs->getData(i, j)){
                            dest->setData(i, j, this->getData(i, j));
                            continue;
                        }
                        T windowSum = 0;
                        int cellsinWin = 0;
                        auto loc = std::make_pair(i, j);
                        for (int y1 = i - winInt; y1 <= i + winInt; ++y1){
                            for (int x1 = j - winInt; x1 <= j + winInt; ++x1){
                                if (obs->getData(y1, x1)){
                                    continue;
                                }
                                auto pos = std::make_pair(y1, x1);
                                float dist = BasicFluidDynamics::Utils::calcCartesianDistance<float>(loc, pos);
                                if (dist <= window){
                                    windowSum += this->getData(y1, x1);
                                    ++cellsinWin;
                                }
                            }
                        }
                        T newVal = windowSum / cellsinWin;
                        dest->setData(i, j, newVal);
                    }
                }
            #pragma omp barrier
        }

        template<typename T>
        void flatStaggeredGrid<T>::movingAverageSmoothing(std::shared_ptr<flatStaggeredGrid<T>> &dest, const float window) const{
            int winInt = std::ceil(window);
            if (winInt < 1){
                dest->copy(*this);
                return;
            }
            #pragma omp parallel for
                for (int i = 0; i < getY(); ++i) {
                    for (int j = 0; j < getX(); ++j){
                        T windowSum = 0;
                        int cellsinWin = 0;
                        auto loc = std::make_pair(i, j);
                        for (int y1 = i - winInt; y1 <= i + winInt; ++y1){
                            for (int x1 = j - winInt; x1 <= j + winInt; ++x1){
                                auto pos = std::make_pair(y1, x1);
                                float dist = BasicFluidDynamics::Utils::calcCartesianDistance<float>(loc, pos);
                                if (dist <= window){
                                    windowSum += this->getData(y1, x1);
                                    ++cellsinWin;
                                }
                            }
                        }
                        T newVal = windowSum / cellsinWin;
                        dest->setData(i, j, newVal);
                    }
                }
            #pragma omp barrier
        }
    }
}

#endif // PWM_FLATSTAGGEREDGRID_H
