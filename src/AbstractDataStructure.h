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
#ifndef BFD_ABSTRACTPWMDATASTRUCTURE_H
#define BFD_ABSTRACTPWMDATASTRUCTURE_H

#include "Coordinate.h"
#include <cstddef>
#include <iostream>
#include <limits>
#include "mathUtils.h"
#include <random>
#include <tuple>
#include <vector>

namespace BasicFluidDynamics {
    namespace Data{
        template <typename T> class AbstractDataStructure;

		template<typename T>
        void swap(AbstractDataStructure<T>& i, AbstractDataStructure<T>& j);

        template <typename T> class AbstractDataStructure{
			private:
                T * data;
            public:
                size_t datasize;
                AbstractDataStructure();
                AbstractDataStructure(const size_t len);
                virtual ~AbstractDataStructure() = 0;
                virtual const Coordinate getCoordinates(const size_t index) const = 0;
                virtual const double gridLength() const = 0;

                void setData(const size_t index, const T& val);
                const T getData(const size_t index) const;
                const T sum() const;
                const T mean() const;
                const T stdDev() const;

                void mul(const T& val);
                void add(const AbstractDataStructure & other);

                void randomInit(const T& min, const T& max);

                void copy(const T& val);
                void copy(const AbstractDataStructure& other);
                bool checkSize(const AbstractDataStructure& other) const;

                const size_t size() const;

                //Functions to help with mapping the dataset to image.
                const T min() const;
                const T max() const;
                const std::pair<T, T> minmax() const;
                const T range() const;

                virtual const std::string print() const = 0;
                virtual const T getInterpolated(const double x, const double y) const = 0;

                AbstractDataStructure& operator=(const AbstractDataStructure & other);
                AbstractDataStructure& operator=(const T& val);
                bool operator==(const AbstractDataStructure & other) const;
                bool operator!=(const AbstractDataStructure & other) const;
                AbstractDataStructure& operator+=(const AbstractDataStructure & other);
                AbstractDataStructure& operator*=(const T& val);

                friend void swap <> (AbstractDataStructure<T>& i, AbstractDataStructure<T>& j);
		};

		template<typename T>
        inline AbstractDataStructure<T>::AbstractDataStructure() : datasize(0){
		}
		
		template<typename T>
        inline AbstractDataStructure<T>::AbstractDataStructure(size_t len) : datasize(len){
			this->data = new T[len];
		}
		
        template<typename T>
        inline AbstractDataStructure<T>::~AbstractDataStructure(){
            delete[] this->data;
        }
		
		template<typename T>
        inline AbstractDataStructure<T>& AbstractDataStructure<T>::operator=(const AbstractDataStructure<T> & other){
			if (*this == other)
				return *this;
			this->data = other.data;
			this->datasize = other.datasize;
			return *this;
		}

		template<typename T>
        inline bool AbstractDataStructure<T>::operator==(const AbstractDataStructure<T> & other) const {
			if (datasize != other.datasize)
				return false;
			for (int i = 0; i < datasize; ++i)
				if (data[i] != other.data[i])
					return false;
			return true;
		}
		
		template<typename T>
        inline bool AbstractDataStructure<T>::operator!=(const AbstractDataStructure<T> & other) const {
			return !(this == other);
		}

		template<typename T>
        inline const size_t AbstractDataStructure<T>::size() const{
			return datasize;
		}

		template<typename T>
        inline bool AbstractDataStructure<T>::checkSize(const AbstractDataStructure<T>& other) const{
			if (size() != other.size())
				return false;
			return true;
		}

		template<typename T>
        inline const T AbstractDataStructure<T>::getData(const size_t index) const {
			if (index < size())
				return this->data[index];
			else{
                std::cerr << "\033[1;31mError: Out of bounds array read access in AbstractDataStructure::getData(const size_t index) with index = " << index << "! Returning default value.\033[0m" << std::endl;
				return T();
			}
		}

		template<typename T>
        inline void AbstractDataStructure<T>::setData(const size_t index, const T& val) {
			if (index < size()){
				this->data[index] = val;
				return;
			}else
                std::cerr << "\033[1;31mError: Out of bounds array read access in AbstractDataStructure::setData(const size_t index, const T& val) with index = " << index << "! Doing nothing.\033[0m" << std::endl;
		}

		template<typename T>
        inline const T AbstractDataStructure<T>::min() const{
			T min = std::numeric_limits<T>::max();
			for (int i = 0; i < size(); ++i)
				if (this->data[i] < min)
					min = this->data[i];
			return min;
		}

		template<typename T>
        inline const T AbstractDataStructure<T>::max() const{
			T max = std::numeric_limits<T>::min();
			for (int i = 0; i < size(); ++i)
				if (this->data[i] > max)
					max = this->data[i];
			return max;
		}

		template<typename T>
        inline const std::pair<T, T> AbstractDataStructure<T>::minmax() const{
			T min = std::numeric_limits<T>::max();
			T max = std::numeric_limits<T>::min();
			for (int i = 0; i < size(); ++i){
				if (this->data[i] < min)
					min = this->data[i];
				if (this->data[i] > max)
					max = this->data[i];
			}
            return std::make_pair(min, max);
        }

        template<typename T>
        inline const T AbstractDataStructure<T>::range() const{
            auto r = this->minmax();
            return r.second - r.first;
        }

		template<typename T>
        inline void AbstractDataStructure<T>::mul(const T& val){
			*this *= val;
		}

		template<typename T>
        inline void AbstractDataStructure<T>::add(const AbstractDataStructure<T>& other){
			*this += other;
		}

		template<typename T>
        inline const T AbstractDataStructure<T>::sum() const{
			auto res = T();
			#pragma omp parallel for reduction(+: res)
				for (int i = 0; i < size(); ++i)
						res += this->data[i];
			#pragma omp barrier
			return res;
		}

		template<typename T>
        inline const T AbstractDataStructure<T>::mean() const{
			T r;
			r = sum() / size();
            return r;
        }

        template<typename T>
        inline const T AbstractDataStructure<T>::stdDev() const{
            T ave = this->mean();
            T var = 0;
            #pragma omp parallel for reduction(+: var)
                for (int i = 0; i < size(); ++i)
                        var += std::pow(this->data[i] - ave, 2);
            #pragma omp barrier
            return std::sqrt(var / size());
        }

		template<typename T>
        inline AbstractDataStructure<T>& AbstractDataStructure<T>::operator*=(const T& val){
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i){
					auto x = getData(i) * val;
					setData(i, x);
				}
			#pragma omp barrier
            return *this;
        }

        template<typename T>
        inline AbstractDataStructure<T>& AbstractDataStructure<T>::operator+=(const AbstractDataStructure<T>& other){
            #pragma omp parallel for
                for (int i = 0; i < size(); ++i){
                    auto x = getData(i) + other.getData(i);
                    setData(i, x);
                }
            #pragma omp barrier
            return *this;
        }

        template<typename T>
        inline AbstractDataStructure<T>& AbstractDataStructure<T>::operator=(const T& val){
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i){
					setData(i, val);
				}
			#pragma omp barrier
			return *this;
		}

		template<typename T>
        inline void AbstractDataStructure<T>::copy(const T& val){
			*this = val;
		}

		template<typename T>
        inline void AbstractDataStructure<T>::copy(const AbstractDataStructure<T>& other){
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i){
					setData(i, other.getData(i));
				}
			#pragma omp barrier
		}

		template<typename T>
        inline void AbstractDataStructure<T>::randomInit(const T& min, const T& max){
			std::mt19937_64 mt(time(NULL));
            auto d = BasicFluidDynamics::Utils::safeMinMax(min, max);
			std::uniform_real_distribution<T> dist(std::get<0>(d), std::get<1>(d));
			for (int i = 0; i < size(); ++i)
				setData(i, (dist(mt) - std::get<2>(d)));
		}

		template<typename T>
        inline void swap(AbstractDataStructure<T>& i, AbstractDataStructure<T>& j){
			std::swap(i.data, j.data);
			std::swap(i.datasize, j.datasize);
		}
	}
}

#endif //PWM_ABSTRACTPWMDATASTRUCTURE_H
