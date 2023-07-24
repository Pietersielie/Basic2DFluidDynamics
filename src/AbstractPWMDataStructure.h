// This class serves as an abstract class for any data structure and should allow the other engines to be general.
// Cilliers Pretorius
// 05 May 2021
#ifndef PWM_ABSTRACTPWMDATASTRUCTURE_H
#define PWM_ABSTRACTPWMDATASTRUCTURE_H

#include "Coordinate.h"
#include <cstddef>
#include <iostream>
#include <limits>
#include "mathUtils.h"
#include <random>
#include <tuple>
#include <vector>

namespace PWM {
	namespace PWMDataStructure{
		template <typename T> class AbstractPWMDataStructure;

		template<typename T>
		void swap(AbstractPWMDataStructure<T>& i, AbstractPWMDataStructure<T>& j);

		template <typename T> class AbstractPWMDataStructure{
			private:
                T * data;
            public:
                size_t datasize;
                AbstractPWMDataStructure();
                AbstractPWMDataStructure(const size_t len);
                virtual ~AbstractPWMDataStructure() = 0;
                /*~AbstractPWMDataStructure();
                AbstractPWMDataStructure(const AbstractPWMDataStructure & other);
                AbstractPWMDataStructure(AbstractPWMDataStructure && other);*/
                virtual const Coordinate getCoordinates(const size_t index) const = 0;
                virtual const double gridLength() const = 0;

                void setData(const size_t index, const T& val);
                const T getData(const size_t index) const;
                const T sum() const;
                const T mean() const;
                const T stdDev() const;

                void mul(const T& val);
                void add(const AbstractPWMDataStructure & other);

                void randomInit(const T& min, const T& max);

                void copy(const T& val);
                void copy(const AbstractPWMDataStructure& other);
                bool checkSize(const AbstractPWMDataStructure& other) const;

                const size_t size() const;

                //Functions to help with mapping the dataset to image.
                const T min() const;
                const T max() const;
                const std::pair<T, T> minmax() const;
                const T range() const;

                virtual const std::string print() const = 0;
                virtual const T getInterpolated(const double x, const double y) const = 0;

                AbstractPWMDataStructure& operator=(const AbstractPWMDataStructure & other);
                AbstractPWMDataStructure& operator=(const T& val);
                bool operator==(const AbstractPWMDataStructure & other) const;
                bool operator!=(const AbstractPWMDataStructure & other) const;
                AbstractPWMDataStructure& operator+=(const AbstractPWMDataStructure & other);
                AbstractPWMDataStructure& operator*=(const T& val);

                friend void swap <> (AbstractPWMDataStructure<T>& i, AbstractPWMDataStructure<T>& j);
		};

		template<typename T>
    inline AbstractPWMDataStructure<T>::AbstractPWMDataStructure() : datasize(0){
		}
		
		template<typename T>
    inline AbstractPWMDataStructure<T>::AbstractPWMDataStructure(size_t len) : datasize(len){
			this->data = new T[len];
		}
		
        template<typename T>
    inline AbstractPWMDataStructure<T>::~AbstractPWMDataStructure(){
            delete[] this->data;
        }
		
		template<typename T>
    inline AbstractPWMDataStructure<T>& AbstractPWMDataStructure<T>::operator=(const AbstractPWMDataStructure<T> & other){
			if (*this == other)
				return *this;
			this->data = other.data;
			this->datasize = other.datasize;
			return *this;
		}

		template<typename T>
    inline bool AbstractPWMDataStructure<T>::operator==(const AbstractPWMDataStructure<T> & other) const {
			if (datasize != other.datasize)
				return false;
			for (int i = 0; i < datasize; ++i)
				if (data[i] != other.data[i])
					return false;
			return true;
		}
		
		template<typename T>
    inline bool AbstractPWMDataStructure<T>::operator!=(const AbstractPWMDataStructure<T> & other) const {
			return !(this == other);
		}

		template<typename T>
    inline const size_t AbstractPWMDataStructure<T>::size() const{
			return datasize;
		}

		template<typename T>
    inline bool AbstractPWMDataStructure<T>::checkSize(const AbstractPWMDataStructure<T>& other) const{
			if (size() != other.size())
				return false;
			return true;
		}

		template<typename T>
    inline const T AbstractPWMDataStructure<T>::getData(const size_t index) const {
			if (index < size())
				return this->data[index];
			else{
				std::cerr << "\033[1;31mError: Out of bounds array read access in AbstractPWMDataStructure::getData(const size_t index) with index = " << index << "! Returning default value.\033[0m" << std::endl;
				return T();
			}
		}

		template<typename T>
    inline void AbstractPWMDataStructure<T>::setData(const size_t index, const T& val) {
			if (index < size()){
				this->data[index] = val;
				return;
			}else
				std::cerr << "\033[1;31mError: Out of bounds array read access in AbstractPWMDataStructure::setData(const size_t index, const T& val) with index = " << index << "! Doing nothing.\033[0m" << std::endl;
		}

		template<typename T>
    inline const T AbstractPWMDataStructure<T>::min() const{
			T min = std::numeric_limits<T>::max();
			for (int i = 0; i < size(); ++i)
				if (this->data[i] < min)
					min = this->data[i];
			return min;
		}

		template<typename T>
    inline const T AbstractPWMDataStructure<T>::max() const{
			T max = std::numeric_limits<T>::min();
			for (int i = 0; i < size(); ++i)
				if (this->data[i] > max)
					max = this->data[i];
			return max;
		}

		template<typename T>
    inline const std::pair<T, T> AbstractPWMDataStructure<T>::minmax() const{
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
    inline const T AbstractPWMDataStructure<T>::range() const{
            auto r = this->minmax();
            return r.second - r.first;
        }

		template<typename T>
    inline void AbstractPWMDataStructure<T>::mul(const T& val){
			*this *= val;
		}

		template<typename T>
    inline void AbstractPWMDataStructure<T>::add(const AbstractPWMDataStructure<T>& other){
			*this += other;
		}

		template<typename T>
    inline const T AbstractPWMDataStructure<T>::sum() const{
			auto res = T();
			#pragma omp parallel for reduction(+: res)
				for (int i = 0; i < size(); ++i)
						res += this->data[i];
			#pragma omp barrier
			return res;
		}

		template<typename T>
    inline const T AbstractPWMDataStructure<T>::mean() const{
			T r;
			r = sum() / size();
            return r;
        }

        template<typename T>
    inline const T AbstractPWMDataStructure<T>::stdDev() const{
            T ave = this->mean();
            T var = 0;
            #pragma omp parallel for reduction(+: var)
                for (int i = 0; i < size(); ++i)
                        var += std::pow(this->data[i] - ave, 2);
            #pragma omp barrier
            return std::sqrt(var / size());
        }

		template<typename T>
    inline AbstractPWMDataStructure<T>& AbstractPWMDataStructure<T>::operator*=(const T& val){
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i){
					auto x = getData(i) * val;
					setData(i, x);
				}
			#pragma omp barrier
            return *this;
        }

        template<typename T>
    inline AbstractPWMDataStructure<T>& AbstractPWMDataStructure<T>::operator+=(const AbstractPWMDataStructure<T>& other){
            #pragma omp parallel for
                for (int i = 0; i < size(); ++i){
                    auto x = getData(i) + other.getData(i);
                    setData(i, x);
                }
            #pragma omp barrier
            return *this;
        }

        template<typename T>
    inline AbstractPWMDataStructure<T>& AbstractPWMDataStructure<T>::operator=(const T& val){
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i){
					setData(i, val);
				}
			#pragma omp barrier
			return *this;
		}

		template<typename T>
    inline void AbstractPWMDataStructure<T>::copy(const T& val){
			*this = val;
		}

		template<typename T>
    inline void AbstractPWMDataStructure<T>::copy(const AbstractPWMDataStructure<T>& other){
			#pragma omp parallel for
				for (int i = 0; i < size(); ++i){
					setData(i, other.getData(i));
				}
			#pragma omp barrier
		}

		template<typename T>
    inline void AbstractPWMDataStructure<T>::randomInit(const T& min, const T& max){
			std::mt19937_64 mt(time(NULL));
			auto d = PWM::Utils::safeMinMax(min, max);
			std::uniform_real_distribution<T> dist(std::get<0>(d), std::get<1>(d));
			for (int i = 0; i < size(); ++i)
				setData(i, (dist(mt) - std::get<2>(d)));
		}

		template<typename T>
    inline void swap(AbstractPWMDataStructure<T>& i, AbstractPWMDataStructure<T>& j){
			std::swap(i.data, j.data);
			std::swap(i.datasize, j.datasize);
		}

//		template<typename R, typename T, typename Function>
//		void map(std::vector<R>& result, const AbstractPWMDataStructure<T>& data, Function mapper){
//			if (result.size() != data.size()){
//				std::cerr << "In map(std::vector<R> result, AbstractPWMDataStructure<T> data, Function mapper), result and data do not have matching sizes!" << std::endl;
//				return;
//			}
//			#pragma omp parallel for
//				for (int i = 0; i < data.size(); ++i){
//					result[i] = mapper(data, i);
//				}
//			#pragma omp barrier
//		}
	}
}

#endif //PWM_ABSTRACTPWMDATASTRUCTURE_H
