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

#ifndef BFD_UTILS_MATH
#define BFD_UTILS_MATH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#if __cplusplus > 201703L
	#include <numbers>
#endif
#include <string>
#include <tuple>
#include <vector>

namespace BasicFluidDynamics{
    namespace Utils{
        /**
         * Linear interpolation
         */
        template<typename T>
        inline const T interpolate(const T& fromPointA, const T& toPointB, float factor){
			return ((1.0 - factor) * fromPointA) + (factor * toPointB);
		}

        /**
         * Bilinear interpolation
         */
        template<typename T>
        inline const T interpolate(const T& fromPointA1, const T& fromPointA2, const T& toPointB1, const T& toPointB2, float factor1, float factor2){
			return interpolate(interpolate(fromPointA1, fromPointA2, factor1), interpolate(toPointB1, toPointB2, factor1), factor2);
		}

        /**
         * @brief smoothInterpolate SmootherStep 4-degree interpolation, see https://en.wikipedia.org/wiki/Perlin_noise#Implementation
         * @param fromPointA first value
         * @param toPointB second value
         * @param factor interpolation factor
         * @return the interpolated value
         */
        template<typename T>
        inline const T smoothInterpolate(const T& fromPointA, const T& toPointB, float factor){
            return (toPointB - fromPointA) * ((factor * (factor * 6.0 - 15.0) + 10.0) * factor * factor * factor) + fromPointA;
        }

        template<typename T, typename T1, typename T2>
        inline const T calcCartesianDistance(const std::pair<T1, T2>& x, const std::pair<T1, T2>& y){
            return std::sqrt(std::pow((y.first - x.first), 2) + std::pow((y.second - x.second), 2));
        }

        template<typename T1, typename T2, typename T3>
        inline const bool calcPointInEllipse(const std::pair<T1, T2>& origin, const std::pair<T1, T2>& point, const T3& height, const T3& width){
            auto a = height * 0.5;
            auto b = width * 0.5;
            auto a2 = std::pow(a, 2);
            auto b2 = std::pow(b, 2);
            return (a2 * std::pow(point.first - origin.first, 2) + b2 * std::pow(point.second - origin.second, 2)) <= a2 * b2;
        }

        template<typename T> inline constexpr
        const int getSign(const T& val, std::false_type is_signed){
            return T(0) < val;
        }

        template<typename T> inline constexpr
        const int getSign(const T& val, std::true_type is_signed){
            return (T(0) < val) - (val < T(0));
        }

        template<typename T> inline constexpr
        const int getSign(const T& val){
            return getSign(val, std::is_signed<T>());
        }

        template<typename T>
        inline const T radToDeg(const T& angleInRad){
            #if __cplusplus > 201703L
				return angleInRad * (180 / std::numbers::pi);
			#else
				return angleInRad * (180 / 3.1415926535);
			#endif
        }

        template<typename T>
        inline const T degToRad(const T& angleInDeg){
            #if __cplusplus > 201703L
				return angleInDeg * (std::numbers::pi / 180);
			#else
				return angleInDeg * (3.1415926535 / 180);
			#endif
        }

		template<typename T>
        inline std::tuple<T, T, T> safeMinMax(const T& min, const T& max){
			if (min < 0){
				return std::make_tuple((T) 0, max + std::abs(min), std::abs(min));
			}
			return std::make_tuple(min, max, (T) 0);
		}

        template<typename T>
        inline T elevThres(const float threshold, const T testElev, const T minH, const T topH, const std::vector<std::pair<T, bool>>& data){
            if ((topH - minH) < 100)
                return testElev;
            int count = 0;
            float thres = threshold * data.size();
            bool t = false;
            for (auto x : data){
                if (x.first >= testElev && !x.second){
                    ++count;
                    if (count >= thres){
                        t = true;
                        break;
                    }
                }
            }
            if (t){
                if (topH - testElev < 100)
                    return testElev;
                else
                    return elevThres(threshold, (topH + testElev) / 2, testElev, topH, data);
            }
            else
                return elevThres(threshold, (testElev + minH) / 2, minH, testElev, data);
        }

        template<typename T, typename U>
        typename std::enable_if<std::is_floating_point<T>::value, T>::type unsignedFloatingModulo(const T x, const U y){
            auto r1 = std::fmod(x, y);
            return (r1 < 0) ? r1 + y : r1;
        }

        template<typename T>
        inline T eigenVecDistance(Eigen::Vector3f& a, Eigen::Vector3f& b){
            T res = sqrt(std::pow(a.x() - b.x(), 2) + std::pow(a.y() - b.y(), 2) + std::pow(a.z() - b.z(), 2));
            return res;
        }

        inline const size_t convert2Dto1DUtil(const size_t rowCount, const size_t colCount, const int i, const int j){
            int i1, j1;
            i1 = (i < 0) ? i + rowCount : i;
            j1 = (j < 0) ? j + colCount : j;
            return (i1 % rowCount) * colCount + (j1 % colCount);
        }

        /**
         * @brief randomGradient Create a pseudorandom direction vector
         *        See https://en.wikipedia.org/wiki/Perlin_noise#Implementation
         */
        template<typename T, typename T1>
        inline std::pair<T, T> randomGradient(T1& ix, T1& iy) {
            // No precomputed gradients mean this works for any number of grid coordinates
            const unsigned w = 8 * sizeof(unsigned);
            const unsigned s = w / 2; // rotation width
            unsigned a = ix, b = iy;
            a *= 3284157443; b ^= a << s | a >> (w-s);
            b *= 1911520717; a ^= b << s | b >> (w-s);
            a *= 2048419325;
            T random = a * (3.14159265 / ~(~0u >> 1)); // in [0, 2*Pi]
            std::pair<T, T> v;
            v.first = std::cos(random); v.second = std::sin(random);
            return v;
        }


        /**
         * @brief dotGridGradient Computes the dot product of the distance and gradient vectors.
         *        See https://en.wikipedia.org/wiki/Perlin_noise#Implementation
         */
        template<typename T, typename T1>
        inline T dotGridGradient(T1 ix, T1 iy, T x, T y) {
            // Get gradient from integer coordinates
            std::pair<T, T> gradient = randomGradient<T, T1>(ix, iy);

            // Compute the distance vector
            T dx = x - (float)ix;
            T dy = y - (float)iy;

            // Compute the dot-product
            return (dx*gradient.first + dy*gradient.second);
        }

        /**
         * @brief perlin A basic perlin noise generator for a given point, with values between 0 and max
         *        See https://en.wikipedia.org/wiki/Perlin_noise#Implementation
         * @param x X coordinate for given location
         * @param y Y coordinate for given location
         * @param max The maximum value from this noise generator
         * @return A value between 0 and max.
         */
        template<typename T, typename T1>
        inline const T perlin(T1& x, T1& y, T& max){
            // Determine grid cell coordinates
            T1 x0 = (int) floor(x);
            T1 x1 = x0 + 1;
            T1 y0 = (int) floor(y);
            T1 y1 = y0 + 1;

            // Determine interpolation weights
            // Could also use higher order polynomial/s-curve here
            T sx = x - (T) x0;
            T sy = y - (T) y0;

            // Interpolate between grid point gradients
            T n0, n1, ix0, ix1, value;

            n0 = dotGridGradient(x0, y0, x, y);
            n1 = dotGridGradient(x1, y0, x, y);
            ix0 = smoothInterpolate(n0, n1, sx);

            n0 = dotGridGradient(x0, y1, x, y);
            n1 = dotGridGradient(x1, y1, x, y);
            ix1 = smoothInterpolate(n0, n1, sx);

            value = smoothInterpolate(ix0, ix1, sy);
            return (value * 0.5 + 0.5) * max; // Will return in range -1 to 1. To make it in range 0 to 1, multiply by 0.5 and add 0.5
        }

        template<typename T>
        const std::vector<std::pair<T, T>>& allocateLayers(const std::vector<T>& heights, std::vector<std::pair<T, T>>& result, const T& maxLayerThickness){
            result.clear();
            for (int i = 0; i < heights.size(); ++i){
                double size;
                if (i == 0)
                    size = std::min(heights[i] - 0, heights[i + 1] - heights[i]) * 0.5;
                else if (i == heights.size() - 1)
                    size = (heights[i] - heights[i - 1]) * 0.5;
                else
                    size = std::min(heights[i] - heights[i - 1], heights[i + 1] - heights[i]) * 0.5;
                result.push_back(std::make_pair(heights[i]-size, heights[i]+size));
            }

            for (int i = 1; i < heights.size(); i += 2){
                T x = result[i].first;
                T y = result[i].second;
                while ((x != result[i - 1].second) && (y != result[i + 1].first)){
                    x -= 1;
                    y += 1;
                    if ((x < result[i - 1].second) || (y > result[i + 1].first)){
                        for (int i = 0; i < result.size(); ++i)
                            std::cerr << "[" << result[i].first << ", " << result[i].second << std::endl;
                        assert(0);
                    }
                }
                result[i].first = x, result[i].second = y;
            }

            auto padder = std::vector<std::pair<T, int>>();
            for (int i = 0; i < result.size(); ++i)
                padder.push_back(std::make_pair(result[i].second - result[i].first, i));
            std::sort(padder.begin(), padder.end(), [](std::pair<T, int> a, std::pair<T, int> b){ return a.first > b.first; });

            for (auto i : padder){
                T j = i.second;
                T upper_diff, lower_diff;
                if (j == 0){
                    upper_diff = result[j + 1].first - result[j].second;
                    if (upper_diff != 0)
                        result[j].second += upper_diff;
                    result[j].first = 0;
                    continue;

                }
                if (j == heights.size() - 1){
                    T lower_diff = result[j].first - result[j - 1].second;
                    if (lower_diff != 0)
                        result[j].first -= lower_diff;
                    continue;
                }
                lower_diff = result[j].first - result[j - 1].second;
                upper_diff = result[j + 1].first - result[j].second;
                if (lower_diff != 0)
                    result[j].first -= lower_diff;
                if (upper_diff != 0)
                    result[j].second += upper_diff;
            }
            auto thicknesses = std::vector<T>();
            for (int i = 0; i < result.size(); ++i)
                thicknesses.push_back(result[i].second - result[i].first);
            auto big = std::max_element(thicknesses.begin(), thicknesses.end());
            bool split = *big > maxLayerThickness;
            while (split){
                std::vector<int> rem = std::vector<int>();
                for (int i = 0; i < result.size(); ++i){
                    T thick = result[i].second - result[i].first;
                    if (thick > maxLayerThickness){
                        T newMid = std::round((result[i].first + result[i].second) * 0.5);
                        auto newBot = std::make_pair(result[i].first, newMid);
                        auto newTop = std::make_pair(newMid, result[i].second);
                        result.insert(result.begin() + i + 1, newBot);
                        result.insert(result.begin() + i + 2, newTop);
                        rem.push_back(i);
                    }
                }
                int rems = 0;
                for (auto i : rem){
                    result.erase(result.begin() + i - rems);
                    ++rems;
                }
                thicknesses.clear();
                for (int i = 0; i < result.size(); ++i)
                    thicknesses.push_back(result[i].second - result[i].first);
                auto newBig = std::max_element(thicknesses.begin(), thicknesses.end());
                split = *newBig > maxLayerThickness;
            }

            return result;
        }
    }
}
#endif //BFD_UTILS_MATH
