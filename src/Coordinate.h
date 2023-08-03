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

#ifndef BFD_COORDINATE_H
#define BFD_COORDINATE_H
namespace BasicFluidDynamics {
    namespace Data{
		class Coordinate{
			private:
				double latitude, longitude;
				static constexpr double approxLimit = 0.001;
			public:
				void setLatitude(const double lat);
				void setLongitude(const double longi);
				double getLatitude() const;
				double getLongitude() const;
				Coordinate(double lat, double longi);
				/*Coordinate(const Coordinate& other);
				Coordinate & operator=(const Coordinate& other);*/
				bool operator==(const Coordinate& other) const;
				bool operator!=(const Coordinate& other) const;
				bool approxEquals(const Coordinate& other) const;
		};
	}
}

#endif //PWM_PWMDATASTRUCTURE_COORDINATE_H
