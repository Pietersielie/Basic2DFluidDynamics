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

#include "Coordinate.h"
#include <cmath>

namespace BasicFluidDynamics{
    namespace Data{
		Coordinate::Coordinate(double lat, double longi): latitude(lat), longitude(longi) {
		}
		
		/*Coordinate::Coordinate(const Coordinate& other): latitude(other.getLatitude()), longitude(other.getLongitude()) {
		}
		
		Coordinate& Coordinate::operator= (const Coordinate& other){
			setLatitude(other.getLatitude());
			setLongitude(other.getLongitude());
			return *this;
		}*/
		
		bool Coordinate::operator==(const Coordinate& other) const{
			return (getLatitude() == other.getLatitude()) && (getLongitude() == other.getLongitude());
		}

		bool Coordinate::operator!=(const Coordinate& other) const{
			return !(*this == other);
		}

		bool Coordinate::approxEquals(const Coordinate& other) const{
			return (std::abs(getLatitude() - other.getLatitude()) < approxLimit) && (std::abs(getLongitude() - other.getLongitude()) < approxLimit);
		}
		
		double Coordinate::getLatitude() const{
			return this->latitude;
		}
		
		double Coordinate::getLongitude() const{
			return this->longitude;
		}
		
		void Coordinate::setLatitude(const double lat){
			this->latitude = lat;
		}
		void Coordinate::setLongitude(const double longi){
			this->longitude = longi;
		}
				
	}
}
