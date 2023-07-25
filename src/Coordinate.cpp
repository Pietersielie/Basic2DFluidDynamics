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
