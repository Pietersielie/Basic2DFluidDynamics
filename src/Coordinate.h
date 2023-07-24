// Wrapper class to deal with all coordinate data (latitude/longitude)
// Cilliers Pretorius
// 07 May 2021

#ifndef PWM_PWMDATASTRUCTURE_COORDINATE_H
#define PWM_PWMDATASTRUCTURE_COORDINATE_H
namespace PWM {
	namespace PWMDataStructure{
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