/////////////////////////////////////////////////
//cVector- by Matthew Bales
//
//Defines 2D and 3D physics vector classes and allows
//normal operator manimulation of them
//in addition to standard vector operations
///////////////////////////////////////////////////
#ifndef CVECTOR_H_INCLUDED
#define CVECTOR_H_INCLUDED
#pragma once

//Standard Libraries
#include <cmath>
#include <math.h>

//Matt Libraries
//#include "constants.h"

class cVector3D{
	friend cVector3D operator * (const double&, const cVector3D&);
public:
	double x, y, z;


	cVector3D(const double inpX, const double inpY, const double inpZ){x=inpX;	y=inpY;	z=inpZ;}
	cVector3D(){x=0; y=0; z=0;}
	void get(double& outX, double& outY, double& outZ);
	const cVector3D operator + (const cVector3D&) const;
	const cVector3D operator - (const cVector3D&) const;
	const cVector3D& operator = (const cVector3D&);
	const cVector3D operator * (const double) const;
	const cVector3D operator / (const double) const;
	const cVector3D& operator += (const cVector3D&);
	const cVector3D& operator -= (const cVector3D&);
	const cVector3D& operator *= (const double);
	const cVector3D& operator /= (const double);
	const bool operator == (const cVector3D&);
	const bool operator != (const cVector3D&);
	double mag();
	double dotProd(const cVector3D&);
	cVector3D crossProd(const cVector3D&);
	cVector3D scale(double);
	double gamma(); //For velocity vectors
	void setVal(const double, const double, const double);
	const cVector3D multComps(cVector3D inp) const;
	const cVector3D divComps(cVector3D inp) const;
	const cVector3D returnTransCoords() const;  //Defined by RDK2
	const cVector3D returnNormCoords() const;  //Defined by RDK2
	cVector3D norm();
	const double minComp();
	const double maxComp();
};

cVector3D zeroVec3D();

class cVector2D{
	friend cVector2D operator * (const double&, const cVector2D&);
public:
	double r,z;
	cVector2D(const double inpR, const double inpZ){r=inpR; z=inpZ;}
    cVector2D(){r=0; z=0;}
	const cVector2D operator + (const cVector2D&) const;
	const cVector2D operator - (const cVector2D&) const;
	const cVector2D& operator = (const cVector2D&);
	const cVector2D operator * (const double) const;
	const cVector2D operator / (const double) const;
	const cVector2D& operator += (const cVector2D&);
	const cVector2D& operator -= (const cVector2D&);
	const cVector2D& operator *= (const double);
	const cVector2D& operator /= (const double);
	const bool operator == (const cVector2D&);
	const bool operator != (const cVector2D&);
	double mag();
	double dotProd(const cVector2D&);
	cVector2D scale(const double);
	double gamma(); //For velocity vectors
	void setVal(const double, const double);
	const cVector2D multComps(cVector2D inp) const;
	const cVector2D divComps(cVector2D inp) const;
};

cVector2D zeroVec2D();


inline void cVector3D::get(double& outX, double& outY, double& outZ){
    outX=x;
    outY=y;
    outZ=z;
}

inline const cVector3D cVector3D::operator+(const cVector3D& rhs) const
{
	cVector3D result = *this;
	result += rhs;
	return result;
}

inline const cVector3D cVector3D::operator-(const cVector3D& rhs) const
{
	cVector3D result = *this;
	result -= rhs;
	return result;
}

inline const cVector3D& cVector3D::operator=(const cVector3D& inpVec)
{
	x = inpVec.x;
	y = inpVec.y;
	z = inpVec.z;
	return *this;
}

inline const cVector3D& cVector3D::operator += (const cVector3D& rhs)
{
	x +=rhs.x;
	y +=rhs.y;
	z +=rhs.z;
	return *this;
}
inline const cVector3D& cVector3D::operator -= (const cVector3D& rhs)
{
	x -=rhs.x;
	y -=rhs.y;
	z -=rhs.z;
	return *this;
}
inline const cVector3D& cVector3D::operator *= (double rhs)
{
	x *= rhs;
	y *= rhs;
	z *= rhs;
	return *this;
}

//Should throw exception if rhs is 0...but doesn't yet
inline const cVector3D& cVector3D::operator /= (double rhs)
{
	x /=rhs;
	y /=rhs;
	z /=rhs;
	return *this;
}


inline const cVector3D cVector3D::operator * (const double rhs) const
{
	cVector3D result = *this;
	result *= rhs;
	return result;
}

inline cVector3D operator* (const double& lhs, const cVector3D& rhs)
{
	cVector3D result = rhs;
	result *= lhs;
	return result;
}

inline const cVector3D cVector3D::operator / (const double rhs) const
{
	cVector3D result = *this;
	result /= rhs;
	return result;
}

inline const bool cVector3D::operator== (const cVector3D& inpVec)
{
	if(x==inpVec.x && y==inpVec.y && z==inpVec.z){
		return true;
	}
	else{
		return false;
	}
}

inline const bool cVector3D::operator != (const cVector3D& rhs)
{
	return !(*this == rhs);
}

inline double cVector3D::mag()
{
	return sqrt(x*x+y*y+z*z);
}

inline cVector3D cVector3D::crossProd(const cVector3D& v)
{
    cVector3D answer(y*v.z-v.y*z,-x*v.z+v.x*z,x*v.y-v.x*y);
    return answer;
}
inline double cVector3D::dotProd(const cVector3D& inpVec)
{
	return x*inpVec.x+y*inpVec.y+z*inpVec.z;
}

inline cVector3D cVector3D::scale(const double scaling)
{
	*this *= scaling;
	return *this;
}

inline double cVector3D::gamma()
{
	return 1./sqrt(1.-(x*x+y*y+z*z)*(8.99875479107368E+18));  //8.99875479107368E+18 = c^2
}
inline void cVector3D::setVal(const double inpX, const double inpY, const double inpZ)
{
	x=inpX;
	y=inpY;
	z=inpZ;
}

//To solve 3D equations easily (this is not a dot product or cross product
inline const cVector3D cVector3D::multComps(cVector3D inp) const
{
	cVector3D result = *this;
	result.x *= inp.x;
	result.y *= inp.y;
	result.z *= inp.z;
	return result;
}
//To solve 3D equations easily (this is not a dot product or cross product
inline const cVector3D cVector3D::divComps(cVector3D inp) const
{
	cVector3D result = *this;
	result.x /= inp.x;
	result.y /= inp.y;
	result.z /= inp.z;
	return result;
}

inline const cVector3D cVector3D::returnTransCoords() const //rotation of 9.5 degrees
{
	cVector3D result;
	result.x = -z*0.16504760586067764838263370078844+x*0.98628560153723140782496568567775;
	result.z = z*0.98628560153723140782496568567775+x*0.16504760586067764838263370078844;
	result.y = y;
	return result;
}
inline const cVector3D cVector3D::returnNormCoords() const  //rotation back 9.5 degrees
{
	cVector3D result;
	result.x = z*0.16504760586067764838263370078844+x*0.98628560153723140782496568567775;
	result.z = z*0.98628560153723140782496568567775-x*0.16504760586067764838263370078844;
	result.y = y;
	return result;
}

inline cVector3D cVector3D::norm()
{
	cVector3D result;
	double oneOverMag = 1/this->mag();

	result.x = x * oneOverMag;
	result.y = y * oneOverMag;
	result.z = z * oneOverMag;
	return result;
}

inline const double cVector3D::minComp()
{
	double min =std::abs(x);
	if(std::abs(y)<min)
		min=std::abs(y);
	if(std::abs(z)<min)
		min=std::abs(z);
	return min;
}

inline const double cVector3D::maxComp()
{
	double max =std::abs(x);
	if(std::abs(y)>max)
		max=std::abs(y);
	if(std::abs(z)>max)
		max=std::abs(z);
	return max;
}

inline cVector3D zeroVec3D(){cVector3D v(0,0,0);return v;}

////////////////////////////////////////////////////////////////////////////////////////////////////////


inline const cVector2D cVector2D::operator+(const cVector2D& rhs) const
{
	cVector2D result = *this;
	result += rhs;
	return result;
}

inline const cVector2D cVector2D::operator-(const cVector2D& rhs) const
{
	cVector2D result = *this;
	result -= rhs;
	return result;
}

inline const cVector2D& cVector2D::operator=(const cVector2D& inpVec)
{
	r = inpVec.r;
	z = inpVec.z;
	return *this;
}

inline const cVector2D& cVector2D::operator += (const cVector2D& rhs)
{
	r +=rhs.r;
	z +=rhs.z;
	return *this;
}
inline const cVector2D& cVector2D::operator -= (const cVector2D& rhs)
{
	r -=rhs.r;
	z -=rhs.z;
	return *this;
}
inline const cVector2D& cVector2D::operator *= (double rhs)
{
	r *= rhs;
	z *= rhs;
	return *this;
}

//Should throw erception if rhs is 0...but doesn't yet
inline const cVector2D& cVector2D::operator /= (double rhs)
{
	r /=rhs;
	z /=rhs;
	return *this;
}


inline const cVector2D cVector2D::operator * (const double rhs) const
{
	cVector2D result = *this;
	result *= rhs;
	return result;
}

inline cVector2D operator* (const double& lhs, const cVector2D& rhs)
{
	cVector2D result = rhs;
	result *= lhs;
	return result;
}

inline const cVector2D cVector2D::operator / (const double rhs) const
{
	cVector2D result = *this;
	result /= rhs;
	return result;
}

inline const bool cVector2D::operator== (const cVector2D& inpVec)
{
	if(r==inpVec.r && z==inpVec.z){
		return true;
	}
	else{
		return false;
	}
}

inline const bool cVector2D::operator != (const cVector2D& rhs)
{
	return !(*this == rhs);
}

inline double cVector2D::mag()
{
	return sqrt(r*r+z*z);
}

inline double cVector2D::dotProd(const cVector2D& inpVec)
{
	return r*inpVec.r+z*inpVec.z;
}

inline cVector2D cVector2D::scale(const double scaling)
{
	*this *= scaling;
	return *this;
}

inline double cVector2D::gamma()
{
	return 1./sqrt(1.-(r*r+z*z)/(8.99875479107368E+18));  //8.99875479107368E+18 = c^2
}
inline void cVector2D::setVal(const double inpR, const double inpZ)
{
	r=inpR;
	z=inpZ;
}

//To solve 2D equations easily (this is not a dot product or cross product
inline const cVector2D cVector2D::multComps(cVector2D inp) const
{
	cVector2D result = *this;
	result.r *= inp.r;
	result.z *= inp.z;
	return result;
}
//To solve 2D equations easily (this is not a dot product or cross product
inline const cVector2D cVector2D::divComps(cVector2D inp) const
{
	cVector2D result = *this;
	result.r /= inp.r;
	result.z /= inp.z;
	return result;
}
inline cVector2D zeroVec2D(){cVector2D v(0,0); return v;}

#endif // CVECTOR_H_INCLUDED

