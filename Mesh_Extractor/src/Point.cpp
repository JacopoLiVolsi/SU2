/*!
 * \file Point.cpp
 * \brief Simple class representing a physical point
 * \author J. Li Volsi
*/

// __INCLUDE__
#include <vector>
#include <iostream>
#include <sstream>
// __INCLUDE.HPP__
#include "../include/Point.hpp"

namespace Geometry
{

double Point::operator[](std::size_t i) const
{
  return coord[i];
}

double & Point::operator[](std::size_t i)
{
  return coord[i];
}

bool Point::operator==(const Point& p) const{
  bool chk = true;
  if (p.get_dim() != Dim){
    std::cerr << "Point dimension not comparable" << '\n';
    return false;
  }
  else{
  for (size_t i = 0; i < coord.size(); i++) {
    if (p[i] != coord[i]) {chk = false;}
  }
}
    return chk;
}


bool Point::operator>(const Point& p) const{
  bool chk = true;
  if (p.get_dim() != Dim){
    std::cerr << "Point dimension not comparable" << '\n';
    return false;
  }
  else{
  for (size_t i = 0; i < coord.size(); i++) {
    if (coord[i] < p[i]) {chk = false;}
  }
}
    return chk;
}

bool Point::operator>=(const Point& p) const{
  bool chk = true;
  if (p.get_dim() != Dim){
    std::cerr << "Point dimension not comparable" << '\n';
    return false;
  }
  else{
  for (size_t i = 0; i < coord.size(); i++) {
    if (coord[i] <= p[i]) {chk = false;}
  }
}
    return chk;
}

bool Point::operator<(const Point& p) const{
  return !(*this>p);
}

double Point::get_x() const{
  return coord[0];
}

double Point::get_y() const{
  return coord[1];
}

double Point::get_z() const{
  if(Dim == 3)
    return coord[2];
  else{
    std::cerr << "No z component: returning 0.0" << '\n';
    return 0.0;
  }
}


void Point::set_coord(std::vector<double> v){
  for (size_t i = 0; i < v.size(); i++)
    coord[i] = v[i];
  Dim = v.size();
}

void Point::print_point() const{
  std::cout << "(";
  std::cout << coord[0] << ", " << coord[1];
  if (Dim == 3)
    std::cout << ", " << coord[2];
  std::cout << ")";
  return;
}

void Point::print_x() const{
  std::cout << coord[0] << ' ';
  return;
}

void Point::print_y() const{
  std::cout << coord[1] << ' ';
  return;
}

void Point::print_z() const{
  if( Dim == 3 )
  std::cout << coord[2] << ' ';
  return;
}



} // end namespace
