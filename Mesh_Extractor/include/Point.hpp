/*!
 * \file Point.hpp
 * \brief Simple class representing a physical point
 * \author J. Li Volsi
*/

#ifndef __POINT_HPP__
#define  __POINT_HPP__
// __INCLUDE__
#include <vector>
#include <iostream>

namespace Geometry
{
  //! A simple extendable class that represents a point.
  /*! The container chosen for the class is simply a vector
      of double initialized by default to the 3 dimensional
      zero vector */

  // template<typename T> possible extension instead of double
  class Point  {
  private:
    int Dim;
    std::vector<double> coord;
  public:

    Point(int N=3,std::vector<double> v={0.0,0.0,0.0}): Dim(N), coord(v) {}; // full constructor
    //Point(std::vector<double> v): coord(v) {Dim = v.size();}; // constructor only coord vector
    Point(std::vector<double> v): coord(v) {this->set_dim(v.size());}; // constructor only coord vector
    Point(Point const &) = default; // default copy constructor - necessary?

    double operator[](std::size_t i) const;

    double & operator[](std::size_t i);

    bool operator==(const Point& p) const;

    bool operator>(const Point& p) const;

    bool operator>=(const Point& p) const;

    bool operator<(const Point& p) const;

    int get_dim() const { return Dim;};

    double get_x() const;

    double get_y() const;

    double get_z() const;

    void set_dim(int d) { Dim = d;};

    void set_coord(std::vector<double> v);

    void print_point() const;

    void print_x() const;

    void print_y() const;

    void print_z() const;


}; // end class

}// end namespace
#endif
