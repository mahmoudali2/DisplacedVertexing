/*
 * Region selection for a vertex, relative to a cylinder of radius R and half
 * length Z centred on the IP (for example the solenoid volume, or a tracker).
 *
 *   outside : r_xy > R  or  |z| > Z
 *   barrel  : r_xy > R
 *   endcap  : r_xy <= R and |z| > Z
 *   inside  : r_xy <= R and |z| <= Z
 *   any     : no requirement
 *
 * Units: mm.
 * Author: Mahmoud Althakeel
 */
#pragma once

#include <TVector3.h>

#include <cmath>
#include <string>

namespace k4vertexing {

class VertexRegion {
public:
  /// returns false if the region name is not one of the five above
  bool configure(const std::string& name, double radius, double halfZ) {
    m_name = name;
    m_radius = radius;
    m_halfZ = halfZ;
    return name == "outside" || name == "barrel" || name == "endcap" || name == "inside" || name == "any";
  }

  bool contains(const TVector3& v) const {
    // where is the point relative to the cylinder?
    const bool beyondRadius = v.Perp() > m_radius;
    const bool beyondHalfZ = std::fabs(v.Z()) > m_halfZ;

    if (m_name == "outside")
      return beyondRadius || beyondHalfZ;
    if (m_name == "barrel")
      return beyondRadius;
    if (m_name == "endcap")
      return !beyondRadius && beyondHalfZ;
    if (m_name == "inside")
      return !beyondRadius && !beyondHalfZ;
    return true; // "any"
  }

private:
  std::string m_name = "any";
  double m_radius = 0;
  double m_halfZ = 0;
};

} // namespace k4vertexing
