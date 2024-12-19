#pragma once

#include "swVec3.h"

namespace sw {

class Material {
  public:
    Material() = default;

    // Constructor with specular added
    Material(const Vec3 &c, float r = 0, float t = 0, float i = 1, int I = 0, float s = 0)
      : color{c}, reflectivity(r), transparency(t), refractiveIndex(i), id(I), specular(s) {}

  public:
    Vec3 color;                  // Base color of the material
    float reflectivity{0};       // Reflectivity of the material (0-1)
    float transparency{0};       // Transparency (0-1)
    float refractiveIndex{1.0f}; // Refractive index for transparent materials
    float specular{0};           // Specular highlight intensity (0-1)
    int id{1};                   // Material identifier
};

} // namespace sw
