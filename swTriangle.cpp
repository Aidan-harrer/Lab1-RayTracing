#include "swTriangle.h"

namespace sw {

bool Triangle::intersect(const Ray &r, Intersection &isect) const {
    // TODO: Implement ray-triangle intersection
    Vec3 e1 = vertices[1] - vertices[0];
    Vec3 e2 = vertices[2] - vertices[0];
    Vec3 n = e1 % e2;
    float m = -n * vertices[0];
    float t = (n * r.origin() + m) / (-n * r.dir);
    if (r.minT > t || t > r.maxT) return false;
    Vec3 Q = r.origin() + t * r.dir;
    Vec3 r_plane = Q - vertices[0];
    float v = (e1%r_plane).norm() / n.norm();
    float w = (r_plane%e2).norm() / n.norm();
    if (v < 0) return false;
    if (w < 0) return false;
    if (v + w >= 1) return false;
    isect.hitT = t;
    isect.normal = n;
    isect.normal.normalize();
    isect.frontFacing = (-r.dir * isect.normal) > 0.0f;
    if (!isect.frontFacing) isect.normal = -isect.normal;
    isect.position = Q;
    isect.material = material;
    isect.ray = r;
    return true;
}

} // namespace sw
