/*
 *  main.cpp
 *  swTracer
 *
 *  Created by Michael Doggett on 2021-09-23.
 *  Copyright (c) 2021 Michael Doggett
 */
#define _USE_MATH_DEFINES
#include <cfloat>
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "swCamera.h"
#include "swIntersection.h"
#include "swMaterial.h"
#include "swRay.h"
#include "swScene.h"
#include "swSphere.h"
#include "swVec3.h"

using namespace sw;

inline float clamp(float x, float min, float max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

float uniform() {
    // Will be used to obtain a seed for the random number engine
    static std::random_device rd;
    // Standard mersenne_twister_engine seeded with rd()
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    return dis(gen);
}

void writeColor(int index, Vec3 p, uint8_t *pixels) {
    // gamma correct for gamma=2.2, x^(1/gamma), more see :
    // https://www.geeks3d.com/20101001/tutorial-gamma-correction-a-story-of-linearity/
    for (int n = 0; n < 3; n++) {
        p.m[n] = pow(p.m[n], 1.0f / 2.2f);
        pixels[index + n] = (uint8_t)(256 * clamp(p.m[n], 0.0f, 0.999f));
    }
}

Vec3 tileNormal(const Vec3 &normal, float tileSize) {
    // Convert normal vector to spherical coordinates
    float u = 0.5f + atan2(normal.m[2], normal.m[0]) / (2.0f * M_PI); // z -> m[2], x -> m[0]
    float v = 0.5f - asin(normal.m[1]) / M_PI;                        // y -> m[1]

    // Map (u, v) to discrete tiles
    u = floor(u / tileSize) * tileSize;
    v = floor(v / tileSize) * tileSize;

    // Convert back to a perturbed normal
    float theta = 2.0f * M_PI * u;
    float phi = M_PI * (1.0f - v);

    Vec3 perturbedNormal;
    perturbedNormal.m[0] = sin(phi) * cos(theta); // x -> m[0]
    perturbedNormal.m[1] = cos(phi);              // y -> m[1]
    perturbedNormal.m[2] = sin(phi) * sin(theta); // z -> m[2]

    return perturbedNormal;
}
Vec3 reflect(const Vec3 &lightDir, const Vec3 &normal) { return 2.0f * (normal * lightDir) * normal - lightDir; }
Vec3 clampColor(const Vec3 &color, float minVal, float maxVal) {
    return Vec3(std::max(minVal, std::min(color.x(), maxVal)), std::max(minVal, std::min(color.y(), maxVal)),
                std::max(minVal, std::min(color.z(), maxVal)));
}

Color calculateLighting(const Vec3 &lightDir, const Vec3 &normal, const Vec3 &viewDir, const Color &lightColor) {
    Vec3 reflectDir = reflect(lightDir, normal);

    // Diffuse shading
    float diffuse = std::max(normal * lightDir, 0.0f);

    // Specular highlights with sharpness control
    float shininess = 16.0f; // Adjust for sharper or softer highlights
    float specular = std::pow(std::max(viewDir * reflectDir, 0.0f), shininess);

    // Combine contributions
    Color lightContribution = lightColor * (diffuse + 0.5f * specular);

    // Clamp light contributions to avoid overexposure
    return clampColor(lightContribution, 0.0f, 1.0f);
}
// Function to generate random directions for scattering rays
Vec3 randomInUnitSphere() {
    Vec3 p;
    do {
        p = Vec3(uniform(), uniform(), uniform()) * 2.0f - Vec3(1.0f, 1.0f, 1.0f); // Random point in cube
    } while (p.lengthSquared() >= 1.0f); // Reject points outside unit sphere
    return p;
}

Color traceRay(const Ray &r, Scene scene, int depth) {
    Color c, directColor, reflectedColor, refractedColor;
    if (depth < 0) return c; // Base case for recursion

    Intersection hit, shadow;
    if (!scene.intersect(r, hit)) return Color(0.0f, 0.0f, 0.0f); // Background color

    // Define 4 small lights inside the disco ball
    const Vec3 lightPositions[] = {//   Vec3(0.7f, 15.5f, -10.2f),
                                   //   Vec3(-0.7f, 15.5f, -9.5f),
                                   //   Vec3(0.5f, 15.5f, -10.5f),
                                   //   Vec3(-0.5f, 15.5f, -9.5f),
                                   // };
                                   Vec3(0.0f, 20.0f, -5.0f)};
    const int numLights = 1;
    const float lightSize = 0.1f; // Small area light size

    Vec3 perturbedNormal = tileNormal(hit.normal, 1.0f); // Adjust `tileSize` to control facet size

    directColor = Color(0.0f, 0.0f, 0.0f); // Accumulate direct lighting from all lights

    // Calculate direct lighting from each light source
    for (int l = 0; l < numLights; ++l) {
        Vec3 lightPos = lightPositions[l];
        Vec3 lightDir = lightPos - hit.position;
        lightDir.normalize();

        // Soft shadows with area light sampling
        int numSamples = 32;
        float shadowIntensity = 0.0f;

        for (int i = 0; i < numSamples; ++i) {
            Vec3 randomPointOnLight =
              lightPos + Vec3((uniform() - 0.5f) * lightSize, 0.0f, (uniform() - 0.5f) * lightSize);
            if (scene.intersect(hit.getShadowRay(randomPointOnLight), shadow)) {
                shadowIntensity += 0.75f;
            }
        }
        // Check if the material is reflective (for disco ball)
        if (hit.material.id == 1) {
            reflectedColor = Color(0.0f, 1.0f, 0.0f);

            // Generate multiple reflection rays for the disco ball
            const int numReflectionRays = 5; // Number of rays to generate
            for (int i = 0; i < numReflectionRays; ++i) {
                // Randomly scatter the reflection rays around the hit point
                Vec3 scatterDirection = hit.normal + randomInUnitSphere();
                scatterDirection.normalize();

                // Create a new ray from the hit point in the scattered direction
                Ray reflectionRay(hit.position, scatterDirection);

                // Trace the scattered ray
                reflectedColor += traceRay(reflectionRay, scene, depth - 1);
            }
        }
        // Handle refractions
        if (hit.material.transparency > 0) {
            refractedColor = traceRay(hit.getRefractedRay(), scene, depth - 1);
            float transparencyFactor = shadow.material.transparency > 0 ? 0.5f : 1.0f;
            shadowIntensity += 0.75f * transparencyFactor;
        }
        shadowIntensity /= numSamples;

        // Calculate direct lighting for this light
        float specularBoost = 2.0f; // Boost for reflected light on the disco ball
        Color lightContribution = hit.material.color * std::max(hit.normal * lightDir, 0.0f);

        // Simulate specular effect
        Vec3 reflectDir = Vec3.normalize(glm::reflect(-lightDir, hit.normal));
        float specFactor =
          pow(std::max(glm::dot(reflectDir, viewDir), 0.0f), 16); // Higher exponent sharpens the highlight
        Color specularHighlight = lightColor * specFactor * specularBoost;

        lightContribution += specularHighlight;

        lightContribution *= (1.0f - shadowIntensity);
        lightContribution.r = std::min(lightContribution.r, 1.0f);
        lightContribution.g = std::min(lightContribution.g, 1.0f);
        lightContribution.b = std::min(lightContribution.b, 1.0f);

        directColor += lightContribution; // Add contribution from this light
    }

    // Combine lighting contributions
    c = (1.0f - hit.material.reflectivity - hit.material.transparency) * directColor +
        hit.material.reflectivity * reflectedColor + hit.material.transparency * refractedColor;

    return c;
}

int main() {
    const int imageWidth = 512;
    const int imageHeight = imageWidth;
    const int numChannels = 3;
    uint8_t *pixels = new uint8_t[imageWidth * imageHeight * numChannels];
    // Define materials
    Material whiteDiffuse = Material(Color(0.9f, 0.9f, 0.9f), 0.0f, 0.0f, 1.0f);
    Material greenDiffuse = Material(Color(0.1f, 0.6f, 0.1f), 0.0f, 0.0f, 1.0f);
    Material redDiffuse = Material(Color(1.0f, 0.1f, 0.1f), 0.0f, 0.0f, 1.0f);
    Material blueDiffuse = Material(Color(0.0f, 0.2f, 0.9f), 0.0f, 0.0f, 1.0f);
    Material yellowReflective = Material(Color(1.0f, 0.6f, 0.1f), 0.2f, 0.0f, 1.0f);
    Material transparent = Material(Color(1.0f, 1.0f, 1.0f), 0.2f, 0.8f, 1.3f);

    // Setup scene
    Scene scene;

    // Add three spheres with diffuse material
    Material discoBallMaterial = Material(Color(0.5f, 0.5f, 0.5f), 0.9f, 0.0f, 1.5f, 1);

    // Add the disco ball to the scene
    scene.push(Sphere(Vec3(0.0f, 15.0f, -10.0f), 2.0f, discoBallMaterial));

    scene.push(Sphere(Vec3(-3.0f, 9.0f, -3.0f), 2.0f, transparent));
    scene.push(Sphere(Vec3(3.0f, 9.0f, -3.0f), 2.0f, transparent));

    // Define vertices for Cornell box
    Vec3 vertices[] = {
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(20.0f, 0.0f, 50.0f),    Vec3(20.0f, 0.0f, -50.0f),   // Floor 1
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(20.0f, 0.0f, -50.0f),   Vec3(-20.0f, 0.0f, -50.0f),  // Floor 2
      Vec3(-20.0f, 0.0f, -50.0f), Vec3(20.0f, 0.0f, -50.0f),   Vec3(20.0f, 40.0f, -50.0f),  // Back wall 1
      Vec3(-20.0f, 0.0f, -50.0f), Vec3(20.0f, 40.0f, -50.0f),  Vec3(-20.0f, 40.0f, -50.0f), // Back wall 2
      Vec3(-20.0f, 40.0f, 50.0f), Vec3(-20.0f, 40.0f, -50.0f), Vec3(20.0f, 40.0f, 50.0f),   // Ceiling 1
      Vec3(20.0f, 40.0f, 50.0f),  Vec3(-20.0f, 40.0f, -50.0f), Vec3(20.0f, 40.0f, -50.0f),  // Ceiling 2
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(-20.0f, 40.0f, -50.0f), Vec3(-20.0f, 40.0f, 50.0f),  // Red wall 1
      Vec3(-20.0f, 0.0f, 50.0f),  Vec3(-20.0f, 0.0f, -50.0f),  Vec3(-20.0f, 40.0f, -50.0f), // Red wall 2
      Vec3(20.0f, 0.0f, 50.0f),   Vec3(20.0f, 40.0f, -50.0f),  Vec3(20.0f, 40.0f, 50.0f),   // Green wall 1
      Vec3(20.0f, 0.0f, 50.0f),   Vec3(20.0f, 0.0f, -50.0f),   Vec3(20.0f, 40.0f, -50.0f)   // Green wall 2
    };

    // TODO: Uncomment to render floor triangles
    scene.push(Triangle(&vertices[0], whiteDiffuse)); // Floor 1
    scene.push(Triangle(&vertices[3], whiteDiffuse)); // Floor 2

    // TODO: Uncomment to render Cornell box
    scene.push(Triangle(&vertices[6], whiteDiffuse));  // Back wall 1
    scene.push(Triangle(&vertices[9], whiteDiffuse));  // Back wall 2
    scene.push(Triangle(&vertices[12], whiteDiffuse)); // Ceiling 1
    scene.push(Triangle(&vertices[15], whiteDiffuse)); // Ceiling 2
    scene.push(Triangle(&vertices[18], redDiffuse));   // Red wall 1
    scene.push(Triangle(&vertices[21], redDiffuse));   // Red wall 2
    scene.push(Triangle(&vertices[24], greenDiffuse)); // Green wall 1
    scene.push(Triangle(&vertices[27], greenDiffuse)); // Green wall 2

    // TODO: Uncomment to render reflective spheres

    // Setup camera
    Vec3 eye(0.0f, 10.0f, 30.0f);
    Vec3 lookAt(0.0f, 10.0f, -5.0f);
    Vec3 up(0.0f, 1.0f, 0.0f);
    Camera camera(eye, lookAt, up, 52.0f, (float)imageWidth / (float)imageHeight);
    camera.setup(imageWidth, imageHeight);
    // Ray trace pixels
    int depth = 3;
    std::cout << "Rendering... ";
    clock_t start = clock();
    for (int j = 0; j < imageHeight; ++j) {
        for (int i = 0; i < imageWidth; ++i) {

            Color pixel;

            // Get center of pixel coordinate
            for (int x = 0; x < 3; ++x) {
                for (int y = 0; y < 3; ++y) {

                    float cx = (i) + (x) * (uniform()) / 3.0f;

                    float cy = (j) + (y) * (uniform()) / 3.0f;

                    // Get a ray and trace it
                    Ray r = camera.getRay(cx, cy);
                    pixel += traceRay(r, scene, depth);
                }
            }
            pixel *= (1.0f / 9.0f);
            // Write pixel value to image
            writeColor((j * imageWidth + i) * numChannels, pixel, pixels);
        }
    }

    // Save image to file
    stbi_write_png("out.png", imageWidth, imageHeight, numChannels, pixels, imageWidth * numChannels);

    // Free allocated memory
    delete[] pixels;

    std::cout << "Done\n";
    std::cout << "Time: " << (float)(clock() - start) / CLOCKS_PER_SEC << " s" << std::endl;
}