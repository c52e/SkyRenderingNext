#pragma once

#include <optional>

#include "Core/Common.h"
#include "Core/Serialization.h"

class Camera
{
  public:
    Camera() = default;

    Camera(Math::float3 position, float yaw, float pitch, float aspect);

    Math::float4x4 ViewMatrix() const;

    Math::float4x4 ProjectionInvertedInfinity() const;

    Math::float4x4 ProjectionFarClip(float zFar) const;

    Math::float4x4 ViewProjection() const
    {
        return ProjectionInvertedInfinity() * ViewMatrix();
    }

    Math::float4x4 PreViewProjection() const
    {
        return m_PreViewProjection.value_or(ViewProjection());
    }

    Math::float3 GetFront() const
    {
        return m_Front;
    }

    float LensDistance() const;

    float ApertureRadius() const;

    void Rotate(float dPitch, float dYaw);

    // Return true if camera view changed
    bool DefaultFrameControl();

    Math::float3 m_Position{};
    float m_zNear = 1e-1f;

    float m_Yaw{};
    float m_Pitch{};
    float m_Aspect{1.0f};
    Math::float2 jitter{};

    // https://borislavkostov.wordpress.com/articles/macro-photography-eng/focal-length-focusing-distance-working-distance/
    float focalLength = 50.0f;     // mm
    float focusingDistance = 2.0f; // m
    float fNumber = 8.0f;          // f/D

    float movingSpeed = 6.0f;

  private:
    void UpdateVectors();

    Math::float3 m_Front{};
    Math::float3 m_Right{};
    Math::float3 m_Up{};

    std::optional<Math::float4x4> m_PreViewProjection{};
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Camera
    , m_Position
    , m_zNear
    , m_Yaw
    , m_Pitch
    , focalLength
    , focusingDistance
    , fNumber
    , movingSpeed)
