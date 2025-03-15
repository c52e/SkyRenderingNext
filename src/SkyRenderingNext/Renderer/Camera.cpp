#include "Camera.h"

#include "imgui.h"

static const Math::float3 kWorldUp{0, 1, 0};

Camera::Camera(Math::float3 position, float yaw, float pitch, float aspect)
    : m_Position(position), m_Up(kWorldUp), m_Yaw(yaw), m_Pitch(pitch), m_Aspect(aspect)
{
    UpdateVectors();
}

Math::float4x4 Camera::ViewMatrix() const
{
    return Math::lookAt(m_Position, m_Position + m_Front * Math::max(Math::length(m_Position), 1.0f), m_Up);
}

static float GetTanHalfFovy(float focalLength)
{
    return 24.0f * 0.5f / focalLength;;
}

Math::float4x4 Camera::ProjectionInvertedInfinity() const
{
    float tanHalfFovy = GetTanHalfFovy(focalLength);
    Math::float4x4 proj{{1.0f / (m_Aspect * tanHalfFovy), 0.0f, 0.0f, 0.0f},
                    {0.0f, 1.0f / tanHalfFovy, 0.0f, 0.0f},
                    {0.0f, 0.0f, 0.0f, -1.0f},
                    {0.0f, 0.0f, m_zNear, 0.0f}};
    proj[2][0] += -jitter.x;
    proj[2][1] += -jitter.y;
    return proj;
}

Math::float4x4 Camera::ProjectionFarClip(float zFar) const
{
    float tanHalfFovy = GetTanHalfFovy(focalLength);
    float zNear = m_zNear;
    Math::float4x4 Result(0.0f);
    Math::float4x4 proj{{1.0f / (m_Aspect * tanHalfFovy), 0.0f, 0.0f, 0.0f},
                    {0.0f, 1.0f / tanHalfFovy, 0.0f, 0.0f},
                    {0.0f, 0.0f, -zFar / (zFar - zNear), -1.0f},
                    {0.0f, 0.0f, -(zFar * zNear) / (zFar - zNear), 0.0f}};
    return proj;
}

float Camera::LensDistance() const
{
    float focalLengthInMeter = focalLength * 0.001f;
    return focalLengthInMeter * focusingDistance / (focusingDistance - focalLengthInMeter);
}

float Camera::ApertureRadius() const
{
    float focalLengthInMeter = focalLength * 0.001f;
    return focalLengthInMeter / fNumber * 0.5f;
}

void Camera::Rotate(float dPitch, float dYaw)
{
    m_Pitch += dPitch;
    m_Yaw += dYaw;

    if (m_Pitch > 89.0f)
        m_Pitch = 89.0f;
    else if (m_Pitch < -89.0f)
        m_Pitch = -89.0f;

    UpdateVectors();
}

bool Camera::DefaultFrameControl()
{
    m_PreViewProjection = ViewProjection();

    float dt = ImGui::GetIO().DeltaTime;

    float dForward = 0.0f;
    float dRight = 0.0f;
    float dUp = 0.0f;
    if (ImGui::IsKeyDown(ImGuiKey_W))
        dForward += movingSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_S))
        dForward -= movingSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_D))
        dRight += movingSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_A))
        dRight -= movingSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_E))
        dUp += movingSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_Q))
        dUp -= movingSpeed * dt;
    auto moveVector = m_Front * dForward + m_Right * dRight + kWorldUp * dUp;
    m_Position += moveVector;

    float rotateSpeed = 180.0f;
    float dPitch = 0.0f;
    float dYaw = 0.0f;
    if (ImGui::IsKeyDown(ImGuiKey_UpArrow))
        dPitch += rotateSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_DownArrow))
        dPitch -= rotateSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_RightArrow))
        dYaw += rotateSpeed * dt;
    if (ImGui::IsKeyDown(ImGuiKey_LeftArrow))
        dYaw -= rotateSpeed * dt;

    if (ImGui::IsMouseDown(ImGuiMouseButton_Right))
    {
        float mouseRotateSpeed = 10.0f;
        dPitch += mouseRotateSpeed * -ImGui::GetIO().MouseDelta.y * dt;
        dYaw += mouseRotateSpeed * ImGui::GetIO().MouseDelta.x * dt;
    }

    Rotate(dPitch, dYaw);

    return Math::length(moveVector) != 0.0f || dPitch != 0 || dYaw != 0;
}

static void FromThetaPhiToDirection(float theta, float phi, float direction[3])
{
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    float cos_phi = cos(phi);
    float sin_phi = sin(phi);
    direction[0] = cos_phi * sin_theta;
    direction[1] = cos_theta;
    direction[2] = sin_phi * sin_theta;
}

void Camera::UpdateVectors()
{
    FromThetaPhiToDirection(Math::radians(90.f - m_Pitch), Math::radians(m_Yaw), Math::value_ptr(m_Front));
    m_Right = Math::normalize(Math::cross(m_Front, kWorldUp));
    m_Up = Math::normalize(Math::cross(m_Right, m_Front));
}
