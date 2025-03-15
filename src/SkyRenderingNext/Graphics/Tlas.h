#pragma once

#include "GltfModel.h"

class Tlas
{
  public:
    Tlas() = default;

    Tlas(GraphicsCommandList *cmd, std::span<GltfModel*> models, std::span<Math::float4x4> matrices);

    const DescriptorHandle &GetHandle() const
    {
        return m_RTASHandleIndex;
    }

    const ShaderResource::RTASRootConstant &GetRootConstant() const
    {
        return m_RTASRootConstant;
    }

  private:
    GpuBuffer m_Tlas;
    GpuBuffer m_TlasScratch;
    DescriptorHandle m_RTASHandleIndex{};
    ShaderResource::RTASRootConstant m_RTASRootConstant{};

    StructuredBuffer m_InstanceBuffer;
    StructuredBuffer m_MaterialBuffer;
    StructuredBuffer m_InstanceOffsetBuffer;
    StructuredBuffer m_GeometryBuffer;
    StructuredBuffer m_ModelBuffer;
};
