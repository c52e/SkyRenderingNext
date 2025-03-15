#include "Tlas.h"

Tlas::Tlas(GraphicsCommandList *cmd, std::span<GltfModel *> models, std::span<Math::float4x4> matrices)
{
    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();

    ASSERT(models.size() == matrices.size());

    {
        std::vector<D3D12_RAYTRACING_INSTANCE_DESC> instanceDescs;

        for (size_t i = 0; i < models.size(); ++i)
        {
            auto &model = models[i];
            const auto &matrix = matrices[i];
            model->ForEachInstanceDesc([&](D3D12_GPU_VIRTUAL_ADDRESS blas, Math::float4x4 m) {
                D3D12_RAYTRACING_INSTANCE_DESC instanceDesc{
                    .InstanceMask = 1,
                    .AccelerationStructure = blas,
                };
                m = matrix * m;
                for (int r = 0; r < 3; ++r)
                {
                    for (int c = 0; c < 4; ++c)
                    {
                        instanceDesc.Transform[r][c] = m[c][r];
                    }
                }
                instanceDescs.push_back(instanceDesc);
            });
        }

        auto instanceBufferSize = instanceDescs.size() * sizeof(instanceDescs[0]);
        if (!instanceDescs.empty())
        {
            m_InstanceBuffer = StructuredBuffer::Create(cmd, instanceDescs);
        }

        D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_INPUTS tlasInput{
            .Type = D3D12_RAYTRACING_ACCELERATION_STRUCTURE_TYPE_TOP_LEVEL,
            .Flags = D3D12_RAYTRACING_ACCELERATION_STRUCTURE_BUILD_FLAG_PREFER_FAST_TRACE,
            .NumDescs = static_cast<UINT>(instanceDescs.size()),
            .DescsLayout = D3D12_ELEMENTS_LAYOUT_ARRAY,
            .InstanceDescs = m_InstanceBuffer.GetResource() ? m_InstanceBuffer.GetGPUVirtualAddress() : 0,
        };

        D3D12_RAYTRACING_ACCELERATION_STRUCTURE_PREBUILD_INFO tlasPrebuildInfo{};
        d3d12Device->GetRaytracingAccelerationStructurePrebuildInfo(&tlasInput, &tlasPrebuildInfo);

        m_TlasScratch = GpuBuffer(tlasPrebuildInfo.ScratchDataSizeInBytes, D3D12_RESOURCE_STATE_COMMON,
                                  D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);

        m_Tlas =
            GpuBuffer(tlasPrebuildInfo.ResultDataMaxSizeInBytes, D3D12_RESOURCE_STATE_RAYTRACING_ACCELERATION_STRUCTURE,
                      D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);

        m_RTASHandleIndex = device.AllocBindlessDescriptor();
        D3D12_SHADER_RESOURCE_VIEW_DESC rtasDesc{
            .ViewDimension = D3D12_SRV_DIMENSION_RAYTRACING_ACCELERATION_STRUCTURE,
            .Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING,
            .RaytracingAccelerationStructure{.Location = m_Tlas.GetGPUVirtualAddress()},
        };
        d3d12Device->CreateShaderResourceView(nullptr, &rtasDesc, device.GetBindlessCpuHandle(m_RTASHandleIndex));

        D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_DESC tlasBuildDesc{
            .DestAccelerationStructureData = m_Tlas.GetGPUVirtualAddress(),
            .Inputs = tlasInput,
            .ScratchAccelerationStructureData = m_TlasScratch.GetGPUVirtualAddress(),
        };

        cmd->BuildRaytracingAccelerationStructure(&tlasBuildDesc, 0, nullptr);
        m_Tlas.BarrierUAV(cmd);
    }

    std::vector<int> materialIndexOffsets;
    {
        std::vector<ShaderResource::Material> materialBuffer;
        for (auto &model : models)
        {
            materialIndexOffsets.push_back(static_cast<int>(materialBuffer.size()));
            model->ForEachMaterial([&](const ShaderResource::Material &m) { materialBuffer.push_back(m); });
        }
        if (!materialBuffer.empty())
        {
            m_MaterialBuffer = StructuredBuffer::Create(cmd, materialBuffer);
        }
    }

    {
        std::vector<UINT> instanceOffsetBuffer;
        std::vector<ShaderResource::Geometry> geometryBuffer;

        UINT modelIndex = 0;
        for (auto &model : models)
        {
            model->ForEachInstanceAndGeometry(
                [&]() { instanceOffsetBuffer.push_back(static_cast<UINT>(geometryBuffer.size())); },
                [&](const GltfModel::Geometry &g) {
                    geometryBuffer.push_back({modelIndex, g.materialIndex + materialIndexOffsets[modelIndex],
                                              g.vertexOffset, g.indexOffset});
                });
            modelIndex++;
        }
        if (!instanceOffsetBuffer.empty() && !geometryBuffer.empty())
        {
            m_InstanceOffsetBuffer = StructuredBuffer::Create(cmd, instanceOffsetBuffer);
            m_GeometryBuffer = StructuredBuffer::Create(cmd, geometryBuffer);
        }
    }

    {
        std::vector<ShaderResource::Model> modelBuffer;
        for (const auto &model : models)
        {
            modelBuffer.push_back({
                model->GetIndexBuffer().GetSrvIndex(),
                model->GetVertexBuffer().GetSrvIndex(),
            });
        }
        if (!modelBuffer.empty())
        {
            m_ModelBuffer = StructuredBuffer::Create(cmd, modelBuffer);
        }
    }
        
    m_RTASRootConstant = {
        .tlas = m_RTASHandleIndex.index(),
        .modelBuffer = m_ModelBuffer.GetSrvIndex(),
        .materialBuffer = m_MaterialBuffer.GetSrvIndex(),
        .instanceOffsetBuffer = m_InstanceOffsetBuffer.GetSrvIndex(),
        .geometryBuffer = m_GeometryBuffer.GetSrvIndex(),
    };
}
