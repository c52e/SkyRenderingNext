#include "RenderBuffer.h"

BufferBase::BufferBase(UINT64 size, D3D12_RESOURCE_STATES initialState, D3D12_RESOURCE_FLAGS flags,
                       D3D12_HEAP_TYPE type)
    : ResourceBarrierBase(initialState), m_Size(size)
{
    CD3DX12_HEAP_PROPERTIES heapProps(type);
    auto desc = CD3DX12_RESOURCE_DESC::Buffer(size, flags);
    GetGraphicsDevice().GetD3D12Device()->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &desc, initialState,
                                                                  nullptr, IID_PPV_ARGS(&m_Resource));
}

void UploadBuffer::UploadData(const void *data, UINT64 size)
{
    void *basePtr{};
    m_Resource->Map(0, nullptr, &basePtr);
    memcpy(basePtr, data, size);
    m_Resource->Unmap(0, nullptr);
}

ConstantBuffer::ConstantBuffer(UINT size) : GpuBuffer(size, D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_FLAG_NONE)
{
    auto &device = GetGraphicsDevice();
    m_Srv = device.AllocBindlessDescriptor();
    D3D12_CONSTANT_BUFFER_VIEW_DESC desc{
        .BufferLocation = GetGPUVirtualAddress(),
        .SizeInBytes = size,
    };
    device.GetD3D12Device()->CreateConstantBufferView(&desc, device.GetBindlessCpuHandle(m_Srv));
}

StructuredBuffer::StructuredBuffer(UINT numElements, UINT structureByteStride, bool allowUav /*= false*/)
    : GpuBuffer(numElements * structureByteStride, D3D12_RESOURCE_STATE_COMMON,
                allowUav ? D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS : D3D12_RESOURCE_FLAG_NONE)
{
    auto &device = GetGraphicsDevice();
    {
        m_Srv = device.AllocBindlessDescriptor();
        D3D12_SHADER_RESOURCE_VIEW_DESC desc{
            .Format = DXGI_FORMAT_UNKNOWN,
            .ViewDimension = D3D12_SRV_DIMENSION_BUFFER,
            .Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING,
            .Buffer{
                .FirstElement = 0,
                .NumElements = numElements,
                .StructureByteStride = structureByteStride,
                .Flags = D3D12_BUFFER_SRV_FLAG_NONE,
            },
        };
        device.GetD3D12Device()->CreateShaderResourceView(GetResource(), &desc, device.GetBindlessCpuHandle(m_Srv));
    }

    if (allowUav)
    {
        m_Uav = device.AllocBindlessDescriptor();
        D3D12_UNORDERED_ACCESS_VIEW_DESC desc{
            .Format = DXGI_FORMAT_UNKNOWN,
            .ViewDimension = D3D12_UAV_DIMENSION_BUFFER,
            .Buffer{
                .FirstElement = 0,
                .NumElements = numElements,
                .StructureByteStride = structureByteStride,
                .CounterOffsetInBytes = 0,
                .Flags = D3D12_BUFFER_UAV_FLAG_NONE,
            },
        };
        device.GetD3D12Device()->CreateUnorderedAccessView(GetResource(), nullptr, &desc,
                                                           device.GetBindlessCpuHandle(m_Uav));
    }
}

ReadbackBuffer::ReadbackBuffer(UINT64 size)
    : BufferBase(size, D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_FLAG_NONE, D3D12_HEAP_TYPE_READBACK)
{
    m_Resource->Map(0, nullptr, &m_Ptr);
}

void ReadbackBuffer::ReadbackFrom(GraphicsCommandList *cmd, StructuredBuffer &from)
{
    m_StartingFrameIndex = GetGraphicsDevice().GetFrameIndex();
    from.BarrierTransitionToCopySrc(cmd);
    cmd->CopyBufferRegion(GetResource(), 0, from.GetResource(), 0, m_Size);
}

void ReadbackBuffer::ReadbackFrom(GraphicsCommandList *cmd, RenderTexture &tex,
                              const D3D12_PLACED_SUBRESOURCE_FOOTPRINT &footprint)
{
    m_StartingFrameIndex = GetGraphicsDevice().GetFrameIndex();
    tex.BarrierTransitionToCopySrc(cmd);

    CD3DX12_TEXTURE_COPY_LOCATION src(tex.GetResource(), 0);
    CD3DX12_TEXTURE_COPY_LOCATION dst(GetResource(), footprint);

    cmd->CopyTextureRegion(&dst, 0, 0, 0, &src, nullptr);
}

bool ReadbackBuffer::Available() const
{
    return GetGraphicsDevice().GetCompleteFrameIndex() >= m_StartingFrameIndex;
}

ReadbackBuffer ReadbackRenderTexture(GraphicsCommandList *cmd, RenderTexture &tex,
                                     D3D12_PLACED_SUBRESOURCE_FOOTPRINT &footprint)
{
    UINT64 totalSize{};
    auto device = GetGraphicsDevice().GetD3D12Device();
    device->GetCopyableFootprints(&tex.GetDesc(), 0, 1, 0, &footprint, nullptr, nullptr, &totalSize);

    ReadbackBuffer buffer(totalSize);
    buffer.ReadbackFrom(cmd, tex, footprint);
    return buffer;
}
