#pragma once

#include <span>

#include "GraphicsDevice.h"
#include "ResourceBarrier.h"
#include "RenderTexture.h"

class BufferBase : public ResourceBarrierBase<BufferBase>
{
  public:
    BufferBase() = default;
    BufferBase(BufferBase &&) = default;
    BufferBase &operator=(BufferBase bufferBase)
    {
        m_Resource = std::move(bufferBase.m_Resource);
        m_Size = bufferBase.m_Size;
        return *this;
    }

    BufferBase(UINT64 size, D3D12_RESOURCE_STATES initialState, D3D12_RESOURCE_FLAGS flags, D3D12_HEAP_TYPE type);

    ~BufferBase() = default;

    ID3D12Resource *GetResource()
    {
        return m_Resource.Get();
    }

    UINT64 size() const
    {
        return m_Size;
    }

    D3D12_GPU_VIRTUAL_ADDRESS GetGPUVirtualAddress()
    {
        return m_Resource->GetGPUVirtualAddress();
    }

    D3D12_GPU_VIRTUAL_ADDRESS GetGPUVirtualAddress(UINT64 offset)
    {
        return GetGPUVirtualAddress() + offset;
    }

    operator bool() const noexcept
    {
        return m_Resource;
    }

  protected:
    ResPtr<ID3D12Resource> m_Resource{};
    UINT64 m_Size{};
};

class UploadBuffer : public BufferBase
{
  public:
    UploadBuffer() = default;

    UploadBuffer(UINT64 size)
        : BufferBase(size, D3D12_RESOURCE_STATE_GENERIC_READ, D3D12_RESOURCE_FLAG_NONE, D3D12_HEAP_TYPE_UPLOAD)
    {
    }

    UploadBuffer(const void *data, UINT64 size) : UploadBuffer(size)
    {
        UploadData(data, size);
    }

    void UploadData(const void *data, UINT64 size);
};

class GpuBuffer : public BufferBase
{
  public:
    GpuBuffer() = default;

    GpuBuffer(UINT64 size, D3D12_RESOURCE_STATES initialState, D3D12_RESOURCE_FLAGS flags)
        : BufferBase(size, initialState, flags, D3D12_HEAP_TYPE_DEFAULT)
    {
    }

    const DescriptorHandle &GetSrv() const noexcept
    {
        return m_Srv;
    }

    UINT GetSrvIndex() const noexcept
    {
        return GetSrv().index();
    }

  protected:
    DescriptorHandle m_Srv{};
};

class ConstantBuffer : public GpuBuffer
{
  public:
    ConstantBuffer() = default;

    ConstantBuffer(UINT size);
};

class StructuredBuffer : public GpuBuffer
{
  public:
    StructuredBuffer() = default;

    StructuredBuffer(UINT numElements, UINT structureByteStride, bool allowUav = false);

    const DescriptorHandle &GetUav() const noexcept
    {
        return m_Uav;
    }

    UINT GetUavIndex() const noexcept
    {
        return GetUav().index();
    }

    template <class T> static StructuredBuffer Create(GraphicsCommandList *cmd, std::span<T> data)
    {
        auto bufferSize = data.size() * sizeof(data[0]);
        auto staging = UploadBuffer(data.data(), bufferSize);
        auto res = StructuredBuffer(static_cast<UINT>(data.size()), sizeof(data[0]));

        cmd->CopyBufferRegion(res.GetResource(), 0, staging.GetResource(), 0, bufferSize);
        res.BarrierTransition(cmd, D3D12_RESOURCE_STATE_COPY_DEST,
                              D3D12_RESOURCE_STATE_NON_PIXEL_SHADER_RESOURCE |
                                  D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
        return res;
    }

    template <class T> static StructuredBuffer Create(GraphicsCommandList *cmd, const std::vector<T> &vector)
    {
        return Create(cmd, std::span(vector));
    }

  protected:
    DescriptorHandle m_Uav{};
};

class ReadbackBuffer : public BufferBase
{
  public:
    ReadbackBuffer() = default;

    ReadbackBuffer(UINT64 size);

    void ReadbackFrom(GraphicsCommandList *cmd, StructuredBuffer &from);

    void ReadbackFrom(GraphicsCommandList *cmd, RenderTexture &tex, const D3D12_PLACED_SUBRESOURCE_FOOTPRINT& footprint);

    void *GetPtr() const
    {
        return m_Ptr;
    }

    bool Available() const;

  private:
    void *m_Ptr{};
    UINT64 m_StartingFrameIndex{};
};

ReadbackBuffer ReadbackRenderTexture(GraphicsCommandList *cmd, RenderTexture &tex, D3D12_PLACED_SUBRESOURCE_FOOTPRINT& footprint);
