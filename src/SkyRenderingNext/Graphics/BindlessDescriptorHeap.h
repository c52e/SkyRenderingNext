#pragma once

#include "Core/Common.h"
#include "CommandList.h"

class BindlessDescriptorHeap
{
  public:
    struct Handle
    {
        bool Empty() const noexcept
        {
            return index == static_cast<UINT>(-1);
        }

        UINT index{static_cast<UINT>(-1)};
    };

    void Init(ID3D12Device *device);

    void SetFrame(GraphicsCommandList *cmd);

    ID3D12RootSignature *GetRootSignature() const
    {
        return m_RootSignature.Get();
    }

    Handle Alloc();

    D3D12_CPU_DESCRIPTOR_HANDLE GetCpuHandle(Handle handle);

    void Free(Handle handle);

  private:
    static constexpr UINT kMaxNumDescriptor = 8192;

    ComPtr<ID3D12RootSignature> m_RootSignature = nullptr;

    UINT kSrvDescriptorSize{};
    std::vector<Handle> m_FreeList;

    CD3DX12_CPU_DESCRIPTOR_HANDLE baseSrvCpuHandle{};
    ComPtr<ID3D12DescriptorHeap> m_SrvDescriptorHeap = nullptr;
};
