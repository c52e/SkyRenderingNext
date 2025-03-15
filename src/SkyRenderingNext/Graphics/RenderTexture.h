#pragma once

#include "GraphicsDevice.h"
#include "ResourceBarrier.h"

class RenderTexture: public ResourceBarrierBase<RenderTexture>
{
  public:
    RenderTexture() = default;

    RenderTexture(Math::uint2 size, DXGI_FORMAT format, bool allowUav, std::wstring name = L"Unnamed texture");

    RenderTexture(DXGI_FORMAT format, bool allowUav, std::wstring name = L"Unnamed texture");

    RenderTexture(Math::uint2 size, const RenderTexture& fromTex);

    RenderTexture(Math::uint3 size, DXGI_FORMAT format, bool allowUav, std::wstring name = L"Unnamed texture");

    ID3D12Resource *GetResource()
    {
        return m_Resource.Get();
    }

    const D3D12_RESOURCE_DESC &GetDesc() const noexcept
    {
        return m_ResourceDesc;
    }

    const DescriptorHandle &GetUav() const noexcept
    {
        return m_Uav;
    }

    const DescriptorHandle &GetSrv() const noexcept
    {
        return m_Srv;
    }

    UINT GetUavIndex() const noexcept
    {
        return GetUav().index();
    }

    UINT GetSrvIndex() const noexcept
    {
        return GetSrv().index();
    }

  private:
    ResPtr<ID3D12Resource> m_Resource{};
    D3D12_RESOURCE_DESC m_ResourceDesc{};
    DescriptorHandle m_Srv{};
    DescriptorHandle m_Uav{};
    std::wstring m_Name{};
};
